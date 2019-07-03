#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <unistd.h>
#include <htslib/kseq.h>

KSEQ_INIT(int, read)

#include "libs/ssw.h"
#include "libs/ssw_cpp.h"

struct breakpoint_t {
    std::string chr;
    bool rev;
    int pos;

    breakpoint_t() {}
    breakpoint_t(std::string& str) {
        char chr[1000], strand;
        sscanf(str.data(), "%[^:]:%c%d", chr, &strand, &pos);
        this->chr = chr;
        rev = (strand == '-');
    }

    bool operator == (const breakpoint_t& other) {
        return chr == other.chr && rev == other.rev && pos == other.pos;
    }

    std::string to_string() {
        std::stringstream ssout;
        ssout << chr << ':' << (rev ? '-' : '+') << pos;
        return ssout.str();
    }
};

struct call_t {
    int id;
    breakpoint_t host_bp, virus_bp;
    int reads, split_reads, reads_w_dups, unique_reads_w_dups, score;
    double pval_mwu, pval_ks, host_pbs, virus_pbs;
    double host_cov, virus_cov;
    bool removed = false;

    call_t(std::string& line) {
        std::stringstream ssin(line);
        std::string host_bp_str, virus_bp_str;
        std::string id_str;
        ssin >> id_str >> host_bp_str >> virus_bp_str >> reads >> split_reads >> score >> pval_mwu >> pval_ks >> host_pbs >> virus_pbs
            >> reads_w_dups >> unique_reads_w_dups >> host_cov >> virus_cov;
        id = std::stoi(id_str.substr(3));
        host_bp = breakpoint_t(host_bp_str);
        virus_bp = breakpoint_t(virus_bp_str);
    }

    double pval() { return std::max(pval_mwu, pval_ks); }

    double coverage() { return (host_cov + virus_cov)/2; }

    std::string to_string() {
        std::stringstream ssout;
        ssout << "ID=" << id << " " << host_bp.to_string() << " " << virus_bp.to_string() << " READS=" << reads << " SPLIT_READS=" << split_reads;
        ssout << " PVAL=" << pval() << " HOST_PBS=" << host_pbs << " COVERAGE=" << coverage();
        return ssout.str();
    }
};

int pair_dist(call_t& c1, call_t& c2) {
    if (c1.host_bp.chr == c2.host_bp.chr && c1.host_bp.rev != c2.host_bp.rev && c1.virus_bp.rev != c2.virus_bp.rev) {
        int rev_pos, fwd_pos;
        if (c1.host_bp.rev) {
            rev_pos = c1.host_bp.pos;
            fwd_pos = c2.host_bp.pos;
        } else {
            rev_pos = c2.host_bp.pos;
            fwd_pos = c1.host_bp.pos;
        }

        if (rev_pos-fwd_pos >= -50 && rev_pos-fwd_pos <= 1000) {
            return rev_pos-fwd_pos;
        }
    }
    return INT32_MAX;
}

void print_calls(std::vector<call_t>& calls) {
    for (int i = 0; i < calls.size(); i++) {
        call_t call = calls[i];
        if (!call.removed) {
            std::cout << call.to_string() << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {

    std::string workspace = argv[1];

    std::ifstream fin(workspace + "/results.txt");

    std::vector<double> args;
    bool print_rejected = false, accept_clipped = false;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--print-rejected") == 0) {
            print_rejected = true;
        } else if (strcmp(argv[i], "--accept-all-clipped") == 0) {
            accept_clipped = true;
        } else {
            args.push_back(std::stod(argv[i]));
        }
    }

    double min_pval = 0.001;
    if (args.size() >= 1) min_pval = args[0];

    double min_host_pbs = 0.8;
    if (args.size() >= 2) min_host_pbs = args[1];

    double min_cov_for_accept = 0.5;
    if (args.size() >= 3) min_cov_for_accept = args[2];

    std::vector<call_t> accepted_calls, rejected_calls;
    std::string line;
    while (std::getline(fin, line)) {
        call_t call(line);

        bool accept = true;
        if (call.pval() < min_pval) accept = false;
        if (call.host_pbs < min_host_pbs) accept = false;

        if (call.coverage() >= min_cov_for_accept) accept = true;
        if (accept_clipped && call.split_reads > 0) accept = true;

        if (accept) {
            accepted_calls.push_back(call);
        } else {
            rejected_calls.push_back(call);
        }
    }


    // find pairs
    std::vector<bool> paired(accepted_calls.size(), false);

    for (int i = 0; i < accepted_calls.size(); i++) {
        call_t a_call1 = accepted_calls[i];

        int pair_j = -1, min_dist = INT32_MAX;
        for (int j = i+1; j < accepted_calls.size(); j++) {
            call_t a_call2 = accepted_calls[j];
            int dist = pair_dist(a_call1, a_call2);
            if (abs(dist) < min_dist) {
                pair_j = j;
                min_dist = abs(dist);
            }
        }

        if (pair_j != -1) {
            paired[i] = paired[pair_j] = true;
        }
    }

    for (int i = 0; i < accepted_calls.size(); i++) {
        if (paired[i]) continue;

        call_t a_call = accepted_calls[i];

        int pair_j = -1, min_dist = INT32_MAX;
        for (int j = 0; j < rejected_calls.size(); j++) {
            call_t r_call = rejected_calls[j];
            int dist = pair_dist(a_call, r_call);
            if (abs(dist) < min_dist) {
                pair_j = j;
                min_dist = abs(dist);
            }
        }

        if (pair_j != -1) {
            accepted_calls.push_back(rejected_calls[pair_j]);
            paired[i] = true;
            paired.push_back(true);
        }
    }


    // find duplicates
    std::string host_bp_seqs_fname = workspace + "/host_bp_seqs.fa", virus_bp_seqs_fname = workspace + "/virus_bp_seqs.fa";
    FILE* host_bp_fasta = fopen(host_bp_seqs_fname.c_str(), "r"), *virus_bp_fasta = fopen(virus_bp_seqs_fname.c_str(), "r");
    kseq_t* hseq = kseq_init(fileno(host_bp_fasta)), *vseq = kseq_init(fileno(virus_bp_fasta));
    std::vector<std::pair<std::string, std::string> > host_virus_bp_seqs;
    while (kseq_read(hseq) >= 0 && kseq_read(vseq) >= 0) {
        host_virus_bp_seqs.push_back({hseq->seq.s, vseq->seq.s});
    }

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;
    for (int j = 0; j < accepted_calls.size(); j++) {
        std::string hseq_j = host_virus_bp_seqs[accepted_calls[j].id].first;
        std::string vseq_j = host_virus_bp_seqs[accepted_calls[j].id].second;

        for (int i = 0; i < j; i++) {
            std::string hseq_i = host_virus_bp_seqs[accepted_calls[i].id].first;
            std::string vseq_i = host_virus_bp_seqs[accepted_calls[i].id].second;

            std::string& h_ref = hseq_i.length() >= hseq_j.length() ? hseq_i : hseq_j;
            std::string& h_query = hseq_i.length() >= hseq_j.length() ? hseq_j : hseq_i;
            std::string& v_ref = vseq_i.length() >= vseq_j.length() ? vseq_i : vseq_j;
            std::string& v_query = vseq_i.length() >= vseq_j.length() ? vseq_j : vseq_i;

            StripedSmithWaterman::Alignment h_aln, v_aln;
            aligner.Align(h_query.c_str(), h_ref.c_str(), h_ref.length(), filter, &h_aln, 0);
            bool same_hseq = h_aln.sw_score >= h_query.length()*0.8;

            aligner.Align(v_query.c_str(), v_ref.c_str(), v_ref.length(), filter, &v_aln, 0);
            bool same_vseq = v_aln.sw_score >= v_query.length()*0.8;

            if (same_hseq && same_vseq) {
                if (hseq_j.length() > hseq_i.length()) {
                    accepted_calls[i].host_bp = accepted_calls[j].host_bp;
                    host_virus_bp_seqs[i].first = host_virus_bp_seqs[j].first;
                }
                if (vseq_j.length() > vseq_i.length()) {
                    accepted_calls[i].virus_bp = accepted_calls[j].virus_bp;
                    host_virus_bp_seqs[i].second = host_virus_bp_seqs[j].second;
                }
                accepted_calls[j].removed = true;
                break;
            }
        }
    }


    if (print_rejected) {
        print_calls(rejected_calls);
    } else {
        print_calls(accepted_calls);
    }
}