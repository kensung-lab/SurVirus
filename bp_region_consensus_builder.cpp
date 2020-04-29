#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <htslib/sam.h>

#include "sam_utils.h"
#include "utils.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"


std::vector<bam1_t*> get_redux_reads(char* reads_fname, std::string& target_chr) {
    std::vector<bam1_t*> reads;
    bam1_t* read = bam_init1();
    open_samFile_t* bam_file = open_samFile(reads_fname, true);
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, target_chr.c_str());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        reads.push_back(bam_dup1(read));
    }
    close_samFile(bam_file);
    bam_destroy1(read);
    return reads;
}

std::string build_consensus(char* ref_region, int start, std::vector<bam1_t*>& reads) {
    int ref_len = strlen(ref_region);
    for (int i = 0; i < ref_len; i++) {
        ref_region[i] = toupper(ref_region[i]);
    }

    auto comp = [](const bam1_t* r1, const bam1_t* r2) {
        if (r1->core.pos == r2->core.pos) {
            return bam_endpos(r1) < bam_endpos(r2);
        }
        return r1->core.pos < r2->core.pos;
    };
    std::sort(reads.begin(), reads.end(), comp);

    std::vector<std::string> read_sequences;
    std::vector<int> read_starts, read_ends, read_pos;
    for (bam1_t* read : reads) {
        std::string temp_sequence = get_sequence(read);
        char read_sequence[10000];
        int pos_dst = 0, pos_src = 0;
        uint32_t* cigar = bam_get_cigar(read);
        for (int i = 0; i < read->core.n_cigar; i++) {
            char op = bam_cigar_opchr(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (op == 'M') {
                strncpy(read_sequence + pos_dst, temp_sequence.c_str() + pos_src, len);
                pos_dst += len;
                pos_src += len;
            }  else if (op == 'D') {
                int len = bam_cigar_oplen(cigar[i]);
                memset(read_sequence+pos_dst, 'D', len);
                pos_dst += len;
            } else if (op == 'I') {
                int len = bam_cigar_oplen(cigar[i]);
                for (int j = 0; j < len; j++) {
                    read_sequence[pos_dst++] = tolower(temp_sequence[pos_src++]);
                }
            } else if (op == 'S') {
                pos_src += bam_cigar_oplen(cigar[i]);
            }
        }
        read_sequence[pos_dst] = '\0';
        read_sequences.push_back(read_sequence);
        read_starts.push_back(read->core.pos-start);
        read_ends.push_back(bam_endpos(read)-start);
        read_pos.push_back(0);
    }

    uint8_t base_to_idx[255] = {};
    base_to_idx['A'] = 1; base_to_idx['C'] = 2; base_to_idx['G'] = 3; base_to_idx['T'] = 4;
    base_to_idx['a'] = 5; base_to_idx['c'] = 6; base_to_idx['g'] = 7; base_to_idx['t'] = 8;
    base_to_idx['D'] = 9;

    char idx_to_base[10] = {'N', 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'D'};

    char consensus[10000];
    int consensus_i = 0;

    int s = 0, e = 1;
    int n_reads = reads.size();
    for (int i = 0; i < ref_len; i++) {
        while (s < n_reads && read_ends[s] < i) s++;
        while (e < n_reads && read_starts[e] <= i) e++;

        int freq[10] = {};
        for (int j = s; j < e; j++) {
            if (read_ends[j] < i) continue;

            char base = read_sequences[j][read_pos[j]];
            freq[base_to_idx[base]]++;
        }

        int max_idx = 0, n_maxs = 0;
        for (int j = 1; j < 10; j++) {
            if (freq[max_idx] < freq[j]) {
                max_idx = j;
            }
        }
        if (freq[max_idx] >= 2) { // we want at least 2 reads for a consensus, otherwise we trust the reference
            for (int j = 1; j < 10; j++) {
                n_maxs += (freq[max_idx] == freq[j]);
            }
        }

        char opt_base = ref_region[i];
        if (n_maxs == 1) { // if undecided on consensus, just choose reference
            opt_base = idx_to_base[max_idx];
        }

        if (!isupper(opt_base)) {
            i--;
        }
        for (int j = s; j < e; j++) {
            if (read_ends[j] < i) continue;

            char base = read_sequences[j][read_pos[j]];
            read_pos[j] += (isupper(base) == isupper(opt_base));
        }

        if (opt_base != 'D') {
            consensus[consensus_i++] = opt_base;
        }
    }
    consensus[consensus_i] = '\0';
    return consensus;
}

bool is_good_alignment(StripedSmithWaterman::Alignment& alignment, int query_len) {
    int query_aln_len = alignment.query_end-alignment.query_begin;
    return query_aln_len >= query_len*0.75; // at least 75% aligned
}

int main(int argc, char* argv[]) {

    std::string host_reference_fname  = argv[1];
    std::string virus_reference_fname  = argv[2];
    std::string workdir = argv[3];

    config_t config = parse_config(workdir + "/config.txt");
    config.min_sc_size = 10;

    // find global max_is
    int max_is = 0;
    for (int i = 4; i < argc; i++) {
        std::string workspace = argv[i];
        stats_t stats = parse_stats(workspace + "/stats.txt");
        max_is = std::max(max_is, stats.max_is);
    }

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(host_reference_fname, false, false);
    chr_seqs.read_fasta_into_map(virus_reference_fname, true, true);

    std::ofstream host_bp_seqs(workdir + "/host_bp_seqs.fa");
    std::ofstream virus_bp_seqs(workdir + "/virus_bp_seqs.fa");

    std::string res_fname = workdir + "/results.txt";
    std::ifstream res_fin(res_fname);

    std::string res_remapped_fname = workdir + "/results.remapped.txt";
    std::ofstream res_remapped_fout(res_remapped_fname);

    std::string line;
    while (std::getline(res_fin, line)) {
        call_t call(line);

        char reads_fname[1024];
        sprintf(reads_fname, "%s/readsx/%d.bam", workdir.c_str(), call.id);
        std::vector<bam1_t*> host_reads = get_redux_reads(reads_fname, call.host_bp.chr);
        std::vector<bam1_t*> virus_reads = get_redux_reads(reads_fname, call.virus_bp.chr);

        host_bp_seqs << ">" << call.id << "_" << (call.host_bp.rev ? "L" : "R") << "\n";
        std::string consensus_host = "A";
        if (!host_reads.empty()) {
            char host_reg[10000];
            strncpy(host_reg, chr_seqs.get_seq(call.host_bp.chr)+call.host_bp.start, call.host_bp.end-call.host_bp.start);
            host_reg[call.host_bp.end-call.host_bp.start] = '\0';
            consensus_host = build_consensus(host_reg, call.host_bp.start, host_reads);
        }
        host_bp_seqs << consensus_host << std::endl;

        virus_bp_seqs << ">" << call.id << "_" << (call.virus_bp.rev ? "L" : "R") << "\n";
        std::string consensus_virus = "A";
        if (!virus_reads.empty()) {
            // if circular call, move reads to second half
            int virus_len = chr_seqs.get_original_len(call.virus_bp.chr);
            if ((!call.virus_bp.rev && call.virus_bp.pos() < max_is) || (call.virus_bp.rev && call.virus_bp.pos() > virus_len-max_is)) {
                for (bam1_t* read : virus_reads) {
                    if (bam_endpos(read) <= max_is) {
                        read->core.pos += virus_len;
                    }
                }
            }

            int virus_start = 1000000000, virus_end = 0;
            for (bam1_t* read : virus_reads) {
                virus_start = std::min(virus_start, read->core.pos);
                virus_end = std::max(virus_end, bam_endpos(read));
            }

            char virus_reg[10000];
            strncpy(virus_reg, chr_seqs.get_seq(call.virus_bp.chr)+virus_start, virus_end-virus_start);
            virus_reg[virus_end-virus_start] = '\0';
            consensus_virus = build_consensus(virus_reg, call.virus_bp.start, virus_reads);
        }
        virus_bp_seqs << consensus_virus << std::endl;



        /* == remap reads == */

        StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
        StripedSmithWaterman::Filter filter;

        if (call.host_bp.rev) get_rc(consensus_host);
        if (!call.virus_bp.rev) get_rc(consensus_virus);
        int junction = consensus_host.length();
        std::string full_region = consensus_host + consensus_virus;

        std::unordered_map<std::string, std::pair<int, int> > pairs_ref_ends;
        for (bam1_t* read : host_reads) { // initialise pairs_ref_ends
            pairs_ref_ends[bam_get_qname(read)] = {full_region.length(), 0};
        }

        uint64_t tot_score = 0, tot_bases = 0;

        // map host reads
        std::unordered_set<std::string> good_host_reads;
        std::vector<std::string> splitreads;
        for (bam1_t* read : host_reads) {
            if (!is_primary(read)) continue;
            if (is_low_complexity(read, false, get_left_clip_len(read), read->core.l_qseq-get_right_clip_len(read)))
                continue;

            StripedSmithWaterman::Alignment alignment;
            std::string read_seq = get_sequence(read);
            if (call.host_bp.rev) get_rc(read_seq);
            aligner.Align(read_seq.c_str(), full_region.c_str(), full_region.length(), filter, &alignment, 0);

            std::string qname = bam_get_qname(read);
            std::pair<int, int>& pre = pairs_ref_ends[qname];
            pre.first = std::min(pre.first, alignment.ref_begin);
            pre.second = std::max(pre.second, alignment.ref_end);

            if (is_good_alignment(alignment, read_seq.length())) {
                good_host_reads.insert(qname);
                tot_score += alignment.sw_score;
                tot_bases += std::max(alignment.ref_end-alignment.ref_begin, alignment.query_end-alignment.query_begin);

                if (alignment.ref_begin <= junction-config.min_sc_size && junction+config.min_sc_size <= alignment.ref_end) {
                    splitreads.push_back(qname);
                }
            }
        }

        // map virus reads
        std::unordered_set<std::string> good_pairs;
        for (bam1_t* read : virus_reads) {
            if (!is_primary(read)) continue;

            if (is_low_complexity(read, false, get_left_clip_len(read), read->core.l_qseq-get_right_clip_len(read)))
                continue;
            if (good_host_reads.count(bam_get_qname(read)) == 0) continue;

            StripedSmithWaterman::Alignment alignment;
            std::string read_seq = get_sequence(read);
            if (!call.virus_bp.rev) get_rc(read_seq);
            aligner.Align(read_seq.c_str(), full_region.c_str(), full_region.length(), filter, &alignment, 0);

            std::string qname = bam_get_qname(read);
            std::pair<int, int>& pre = pairs_ref_ends[qname];
            pre.first = std::min(pre.first, alignment.ref_begin);
            pre.second = std::max(pre.second, alignment.ref_end);

            if (is_good_alignment(alignment, read_seq.length())) {
                std::string qname = bam_get_qname(read);
                if (pre.second-pre.first <= max_is) {
                    good_pairs.insert(qname);
                }
                tot_score += alignment.sw_score;
                tot_bases += std::max(alignment.ref_end-alignment.ref_begin, alignment.query_end-alignment.query_begin);

                if (alignment.ref_begin <= junction-config.min_sc_size && junction+config.min_sc_size <= alignment.ref_end) {
                    splitreads.push_back(qname);
                }
            }
        }

        call.good_pairs = good_pairs.size();
        call.split_reads = splitreads.size();
        call.host_pbs = call.virus_pbs = tot_score/double(tot_bases);

        res_remapped_fout << call.to_string() << std::endl;
    }

    host_bp_seqs.close();
    virus_bp_seqs.close();

    res_remapped_fout.close();
}
