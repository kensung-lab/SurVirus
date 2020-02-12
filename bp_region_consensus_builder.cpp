#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <cassert>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "utils.h"

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;

std::unordered_map<std::string, int> virus_names;


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

int main(int argc, char* argv[]) {

    std::string host_reference_fname  = argv[1];
    std::string virus_reference_fname  = argv[2];
    std::string workdir = argv[3];

    // find global max_is
    int max_is = 0;
    for (int i = 4; i < argc; i++) {
        std::string workspace = argv[i];
        stats_t stats = parse_stats(workspace + "/stats.txt");
        max_is = std::max(max_is, stats.max_is);
    }

    chr_seqs_map_t chr_seqs;
    read_fasta_into_map(chr_seqs, host_reference_fname);
    read_fasta_into_map(chr_seqs, virus_reference_fname, true);

    std::ofstream host_bp_seqs(workdir + "/host_bp_seqs2.fa");
    std::ofstream virus_bp_seqs(workdir + "/virus_bp_seqs2.fa");

    std::string res_fname = workdir + "/results.txt";
    std::ifstream res_fin(res_fname);
    std::string line;
    while (std::getline(res_fin, line)) {
        call_t call(line);

        char reads_fname[1024];
        sprintf(reads_fname, "%s/readsx/%d.sorted.bam", workdir.c_str(), call.id);
        std::vector<bam1_t*> host_reads = get_redux_reads(reads_fname, call.host_bp.chr);
        std::vector<bam1_t*> virus_reads = get_redux_reads(reads_fname, call.virus_bp.chr);

        char host_reg[10000];
        strncpy(host_reg, chr_seqs[call.host_bp.chr]->seq+call.host_bp.start, call.host_bp.end-call.host_bp.start);
        host_reg[call.host_bp.end-call.host_bp.start] = '\0';

        std::string consensus_host = build_consensus(host_reg, call.host_bp.start, host_reads);
        host_bp_seqs << ">" << call.id << "\n";
        host_bp_seqs << consensus_host << std::endl;

        // if circular call, move reads to second half
        int virus_len = chr_seqs[call.virus_bp.chr]->len;
        if (chr_seqs[call.virus_bp.chr]->circular) virus_len = 2*virus_len/3;
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
        strncpy(virus_reg, chr_seqs[call.virus_bp.chr]->seq+virus_start, virus_end-virus_start);
        virus_reg[virus_end-virus_start] = '\0';

        std::string consensus_virus = build_consensus(virus_reg, call.virus_bp.start, virus_reads);
        virus_bp_seqs << ">" << call.id << "\n";
        virus_bp_seqs << consensus_virus << std::endl;
    }

    host_bp_seqs.close();
    virus_bp_seqs.close();
}
