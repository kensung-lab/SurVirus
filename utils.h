#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include "config.h"


enum strand_t {
    F, R
};
struct anchor_t {
    strand_t strand;
    int contig_id, start, end;
    bool dead = false;

    anchor_t() {}

    anchor_t(strand_t strand, int contig_id, int start, int end) : strand(strand), contig_id(contig_id), start(start),
                                                                          end(end) {}

    static bool can_merge(anchor_t& a1, anchor_t& a2, int max_is) {
        if (a1.contig_id != a2.contig_id || a1.strand != a2.strand) return false;
        return std::max(a1.end,  a2.end)-std::min(a1.start, a2.start) <= max_is;
    }

    static anchor_t* merge(anchor_t* a1, anchor_t* a2) {
        return new anchor_t(a1->strand, a1->contig_id, std::min(a1->start, a2->start), std::max(a1->end,  a2->end));
    }

    static int distance(anchor_t& a1, anchor_t& a2) {
        int overlap = std::min(a1.end, a2.end) - std::max(a1.start, a2.start);
        if (overlap > 0) {
            int l2 = std::min(a1.len(), a2.len());
            return -100*overlap/l2;
        } else {
            return std::min(std::abs(a1.start-a2.end), std::abs(a2.start-a1.end));
        }
    }

    int len() {
        return end-start+1;
    }
};
bool operator < (const anchor_t& a1, const anchor_t& a2) {
    if (a1.start != a2.start) return a1.start < a2.start;
    return a1.end < a2.end;
}


struct aa_distance_t {
    int distance;
    anchor_t* a1, * a2;

    aa_distance_t(int distance, anchor_t *a1, anchor_t *a2) : distance(distance), a1(a1), a2(a2) {}
};
bool operator < (const aa_distance_t& aad1, const aa_distance_t& aad2) { // reverse op for priority queue
    return aad1.distance > aad2.distance;
}

struct breakpoint_t {
    std::string chr;
    bool rev;
    int start, end;

    int pos() {return rev ? start : end;}

    breakpoint_t() {}
    breakpoint_t(std::string& str) {
        char chr[1000], strand;
        sscanf(str.data(), "%[^:]:%c:%d:%d", chr, &strand, &start, &end);
        this->chr = chr;
        rev = (strand == '-');
    }

    bool operator == (breakpoint_t& other) {
        return chr == other.chr && rev == other.rev && pos() == other.pos();
    }

    std::string to_string() {
        std::stringstream ssout;
        ssout << chr << ':' << (rev ? '-' : '+') << ":" << start << ":" << end;
        return ssout.str();
    }

    std::string to_human_string() {
        std::stringstream ssout;
        ssout << chr << ':' << (rev ? '-' : '+') << pos();
        return ssout.str();
    }
};

struct call_t {
    int id;
    breakpoint_t host_bp, virus_bp;
    int reads, good_reads, split_reads, reads_w_dups, unique_reads_w_dups, score;
    double host_pbs, virus_pbs;
    double host_cov, virus_cov;
    bool removed = false;
    int paired_with = -1;

    call_t(std::string& line) {
        std::stringstream ssin(line);
        std::string host_bp_str, virus_bp_str;
        ssin >> id >> host_bp_str >> virus_bp_str >> reads >> good_reads >> split_reads >> score >> host_pbs >> virus_pbs
             >> reads_w_dups >> unique_reads_w_dups >> host_cov >> virus_cov;
        host_bp = breakpoint_t(host_bp_str);
        virus_bp = breakpoint_t(virus_bp_str);
    }

    double coverage() { return (host_cov + virus_cov)/2; }

    bool is_paired() { return paired_with >= 0; }

    std::string to_string() {
        std::stringstream ssout;
        ssout << id << " " << host_bp.to_string() << " " << virus_bp.to_string() << " " << reads << " " << good_reads << " ";
        ssout << split_reads << " " << score << " " << host_pbs << " " << virus_pbs << " " << reads_w_dups << " ";
        ssout << unique_reads_w_dups << " " << host_cov << " " << virus_cov;
        return ssout.str();
    }

    std::string to_human_string() {
        char buffer[10000];
        sprintf(buffer, "ID=%d %s %s GOOD_READS=%d TOT_READS=%d SPLIT_READS=%d HOST_PBS=%lf COVERAGE=%lf",
                id, host_bp.to_human_string().c_str(), virus_bp.to_human_string().c_str(), good_reads, reads, split_reads, host_pbs, coverage());
        std::string sout = buffer;
        if (is_paired()) sout += " PAIRED WITH ID=" + std::to_string(paired_with);
        return sout;
    }
};

int pair_dist(call_t& c1, call_t& c2) {
    if (c1.host_bp.chr == c2.host_bp.chr && c1.host_bp.rev != c2.host_bp.rev && c1.virus_bp.rev != c2.virus_bp.rev) {
        int rev_pos, fwd_pos;
        if (c1.host_bp.rev) {
            rev_pos = c1.host_bp.pos();
            fwd_pos = c2.host_bp.pos();
        } else {
            rev_pos = c2.host_bp.pos();
            fwd_pos = c1.host_bp.pos();
        }

        if (rev_pos-fwd_pos >= -50 && rev_pos-fwd_pos <= 1000) {
            return rev_pos-fwd_pos;
        }
    }
    return INT32_MAX;
}

struct chr_seq_t {
    char* seq;
    int len;
    bool circular;

    chr_seq_t(char* seq, int len, bool circular) : seq(seq), len(len), circular(circular) {}
    ~chr_seq_t() {delete[] seq;}
};
typedef std::unordered_map<std::string, chr_seq_t*> chr_seqs_map_t;

void read_fasta_into_map(chr_seqs_map_t& chr_seqs, std::string& reference_fname, bool circular_seq = false) {
    FILE* fasta = fopen(reference_fname.c_str(), "r");
    kseq_t* seq = kseq_init(fileno(fasta));
    while (kseq_read(seq) >= 0) {
        std::string seq_name = seq->name.s;
        if (circular_seq) {
            char* chr_seq = new char[seq->seq.l + seq->seq.l/2 + 1];
            strcpy(chr_seq, seq->seq.s);
            strncpy(chr_seq+seq->seq.l, seq->seq.s, seq->seq.l/2);
            chr_seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l + seq->seq.l/2, true);
        } else {
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            chr_seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l, false);
        }
    }
    kseq_destroy(seq);
    fclose(fasta);
}

#endif //SURVEYOR_CLUSTER_H
