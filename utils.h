#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include <htslib/kseq.h>
#include <unistd.h>
#include <vector>
#include <map>

KSEQ_INIT(int, read)

#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"


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
    if (a1.contig_id != a2.contig_id) return a1.contig_id < a2.contig_id;
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
    breakpoint_t(std::string chr, int start, int end, bool rev) : chr(chr), start(start), end(end), rev(rev) {}

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
    int reads, good_pairs, split_reads, reads_w_dups, unique_reads_w_dups, score;
    double host_pbs, virus_pbs;
    double host_cov, virus_cov;
    bool removed = false;
    int paired_with = -1;

    call_t(std::string& line) {
        std::stringstream ssin(line);
        std::string host_bp_str, virus_bp_str;
        ssin >> id >> host_bp_str >> virus_bp_str >> reads >> good_pairs >> split_reads >> score >> host_pbs >> virus_pbs
             >> reads_w_dups >> unique_reads_w_dups >> host_cov >> virus_cov;
        host_bp = breakpoint_t(host_bp_str);
        virus_bp = breakpoint_t(virus_bp_str);
    }

    call_t(int id, breakpoint_t& host_bp, breakpoint_t& virus_bp, int reads, int good_pairs, int split_reads,
        int reads_w_dups, int uniquer_reads_w_dups, int score, double host_pbs, double virus_pbs,
        double host_cov, double virus_cov) :
    id(id), host_bp(host_bp), virus_bp(virus_bp), reads(reads), good_pairs(good_pairs), split_reads(split_reads),
    reads_w_dups(reads_w_dups), unique_reads_w_dups(uniquer_reads_w_dups), score(score),
    host_pbs(host_pbs), virus_pbs(virus_pbs),
    host_cov(host_cov), virus_cov(virus_cov) {}

    double coverage() { return (host_cov + virus_cov)/2; }

    bool is_paired() { return paired_with >= 0; }

    std::string to_string() {
        char buffer[4096];
        sprintf(buffer, "%d %s %s %d %d %d %d %lf %lf %d %d %lf %lf", id,
                host_bp.to_string().c_str(), virus_bp.to_string().c_str(), reads, good_pairs, split_reads, score,
                host_pbs, virus_pbs, reads_w_dups, unique_reads_w_dups, host_cov, virus_cov);
        return buffer;
    }

    std::string to_human_string() {
        char buffer[10000];
        sprintf(buffer, "ID=%d %s %s SUPPORTING_PAIRS=%d SPLIT_READS=%d HOST_PBS=%lf COVERAGE=%lf",
                id, host_bp.to_human_string().c_str(), virus_bp.to_human_string().c_str(), good_pairs, split_reads, host_pbs, coverage());
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

struct region_t {
    int id;
    int contig_id; // id in our own mapping
    int start, end;

    region_t(int contig_id, int start, int end) : id(-1), contig_id(contig_id), start(start), end(end) {}

    int len() {
        return end-start+1;
    }

    bool overlaps_with(region_t* other) {
        if (contig_id != other->contig_id) return false;
        return std::max(start, other->start) < std::min(end, other->end);
    }
};


struct chr_seq_t {
    char* seq;
    int len;
    bool circular, is_virus;

    chr_seq_t(char* seq, int len, bool is_virus, bool circular) : seq(seq), len(len), is_virus(is_virus), circular(circular) {}
    ~chr_seq_t() {delete[] seq;}
};
struct chr_seqs_map_t {
    std::unordered_map<std::string, chr_seq_t*> seqs;

    void read_fasta_into_map(std::string& reference_fname, bool is_virus = false, bool circular_seq = false) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            if (circular_seq) {
                char* chr_seq = new char[seq->seq.l + seq->seq.l/2 + 1];
                strcpy(chr_seq, seq->seq.s);
                strncpy(chr_seq+seq->seq.l, seq->seq.s, seq->seq.l/2);
                chr_seq[seq->seq.l + seq->seq.l/2] = '\0';
                seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l + seq->seq.l/2, true, is_virus);
            } else {
                char* chr_seq = new char[seq->seq.l + 1];
                strcpy(chr_seq, seq->seq.s);
                chr_seq[seq->seq.l] = '\0';
                seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l, false, is_virus);
            }
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        return seqs[seq_name]->seq;
    }

    int get_len(std::string seq_name) {
        return seqs[seq_name]->len;
    }
    int get_original_len(std::string seq_name) {
        auto& e = seqs[seq_name];
        if (e->circular) return e->len - e->len/3;
        return e->len;
    }

    bool is_virus(std::string seq_name) {
        return seqs[seq_name]->is_virus;
    }

    ~chr_seqs_map_t() {
        for (auto& e : seqs) {
            delete e.second;
        }
    }
};


struct contig_map_t {
    std::vector<std::string> contig_id2name;
    std::unordered_map<std::string, int> contig_name2id;

    contig_map_t() {}
    contig_map_t(std::string& workdir) {
        std::ifstream contig_map_fin(workdir + "/contig_map");
        std::string contig_name; int contig_id;
        contig_id2name.push_back("");
        while (contig_map_fin >> contig_name >> contig_id) {
            contig_id2name.push_back(contig_name);
            contig_name2id[contig_name] = contig_id;

        }
    }

    std::string id2name(int id) { return contig_id2name[id]; }
    int name2id(std::string name) { return contig_name2id[name]; }
};



void remove_anchor_from_mm(std::multimap<int, anchor_t *> &mm, anchor_t *a, int pos) {
    auto bounds = mm.equal_range(pos);
    for (auto it = bounds.first; it != bounds.second; it++) {
        if (it->second == a) {
            mm.erase(it);
            break;
        }
    }
}
void remove_anchor_from_mm(std::multimap<int, anchor_t *> &mm, anchor_t *a) {
    remove_anchor_from_mm(mm, a, a->start);
    remove_anchor_from_mm(mm, a, a->end);
}
void extract_regions(std::vector<bam1_t *> &reads, std::vector<region_t *> &regions, int max_is, contig_map_t& contig_map,
                     std::unordered_map<int, int>& contig_tid2id, chr_seqs_map_t& chr_seqs, bool process_xa = true) {

    if (reads.empty()) return;

    std::vector<anchor_t*> anchors;
    for (bam1_t* read : reads) {
        anchor_t* a = new anchor_t(bam_is_rev(read) ? strand_t::R : strand_t::F, contig_tid2id[read->core.tid],
                                   read->core.pos, bam_endpos(read));
        anchors.push_back(a);

        uint8_t* xa = bam_aux_get(read, "XA");
        if (xa != NULL && process_xa) {
            std::string xa_s = bam_aux2Z(xa);
            size_t pos = 0, prev = 0;
            while ((pos = xa_s.find(';', pos+1)) != std::string::npos) { //TODO: crashes if virus name contains ";"
                std::string alt_align = xa_s.substr(prev, pos-prev);
                size_t pos2 = 0;

                std::string xa_tname = alt_align.substr(0, (pos2 = alt_align.find(',')));

                std::string xa_dir_pos = alt_align.substr(pos2+1, alt_align.find(',', pos2+1)-pos2-1);
                pos2 = alt_align.find(',', pos2+1);
                char xa_dir = xa_dir_pos[0];
                int xa_pos = std::stoi(xa_dir_pos.substr(1));

                std::string xa_cigar = alt_align.substr(pos2+1, alt_align.find(',', pos2+1)-pos2-1);
                std::pair<int, const uint32_t*> cigar_array = cigar_str_to_array(xa_cigar);

                int qlen = bam_cigar2qlen(cigar_array.first, cigar_array.second);
                delete[] cigar_array.second;

                prev = pos+1;

                anchor_t* a = new anchor_t(xa_dir == '-' ? strand_t::R : strand_t::F, contig_map.name2id(xa_tname), xa_pos, xa_pos+qlen);
                anchors.push_back(a);
            };
        }
    }

    sort(anchors.begin(), anchors.end(), [](const anchor_t* a1, const anchor_t* a2) {
        return *a1 < *a2;
    });


    // merge first similar anchors (essentially downsample reads)
    int prev;
    do {
        prev = anchors.size();
        for (int i = 0; i < anchors.size()-1; i++) {
            anchor_t* a1 = anchors[i];
            anchor_t* a2 = anchors[i+1];
            if (a1 != nullptr && std::abs(a1->start-a2->start) <= 10 && std::abs(a1->end-a2->end) <= 10) {
                anchors[i] = anchor_t::merge(a1, a2);
                anchors[i+1] = nullptr;
                delete a1;
                delete a2;
            }
        }
        anchors.erase(std::remove(anchors.begin(), anchors.end(), nullptr), anchors.end());
    } while (prev != anchors.size());


    std::vector<int> max_dists;
    for (int i = 0; i < 10; i++) max_dists.push_back(i);
    for (int i = 10; i < 100; i+=10) max_dists.push_back(i);
    max_dists.push_back(max_is);

    for (int max_dist : max_dists) {
        std::multimap<int, anchor_t*> anchors_map;
        for (anchor_t* a : anchors) {
            anchors_map.insert({a->start, a});
            anchors_map.insert({a->end, a});
        }

        std::priority_queue<aa_distance_t> pq;
        for (anchor_t* a1 : anchors) {
            if (a1->dead) continue;
            auto end = anchors_map.upper_bound(a1->end+max_dist);
            for (auto map_it = anchors_map.lower_bound(a1->start); map_it != end; map_it++) {
                anchor_t* a2 = map_it->second;
                if (a1 != a2 && anchor_t::can_merge(*a1, *a2, max_is) &&
                    (a1->start <= a2->start)) {
                    pq.push(aa_distance_t(anchor_t::distance(*a1, *a2), a1, a2));
                }
            }
        }

        while (!pq.empty()) {
            aa_distance_t aad = pq.top();
            pq.pop();

            if (aad.a1->dead || aad.a2->dead) continue;

            anchor_t* new_anchor = anchor_t::merge(aad.a1, aad.a2);
            anchors.push_back(new_anchor);

            aad.a1->dead = true;
            remove_anchor_from_mm(anchors_map, aad.a1);
            aad.a2->dead = true;
            remove_anchor_from_mm(anchors_map, aad.a2);

            auto end = anchors_map.upper_bound(new_anchor->end + max_dist);
            for (auto map_it = anchors_map.lower_bound(new_anchor->start - max_dist); map_it != end; map_it++) {
                if (anchor_t::can_merge(*new_anchor, *map_it->second, max_is)) {
                    pq.push(aa_distance_t(anchor_t::distance(*new_anchor, *map_it->second), new_anchor,
                                          map_it->second));
                }
            }
            anchors_map.insert({new_anchor->start, new_anchor});
            anchors_map.insert({new_anchor->end, new_anchor});
        }

        for (int i = 0; i < anchors.size(); i++) {
            if (anchors[i]->dead) {
                delete anchors[i];
                anchors[i] = nullptr;
            }
        }
        anchors.erase(std::remove(anchors.begin(), anchors.end(), nullptr), anchors.end());
    }

    sort(anchors.begin(), anchors.end(), [](const anchor_t* a1, const anchor_t* a2) {
        return *a1 < *a2;
    });


    // remove anchors entirely contained in another anchor
    do {
        prev = anchors.size();
        for (int i = 0; i < anchors.size()-1; i++) {
            anchor_t* a1 = anchors[i];
            anchor_t* a2 = anchors[i+1];
            if (a1 != nullptr && a1->contig_id == a2->contig_id && a1->end+10 > a2->end) {
                anchors[i] = anchor_t::merge(a1, a2);
                anchors[i+1] = nullptr;
                delete a1;
                delete a2;
            }
        }
        anchors.erase(std::remove(anchors.begin(), anchors.end(), nullptr), anchors.end());
    } while (prev != anchors.size());

    for (anchor_t* a : anchors) {
        if (!a->dead) {
            std::string contig_name = contig_map.id2name(a->contig_id);
            int len = (int) chr_seqs.get_len(contig_name);
            int start = std::max(0, a->start - max_is + a->len());
            int end = std::min(a->end + max_is - a->len(), len);
            region_t* region = new region_t(a->contig_id, start, end);
            regions.push_back(region);
        };
        delete a;
    }
}


struct read_realignment_t {
    bam1_t* read;
    uint16_t offset_start, offset_end;
    uint8_t cigar_len;
    const uint32_t* cigar;
    uint16_t score;
    bool rc;
    bool suspicious = false, too_long = false;

    read_realignment_t() : read(NULL), offset_start(0), offset_end(0), cigar_len(0), cigar(NULL), score(0), rc(false) {}
    read_realignment_t(bam1_t* read) : read(read), offset_start(0), offset_end(0), cigar_len(0), cigar(NULL), score(0), rc(false) {};
    read_realignment_t(bam1_t* read, uint16_t offset_start, uint16_t offset_end, uint8_t cigar_len, const uint32_t* cigar,
                       uint16_t score, bool rc)
            : read(read), offset_start(offset_start), offset_end(offset_end), cigar_len(cigar_len), cigar(cigar), score(score), rc(rc) {}

    bool accepted() { return cigar != NULL; }

    int left_clip_len() {
        return bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]) : 0;
    }
    int left_clip_avg_qual() {
        // if bam_is_rev, get_sequence(read, true) returns a sequence that is RC wrt the BAM file
        // if read_realignment.rc, then we RCed the result of get_sequence(read, true) before aligning to the region
        // this means, if bam_is_rev == read_realignment.rc, the sequence is the same as the BAM file
        // otherwise, we RCed it, and the qualities are now reversed
        return bam_is_rev(read) != rc ? get_avg_qual(read, read->core.l_qseq-left_clip_len(), read->core.l_qseq) :
               get_avg_qual(read, 0, left_clip_len());
    }

    int right_clip_len() {
        return bam_cigar_opchr(cigar[cigar_len-1]) == 'S' ? bam_cigar_oplen(cigar[cigar_len-1]) : 0;
    }
    int right_clip_avg_qual() {
        return bam_is_rev(read) != rc ? get_avg_qual(read, 0, right_clip_len()) :
               get_avg_qual(read, read->core.l_qseq-right_clip_len(), read->core.l_qseq);
    }
};

char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
}

std::string alignment_cigar_to_bam_cigar(std::vector<uint32_t> cigar) {
    std::stringstream ss;
    char op = ' '; int len = 0;
    for (uint32_t c : cigar) {
        if (op != _cigar_int_to_op(c)) {
            if (op != ' ') ss << len << op;
            op = _cigar_int_to_op(c);
            len = cigar_int_to_len(c);
        } else {
            len += cigar_int_to_len(c);
        }
    }
    ss << len << op;
    return ss.str();
}

bool accept_alignment(StripedSmithWaterman::Alignment& alignment, std::string& query, int min_sc_size) {
    bool long_enough = alignment.ref_end-alignment.ref_begin+1 >= 30;
    bool score_enough = alignment.sw_score >= 30;
    uint32_t c0 = alignment.cigar[0], cl = alignment.cigar[alignment.cigar.size()-1];
    bool left_clipped = cigar_int_to_op(c0) == 'S' && cigar_int_to_len(c0) >= min_sc_size;
    bool right_clipped = cigar_int_to_op(cl) == 'S' && cigar_int_to_len(cl) >= min_sc_size;
    bool clipped_both_sides = left_clipped && right_clipped;
    return long_enough && score_enough && !clipped_both_sides;
}


void add_realignment_to_region(read_realignment_t& rr, int region_id, int min_sc_size,
                               std::vector<read_realignment_t>* reads_per_region, std::vector<read_realignment_t>* reads_per_region_rc) {
    if (rr.read->core.flag & BAM_FSECONDARY) { // good clip
        if (rr.left_clip_len() >= min_sc_size) {
            reads_per_region_rc[region_id].push_back(rr);
        }
        if (rr.right_clip_len() >= min_sc_size) {
            reads_per_region[region_id].push_back(rr);
        }
    } else { // the read points toward the bp
        std::vector<read_realignment_t>& reads_per_region_v = rr.rc ? reads_per_region_rc[region_id] : reads_per_region[region_id];
        reads_per_region_v.push_back(rr);
    }
}

#endif //SURVEYOR_CLUSTER_H
