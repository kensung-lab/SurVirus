#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <queue>
#include <cassert>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <htslib/hts_endian.h>
#include <cstring>
#include <algorithm>
#include <random>
#include <chrono>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "cluster.h"
#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

typedef unsigned long long ull;

int KMER_LEN = 13;
int KMER_BITS = KMER_LEN * 2;
ull KMER_MASK = (1ll << KMER_BITS)-1;

int POP_SIZE = 10000;

config_t config;
std::vector<stats_t> bam_stats;
stats_t global_stats;
int host_maxIS = 0;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;
std::unordered_map<int, int> contig_id2tid, contig_tid2id;

std::unordered_map<std::string, std::pair<char*, size_t> > chrs;
std::unordered_map<std::string, int> virus_names;

google::dense_hash_map<std::string, std::pair<int, const uint32_t*> > cstring_to_cigar;
std::vector<std::pair<int, const uint32_t*> > cid_to_cigar;

google::dense_hash_map<std::string, std::vector<std::string> > dedup_pairs_by_head_qname;

ull nucl_bm[256] = { 0 };
char bm_nucl[4] = { 'A', 'C', 'G', 'T' };
char nucl2chr[16];


std::atomic<int> virus_integration_id(0);

struct int_breakpoint_t {
    int contig_id;
    int min_pos, max_pos;
    bool fwd;

    int_breakpoint_t(int contig_id, int min_pos, int max_pos, bool fwd) : contig_id(contig_id), min_pos(min_pos),
                                                                          max_pos(max_pos), fwd(fwd) {}

    int bp() {
        return fwd ? max_pos : min_pos;
    }

    std::string to_str() {
        char buffer[4096];
        int start = min_pos, end = max_pos;
        if (virus_names.count(contig_id2name[contig_id])) {
            start %= virus_names[contig_id2name[contig_id]];
            end %= virus_names[contig_id2name[contig_id]];
        }
        sprintf(buffer, "%s:%c:%d:%d", contig_id2name[contig_id].c_str(), fwd ? '+' : '-', start, end);
        return buffer;
    }
};
struct virus_integration_t {

    int id;
    int_breakpoint_t host_bp, virus_bp;
    int reads, good_reads, split_reads, reads_w_dups, unique_reads_w_dups;
    int score;
    double h_score_ratio, v_score_ratio;
    double host_coverage, virus_coverage;

    virus_integration_t(int contig_id, int_breakpoint_t& host_bp, int_breakpoint_t& virus_bp, int reads, int good_reads, int split_reads,
            int reads_w_dups, int uniquer_reads_w_dups, int score, double h_score_ratio,
            double v_score_ratio, double host_coverage, double virus_coverage) :
    id(virus_integration_id++), host_bp(host_bp), virus_bp(virus_bp), reads(reads), good_reads(good_reads), split_reads(split_reads),
            reads_w_dups(reads_w_dups), unique_reads_w_dups(uniquer_reads_w_dups), score(score),
            h_score_ratio(h_score_ratio), v_score_ratio(v_score_ratio),
            host_coverage(host_coverage), virus_coverage(virus_coverage) {}

    std::string to_str() {
        char buffer[4096];
        sprintf(buffer, "%d %s %s %d %d %d %d %lf %lf %d %d %lf %lf", id,
                host_bp.to_str().c_str(), virus_bp.to_str().c_str(), reads, good_reads, split_reads, score,
                h_score_ratio, v_score_ratio, reads_w_dups, unique_reads_w_dups, host_coverage, virus_coverage);
        return buffer;
    }
};


struct region_t {
    int id;
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : id(-1), contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}

    std::string to_str() {
        return contig_id2name[contig_id] + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    int len() {
        return end-start+1;
    }

    bool overlaps_with(region_t* other) {
        if (contig_id != other->contig_id) return false;
        return std::max(start, other->start) < std::min(end, other->end);
    }
};

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

struct region_score_t {
    int id;
    region_t* region;
    bool fwd;
    int score, reads, good_reads;

    region_score_t(int id, region_t* region, bool fwd, int score, int reads, int good_reads)
    : id(id), region(region), fwd(fwd), score(score), reads(reads), good_reads(good_reads) {};

    std::string to_str() {
        return std::to_string(id) + " " + region->to_str() + " READS=" + std::to_string(reads) + " SCORE=" + std::to_string(score)
        + (fwd ? " FWD" : " REV");
    }
};

struct read_seq_t {
    bam1_t* read;
    std::string seq, seq_rc;

    read_seq_t() : read(NULL), seq(""), seq_rc("") {}
    read_seq_t(bam1_t* read, std::string& seq) : read(read), seq(seq), seq_rc(seq) {
        get_rc(seq_rc);
    }
};

struct window_t {
    int start, end;

    window_t() : start(0), end(0) {}
    window_t(int start, int end) : start(start), end(end) {}
};

char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
};

bool is_host_region(const region_t *region) {
    return virus_names.count(contig_id2name[region->contig_id]) == 0;
};
bool is_virus_region(const region_t *region) {
    return virus_names.count(contig_id2name[region->contig_id]);
};

bool reads_contain_qname(std::vector<bam1_t*>& reads, std::string qname) {
    for (bam1_t* r : reads) {
        if (std::string(bam_get_qname(r)) == qname) {
            return true;
        }
    }
    return false;
}
bool rr_contain_qname(std::vector<read_realignment_t>& reads, std::string qname) {
    for (read_realignment_t& rr : reads) {
        if (std::string(bam_get_qname(rr.read)) == qname) {
            return true;
        }
    }
    return false;
}


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

void extract_regions(std::vector<bam1_t *> &reads, std::vector<region_t *> &regions, bool process_xa = true) {

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

                anchor_t* a = new anchor_t(xa_dir == '-' ? strand_t::R : strand_t::F, contig_name2id[xa_tname], xa_pos, xa_pos+qlen);
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
    max_dists.push_back(global_stats.max_is);

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
                if (a1 != a2 && anchor_t::can_merge(*a1, *a2, global_stats.max_is) &&
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
                if (anchor_t::can_merge(*new_anchor, *map_it->second, global_stats.max_is)) {
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
            if (a1 != nullptr && a1->end+10 > a2->end) {
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
            regions.push_back(new region_t(a->contig_id, contig_id2tid[a->contig_id],
                                           std::max(0, a->start - global_stats.max_is + a->len()),
                                           std::min(a->end + global_stats.max_is - a->len(),
                                                   (int) chrs[contig_id2name[a->contig_id]].second)));
        };
        delete a;
    }
}


// TODO: unify with isolate_relevant_pairs
std::vector<region_t*>* kmer_index, * kmer_index_rc;
std::vector<read_realignment_t>* reads_per_region, * reads_per_region_rc;
std::mutex* mtx_kmers, * mtx_regions;


inline void insert(ull masked_kmer, region_t* region, bool rc) {
    mtx_kmers[masked_kmer].lock();
    if (rc) {
        kmer_index_rc[masked_kmer].push_back(region);
    } else {
        kmer_index[masked_kmer].push_back(region);
    }
    mtx_kmers[masked_kmer].unlock();
}

inline bool valid_kmer(ull kmer, int len) {
    int count[256];
    count['A'] = count['C'] = count['G'] = count['T'] = 0;
    for (int i = 0; i < len; i++) {
        count[bm_nucl[kmer%4]]++;
        kmer /= 4;
    }

    // filter poly-(ACGT)
    int max_freq = std::max(std::max(count['A'], count['C']), std::max(count['G'], count['T']));
    if (max_freq >= len-2) return false;
    return true;
}


std::mutex mtx;

void index_seq(char* seq, size_t len, region_t* region, bool rc) {
    ull kmer = 0;
    for (int i = 0; i < len; i++) {
        ull nv = nucl_bm[seq[i]];
        kmer = ((kmer << 2) | nv) & KMER_MASK;

        if (i >= KMER_LEN-1) {
            if (valid_kmer(kmer, KMER_LEN)) {
                insert(kmer, region, rc);
            }
        }
    }
}

void index_regions(int id, std::vector<region_t*>* regions, int start, int end) {
    mtx.lock();
    std::cerr << "Thread " << id << ": " << "indexing regions " << start << " to " << end << " " << std::endl;
    mtx.unlock();

    char region_str[100000];
    for (int r = start; r < end; r++) {
        region_t* region = (*regions)[r];
        for (int i = 0; i < region->len(); i++) {
            region_str[i] = toupper(chrs[contig_id2name[region->contig_id]].first[region->start+i]);
        }
        region_str[region->len()] = '\0';
        index_seq(region_str, region->len(), region, false);
        get_rc(region_str, strlen(region_str));
        index_seq(region_str, region->len(), region, true);
    }
}

std::atomic<int> remappings(0);
FILE* scores_file_out_bin;

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


bool accept_alignment(StripedSmithWaterman::Alignment& alignment, std::string& query) {
//    bool poly_ACGT = is_poly_ACGT(query.c_str()+alignment.query_begin, alignment.query_end-alignment.query_begin+1);
    bool long_enough = alignment.ref_end-alignment.ref_begin+1 >= 30;
    bool score_enough = alignment.sw_score >= 30;
    uint32_t c0 = alignment.cigar[0], cl = alignment.cigar[alignment.cigar.size()-1];
    bool left_clipped = cigar_int_to_op(c0) == 'S' && cigar_int_to_len(c0) >= config.min_sc_size;
    bool right_clipped = cigar_int_to_op(cl) == 'S' && cigar_int_to_len(cl) >= config.min_sc_size;
    bool clipped_both_sides = left_clipped && right_clipped;
    return /*!poly_ACGT &&*/ long_enough && score_enough && !clipped_both_sides;
};


void write_alignment_to_bytes(char bytes[16], uint32_t region_id, uint32_t read_id, uint32_t cigar_id, uint16_t ref_begin,
        bool is_rc, uint16_t score) {

    cigar_id = (is_rc ? 0x80000000 : 0) | (cigar_id & 0x7FFFFFFF); // forcing is_rc inside the cigar_id to make 16 bytes

    memcpy(bytes, &region_id, 4);
    memcpy(bytes+4, &read_id, 4);
    memcpy(bytes+8, &cigar_id, 4);
    memcpy(bytes+12, &ref_begin, 2);
    memcpy(bytes+14, &score, 2);
}


void add_realignment_to_region(read_realignment_t& rr, int region_id, int max_score = 0) {
    if (rr.read->core.flag & BAM_FSECONDARY) { // good clip
        if (rr.left_clip_len() >= config.min_sc_size) {
            reads_per_region_rc[region_id].push_back(rr);
        }
        if (rr.right_clip_len() >= config.min_sc_size) {
            reads_per_region[region_id].push_back(rr);
        }
    } else { // the read points toward the bp
        std::vector<read_realignment_t>& reads_per_region_v = rr.rc ? reads_per_region_rc[region_id] : reads_per_region[region_id];
        reads_per_region_v.push_back(rr);
    }
}

inline read_realignment_t remap_read(std::string seq, bam1_t* read, region_t* region, bool is_rc, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
                  region->end - region->start, filter, &alignment, 0);

    if (!accept_alignment(alignment, seq)) return read_realignment_t();

    std::string cigar = alignment_cigar_to_bam_cigar(alignment.cigar);
    mtx.lock();
    if (!cstring_to_cigar.count(cigar)) {
        cstring_to_cigar[cigar] = cigar_str_to_array(cigar);
    }
    auto cigar_v = cstring_to_cigar[cigar];
    mtx.unlock();

    remappings++;
    if (remappings % 1000000 == 0) {
        std::cerr << remappings << " remappings done." << std::endl;
    }

    return read_realignment_t(read, alignment.ref_begin, alignment.ref_end, cigar_v.first, cigar_v.second, alignment.sw_score, is_rc);
}

void associate_reads_to_regions(int id, std::vector<bam1_t*>* host_reads, int start, int end) {
    mtx.lock();
    std::cerr << "Thread " << id << ": " << "considering reads " << start << " to " << end << " " << std::endl;
    mtx.unlock();

    const int MAX_REGIONS = 1000000;
    int* last_read_for_region = new int[MAX_REGIONS];
    int* last_pos = new int[MAX_REGIONS];
    int* last_read_for_region_rc = new int[MAX_REGIONS];
    int* last_pos_rc = new int[MAX_REGIONS];
    std::fill(last_read_for_region, last_read_for_region+MAX_REGIONS, end+1);
    std::fill(last_read_for_region_rc, last_read_for_region_rc+MAX_REGIONS, end+1);

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    for (int r = start+1; r <= end; r++) {
        bam1_t* read = (*host_reads)[r-1];
        std::vector<std::pair<read_realignment_t, int> > realignments;

        std::string seq = get_sequence(read, true);
        std::string seq_rc = seq;
        get_rc(seq_rc);

        ull kmer = 0;
        int len = 0;
        for (int i = 0; i < read->core.l_qseq; i++) {
            len++;

            ull nv = nucl_bm[seq[i]];
            kmer = ((kmer << 2) | nv) & KMER_MASK;

            if (len >= KMER_LEN) {
                /* abs(last_read_for_region[region->id]) != r : we haven't yet matched read and region
                 * last_read_for_region[region->id] == r : we have matched read and region, but we have not aligned them yet
                 * last_read_for_region[region->id] == -r : we have already aligned region and read */
                for (region_t* region : kmer_index[kmer]) {
                    if (last_read_for_region[region->id] == r && i - last_pos[region->id] >= KMER_LEN) {
                        last_read_for_region[region->id] = -r;

                        read_realignment_t rr = remap_read(seq, read, region, false, aligner, filter);
                        if (rr.accepted()) realignments.push_back({rr, region->id});
                    } else if (abs(last_read_for_region[region->id]) != r) {
                        last_read_for_region[region->id] = r;
                        last_pos[region->id] = i;
                    }
                }

                for (region_t* region : kmer_index_rc[kmer]) {
                    if (last_read_for_region_rc[region->id] == r && i - last_pos_rc[region->id] >= KMER_LEN) {
                        last_read_for_region_rc[region->id] = -r;

                        read_realignment_t rr = remap_read(seq_rc, read, region, true, aligner, filter);
                        if (rr.accepted()) realignments.push_back({rr, region->id});
                    } else if (abs(last_read_for_region_rc[region->id]) != r) {
                        last_read_for_region_rc[region->id] = r;
                        last_pos_rc[region->id] = i;
                    }
                }
            }
        }

        uint16_t max_score = 0;
        for (std::pair<read_realignment_t, int>& rr : realignments) {
            max_score = std::max(max_score, rr.first.score);
        }
        for (std::pair<read_realignment_t, int>& rr_and_region : realignments) {
            if (rr_and_region.first.score >= max_score*0.75) {
                mtx_regions[rr_and_region.second].lock();
                add_realignment_to_region(rr_and_region.first, rr_and_region.second, max_score);
                mtx_regions[rr_and_region.second].unlock();
            }
        }
    }

    delete[] last_read_for_region;
    delete[] last_pos;
    delete[] last_read_for_region_rc;
    delete[] last_pos_rc;
}


std::pair<int, int> compute_region_score(region_t* region, bool rc, google::dense_hash_set<std::string>& already_used) {
    int score = 0, reads = 0;
    std::vector<read_realignment_t>& read_realignments = rc ? reads_per_region_rc[region->id] : reads_per_region[region->id];
    for (read_realignment_t& rr : read_realignments) {
        if (!already_used.count(bam_get_qname(rr.read))) {
            score += rr.score;
            reads++;
        }
    }
    return {score, reads};
};

void del_aux(bam1_t* read, const char* tag) {
    uint8_t* b = bam_aux_get(read, tag);
    if (b == NULL) return;
    bam_aux_del(read, b);
}

static int do_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
    new_data = (uint8_t*) realloc(b->data, new_m_data);
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}
static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    uint32_t new_len = b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < b->l_data) {
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return do_realloc_bam_data(b, new_len);
}
int bam_aux_update_int(bam1_t *b, const char tag[2], int64_t val)
{
    uint32_t sz, old_sz = 0, _new = 0;
    uint8_t *s, type;

    if (val < INT32_MIN || val > UINT32_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    if (val < INT16_MIN)       { type = 'i'; sz = 4; }
    else if (val < INT8_MIN)   { type = 's'; sz = 2; }
    else if (val < 0)          { type = 'c'; sz = 1; }
    else if (val < UINT8_MAX)  { type = 'C'; sz = 1; }
    else if (val < UINT16_MAX) { type = 'S'; sz = 2; }
    else                       { type = 'I'; sz = 4; }

    s = bam_aux_get(b, tag);
    if (s) {  // Tag present - how big was the old one?
        switch (*s) {
            case 'c': case 'C': old_sz = 1; break;
            case 's': case 'S': old_sz = 2; break;
            case 'i': case 'I': old_sz = 4; break;
            default: errno = EINVAL; return -1;  // Not an integer
        }
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            s = b->data + b->l_data;
            _new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    if (_new || old_sz < sz) {
        // Make room for new tag
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, (_new ? 3 : 0) + sz - old_sz) < 0)
            return -1;
        s =  b->data + s_offset;
        if (_new) { // Add tag id
            *s++ = tag[0];
            *s++ = tag[1];
        } else {   // Shift following data so we have space
            memmove(s + sz, s + old_sz, b->l_data - s_offset - old_sz);
        }
    } else {
        // Reuse old space.  Data value may be bigger than necessary but
        // we avoid having to move everything else
        sz = old_sz;
        type = (val < 0 ? "\0cs\0i" : "\0CS\0I")[old_sz];
        assert(type > 0);
    }
    *s++ = type;
#ifdef HTS_LITTLE_ENDIAN
    memcpy(s, &val, sz);
#else
    switch (sz) {
        case 4:  u32_to_le(val, s); break;
        case 2:  u16_to_le(val, s); break;
        default: *s = val; break;
    }
#endif
    b->l_data += (_new ? 3 : 0) + sz - old_sz;
    return 0;
}

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

void rc_sequence(bam1_t* read) {
    std::string seq = get_sequence(read);
    get_rc(seq);
    uint8_t* s = bam_get_seq(read);
    for (int i = 0; i < read->core.l_qseq; ++i){
        bam1_seq_seti(s, i, seq_nt16_table[seq[i]]);
    }
}
void set_to_forward(bam1_t* read) {
    if (bam_is_rev(read)) {
        rc_sequence(read);
    }
    read->core.flag &= ~BAM_FREVERSE;
}
void set_to_reverse(bam1_t* read) {
    if (!bam_is_rev(read)) {
        rc_sequence(read);
    }
    read->core.flag |= BAM_FREVERSE;
}
bam1_t* rc_read(bam1_t* read) {
    if (bam_is_rev(read)) set_to_forward(read);
    else set_to_reverse(read);
    return read;
}

void edit_remapped_reads(region_t* region, std::vector<read_realignment_t>& read_realignments, bool rc) {
    for (read_realignment_t rr : read_realignments) {
        bam1_t* read = rr.read;
        read->core.tid = region->original_bam_id;
        read->core.pos = region->start + rr.offset_start;

        if (rr.rc) set_to_reverse(read);
        else set_to_forward(read);

        del_aux(read, "AS");
        del_aux(read, "XS");
        del_aux(read, "XA");
        del_aux(read, "SA");
        del_aux(read, "NM");
        del_aux(read, "MD");

        int l_aux = bam_get_l_aux(read);
        int l_data = read->core.l_qname + 4*rr.cigar_len + (read->core.l_qseq+1)/2
                + read->core.l_qseq + l_aux;
        uint32_t m_data = l_data;
        kroundup32(m_data);
        uint8_t* data = new uint8_t[m_data];
        memset(data, 0, m_data);

        uint8_t* mov_data = data;
        memcpy(mov_data, (uint8_t*) bam_get_qname(read), read->core.l_qname);
        mov_data += read->core.l_qname;
        memcpy(mov_data, rr.cigar, 4*rr.cigar_len);
        mov_data += 4*rr.cigar_len;
        memcpy(mov_data, (uint8_t*) bam_get_seq(read), (read->core.l_qseq+1)/2);
        mov_data += (read->core.l_qseq+1)/2;
        memcpy(mov_data, (uint8_t*) bam_get_qual(read), read->core.l_qseq);
        mov_data += read->core.l_qseq;
        memcpy(mov_data, (uint8_t*) bam_get_aux(read), l_aux);

        read->l_data = l_data;
        read->m_data = m_data;
        read->data = data;
        read->core.n_cigar = rr.cigar_len;
    }
}


void write_qnames_indices(std::string& workspace, std::vector<bam1_t *>& reads) {
    std::ofstream qnames_map_out(workspace + "/qnames-map");
    for (int i = 0; i < reads.size(); i++) {
        qnames_map_out << i << " " << bam_get_qname(reads[i]) << " " << get_sequence(reads[i]) << std::endl;
    }
    qnames_map_out.close();
}
void load_qnames_indices(std::string& workspace, std::vector<bam1_t*> reads, std::vector<bam1_t*>& id_to_read) {
    std::ifstream qnames_map_in(workspace + "/qnames-map");

    std::unordered_map<std::string, std::vector<int> > qname_to_id;
    int i; std::string qname, seq;
    while (qnames_map_in >> i >> qname >> seq) {
        qname_to_id[qname + " " + seq].push_back(i);
    }
    qnames_map_in.close();

    id_to_read.resize(i+1);

    for (bam1_t* r : reads) {
        std::string qname = bam_get_qname(r), seq = get_sequence(r);
        for (int i : qname_to_id[qname + " " + seq]) {
            id_to_read[i] = r;
        }
    }
}
void load_cigars_indices(std::string& workspace, std::vector<std::pair<int, const uint32_t*> >& id_to_cigar) {
    std::ifstream cigars_map_in(workspace + "/cigars-map");
    int i; std::string cigar_str;

    std::vector<std::pair<int, std::string> > temp;
    while (cigars_map_in >> i >> cigar_str) {
        temp.push_back({i, cigar_str});
    }

    id_to_cigar.resize(temp.size()+1);
    for (std::pair<int, std::string>& p : temp) {
        id_to_cigar[p.first] = cigar_str_to_array(p.second);
    }
}


int get_first_max_element(int* v, int len) {
    int max_pos = 0;
    for (int i = 1; i < len; i++) {
        if (v[max_pos] < v[i]) max_pos = i;
    }
    return max_pos;
}
int get_last_max_element(int* v, int len) {
    int max_pos = 0;
    for (int i = 1; i < len; i++) {
        if (v[max_pos] <= v[i]) max_pos = i;
    }
    return max_pos;
}
window_t get_accept_window(std::vector<read_realignment_t>& realigned_reads, bool rc, int region_len,
        google::dense_hash_map<std::string, std::vector<std::string> >& deduped_qnames) {

    int clip_positions[100000];
    uint16_t min_offset = region_len, max_offset = 0;
    std::fill(clip_positions, clip_positions+region_len, 0);
    for (read_realignment_t& rr : realigned_reads) {
        if (!deduped_qnames.count(bam_get_qname(rr.read))) continue;
        if (rc && rr.left_clip_len()) {
            clip_positions[rr.offset_start]++;
        }
        if (!rc && rr.right_clip_len()) {
            clip_positions[rr.offset_end]++;
        }

        min_offset = std::min(min_offset, rr.offset_start);
        max_offset = std::max(max_offset, rr.offset_end);
    }

    // if it is FWD, we want to get the last maximum element
    int bp_pos = rc ? get_first_max_element(clip_positions, region_len) : get_last_max_element(clip_positions, region_len);
    bool clip_at_extremity = bp_pos == (rc ? min_offset : max_offset);

    if (clip_positions[bp_pos] < 2 && !clip_at_extremity) { // cannot get bp from clips, find point of max score density
        int score_points[100000];
        std::fill(score_points, score_points+region_len, 0);
        for (read_realignment_t& rr : realigned_reads) {
            if (!deduped_qnames.count(bam_get_qname(rr.read))) continue;
            score_points[rc ? rr.offset_start : rr.offset_end] += rr.score;
        }
        if (rc) {
            for (int i = region_len-2; i >= 0; i--) {
                score_points[i] += score_points[i+1];
            }
        } else {
            for (int i = 1; i < region_len; i++) {
                score_points[i] += score_points[i-1];
            }
        }

        int score_windows[100000]; // i: total scores of the maxIS-window ending at i
        std::copy(score_points, score_points+region_len, score_windows);
        if (rc) {
            for (int i = region_len-1-global_stats.max_is; i >= 0; i--) {
                score_windows[i] -= score_points[i+global_stats.max_is];
            }
        } else {
            for (int i = global_stats.max_is; i < region_len; i++) {
                score_windows[i] -= score_points[i-global_stats.max_is];
            }
        }

        bp_pos = rc ? get_last_max_element(score_windows, region_len) : get_first_max_element(score_windows, region_len);
    }

    int accept_window_start, accept_window_end;
    if (rc) {
        accept_window_start = bp_pos;
        accept_window_end = accept_window_start + global_stats.max_is;
    } else {
        accept_window_end = bp_pos;
        accept_window_start = accept_window_end - global_stats.max_is;
    }
    return window_t(accept_window_start, accept_window_end);
}


void dedup_reads(std::vector<bam1_t*>& host_reads, std::vector<bam1_t*>& virus_reads) {

    if (host_reads.empty() || virus_reads.empty()) return;

    google::dense_hash_map<std::string, bam1_t*> virus_reads_by_name;
    virus_reads_by_name.set_empty_key("");
    for (bam1_t* read : virus_reads) {
        virus_reads_by_name[bam_get_qname(read)] = read;
    }

    // DEDUP BY SEQUENCE

    // dedup pairs
    dedup_pairs_by_head_qname.set_empty_key("");
    for (bam1_t* hread : host_reads) {
        if (!(hread->core.flag & BAM_FSECONDARY) && virus_reads_by_name.count(bam_get_qname(hread))) {
            bam1_t* vread = virus_reads_by_name[bam_get_qname(hread)];
            std::string joint_string = get_sequence(hread, true) + "-" + get_sequence(vread, true);
            dedup_pairs_by_head_qname[joint_string].push_back(bam_get_qname(hread));
        }
    }
    for (auto& e : dedup_pairs_by_head_qname) {
        std::sort(e.second.begin(), e.second.end());
    }

    google::dense_hash_set<std::string> deduped_qnames;
    deduped_qnames.set_empty_key("");
    for (auto& e : dedup_pairs_by_head_qname) {
        deduped_qnames.insert(e.second[0]);
    }

    for (int i = 0; i < host_reads.size(); i++) {
        if (!deduped_qnames.count(bam_get_qname(host_reads[i]))) {
            bam_destroy1(host_reads[i]);
            host_reads[i] = nullptr;
        }
    }
    host_reads.erase(std::remove(host_reads.begin(), host_reads.end(), nullptr), host_reads.end());
    for (int i = 0; i < virus_reads.size(); i++) {
        if (!deduped_qnames.count(bam_get_qname(virus_reads[i]))) {
            bam_destroy1(virus_reads[i]);
            virus_reads[i] = nullptr;
        }
    }
    virus_reads.erase(std::remove(virus_reads.begin(), virus_reads.end(), nullptr), virus_reads.end());

    std::cerr << "AFTER DEDUP BY SEQ" << std::endl;
    std::cerr << "HOST READS " << host_reads.size() << std::endl;
    std::cerr << "VIRUS READS " << virus_reads.size() << std::endl;


    // DEDUP BY POSITION AND CIGAR

    deduped_qnames.clear();

    // dedup pairs
    std::vector<std::pair<bam1_t*, bam1_t*> > paired_reads;
    for (bam1_t* hread : host_reads) {
        if (virus_reads_by_name.count(bam_get_qname(hread))) {
            bam1_t* vread = virus_reads_by_name[bam_get_qname(hread)];
            paired_reads.push_back({hread, vread});
        }
    }

    // order by host and virus pos, and then by average base quality
    auto pair_cmp = [] (const std::pair<bam1_t*,bam1_t*> p1, const std::pair<bam1_t*,bam1_t*> p2) {
        bam1_t* h1 = p1.first, * h2 = p2.first;
        bam1_t* v1 = p1.second, * v2 = p2.second;
        int avg_qual1 = get_avg_qual(h1, false)+get_avg_qual(v1, false);
        int avg_qual2 = get_avg_qual(h2, false)+get_avg_qual(v2, false);
        return std::tie(h1->core.tid, h1->core.pos, v1->core.tid, v1->core.pos, avg_qual1) <
               std::tie(h2->core.tid, h2->core.pos, v2->core.tid, v2->core.pos, avg_qual2);
    };
    std::sort(paired_reads.begin(), paired_reads.end(), pair_cmp);

    int last_valid_i = 0;
    for (int i = 1; i < paired_reads.size(); i++) {
        bam1_t* h1 = paired_reads[last_valid_i].first,  * h2 = paired_reads[i].first;
        bam1_t* v1 = paired_reads[last_valid_i].second, * v2 = paired_reads[i].second;
        std::string h1_cigar = get_cigar_code(h1), v1_cigar = get_cigar_code(v1);
        std::string h2_cigar = get_cigar_code(h2), v2_cigar = get_cigar_code(v2);

        if (std::tie(h1->core.tid, h1->core.pos, h1_cigar, v1->core.tid, v1->core.pos, v1_cigar) ==
            std::tie(h2->core.tid, h2->core.pos, h2_cigar, v2->core.tid, v2->core.pos, v2_cigar)) {
            // keep last_valid_i
            std::string last_valid_seq = get_sequence(h1, true) + "-" + get_sequence(v1, true);
            std::string curr_seq = get_sequence(h2, true) + "-" + get_sequence(v2, true);
            if (dedup_pairs_by_head_qname.count(curr_seq) && dedup_pairs_by_head_qname.count(last_valid_seq)) {
                std::vector<std::string>& pairs_by_curr_seq = dedup_pairs_by_head_qname[curr_seq];
                std::vector<std::string>& pairs_by_last_valid_seq = dedup_pairs_by_head_qname[last_valid_seq];
                pairs_by_last_valid_seq.insert(pairs_by_last_valid_seq.end(), pairs_by_curr_seq.begin(), pairs_by_curr_seq.end());
                pairs_by_curr_seq.clear();
            }
        } else {
            deduped_qnames.insert(bam_get_qname(h1));
            last_valid_i = i;
        }
    }
    std::cerr << last_valid_i << std::endl;
    deduped_qnames.insert(bam_get_qname(paired_reads[last_valid_i].first));

    for (int i = 0; i < host_reads.size(); i++) {
        if (!deduped_qnames.count(bam_get_qname(host_reads[i]))) {
            bam_destroy1(host_reads[i]);
            host_reads[i] = nullptr;
        }
    }
    host_reads.erase(std::remove(host_reads.begin(), host_reads.end(), nullptr), host_reads.end());
    for (int i = 0; i < virus_reads.size(); i++) {
        if (!deduped_qnames.count(bam_get_qname(virus_reads[i]))) {
            bam_destroy1(virus_reads[i]);
            virus_reads[i] = nullptr;
        }
    }
    virus_reads.erase(std::remove(virus_reads.begin(), virus_reads.end(), nullptr), virus_reads.end());

    std::cerr << "AFTER DEDUP BY POS" << std::endl;
    std::cerr << "HOST READS " << host_reads.size() << std::endl;
    std::cerr << "VIRUS READS " << virus_reads.size() << std::endl;
}


google::dense_hash_map<std::string, int> get_distances(std::vector<read_realignment_t>& realignments, int_breakpoint_t& bp, region_t region) {
    google::dense_hash_map<std::string, int> dists;
    dists.set_empty_key("");
    for (read_realignment_t& rr : realignments) {
        if (rr.read->core.flag & BAM_FSECONDARY) continue;
        std::string qname = bam_get_qname(rr.read);

        int read_endpos;
        if (bp.fwd) { // include spurious soft-clipping as part of the distance
            read_endpos = region.start + rr.offset_start - rr.left_clip_len();
        } else {
            read_endpos = region.start + rr.offset_end + rr.right_clip_len();
        }

        int dist = abs(read_endpos - bp.bp());
        if (!dists.count(qname) || dists[qname] < dist) {
            dists[qname] = dist;
        }
    }
    return dists;
}
google::dense_hash_map<std::string, int> get_distances(std::vector<read_realignment_t>& realignments, int_breakpoint_t& bp) {
    return get_distances(realignments, bp, region_t(0,0,0,0));
}
google::dense_hash_map<std::string, std::vector<std::string> > dedup_reads(
        std::vector<read_realignment_t>& host_realignments, std::vector<read_realignment_t> virus_realignments,
        region_t* host_region, bool host_fwd, region_t* virus_region, bool virus_fwd) {
    if (host_realignments.empty() || virus_realignments.empty()) {
        google::dense_hash_map<std::string, std::vector<std::string> > empty_set;
        empty_set.set_empty_key("");
        return empty_set;
    }

    // remove virus realignments without a host realignment (TODO: necessary?)
    std::unordered_map<std::string, read_realignment_t> host_reads_by_name;
    for (read_realignment_t& rr : host_realignments) {
        if (rr.read->core.flag & BAM_FSECONDARY) continue;
        host_reads_by_name[bam_get_qname(rr.read)] = rr;
    }
    virus_realignments.erase(
            std::remove_if(virus_realignments.begin(), virus_realignments.end(),
                           [&host_reads_by_name] (const read_realignment_t& rr) {
                               return !host_reads_by_name.count(bam_get_qname(rr.read));
                           }), virus_realignments.end());

    // in case there are two reads with the same qname (i.e. a mate and an anchor/clip), choose the farthest from the bp
    int_breakpoint_t host_bp(host_region->contig_id, 0, host_region->len(), host_fwd);
    google::dense_hash_map<std::string, int> h_dists = get_distances(host_realignments, host_bp);
    int_breakpoint_t virus_bp(virus_region->contig_id, 0, virus_region->len(), virus_fwd);
    google::dense_hash_map<std::string, int> v_dists = get_distances(virus_realignments, virus_bp);

    std::sort(virus_realignments.begin(), virus_realignments.end(),
              [] (const read_realignment_t& rr1, const read_realignment_t& rr2) {
                    return rr1.score > rr2.score;
              });

    std::vector<std::string> h_v_seqs;
    for (read_realignment_t& vrr : virus_realignments) {
        std::string qname = bam_get_qname(vrr.read);

        read_realignment_t& hrr = host_reads_by_name[qname];
        std::string seq = get_sequence(hrr.read, true) + "-" + get_sequence(vrr.read, true);
        h_v_seqs.push_back(seq);
    }

    google::dense_hash_map<std::string, std::vector<std::string> > qnames_to_keep;
    qnames_to_keep.set_empty_key("");
    for (int i = 0; i < virus_realignments.size(); i++) {
        std::string qname_i = bam_get_qname(virus_realignments[i].read);
        int host_i_is = h_dists[qname_i];
        int virus_i_is = v_dists[qname_i];

        bool has_prev_dup = false;
        std::string seq_i = h_v_seqs[i];
        for (int j = 0; j < i; j++) {
            std::string qname_j = bam_get_qname(virus_realignments[j].read);
            int host_j_is = h_dists[qname_j];
            int virus_j_is = v_dists[qname_j];

            if (host_i_is == host_j_is && virus_i_is == virus_j_is) {
                has_prev_dup = true;
                qnames_to_keep[qname_j].insert(qnames_to_keep[qname_j].end(),
                        dedup_pairs_by_head_qname[seq_i].begin(), dedup_pairs_by_head_qname[seq_i].end());
                break;
            }
        }

        if (!has_prev_dup) {
            qnames_to_keep[qname_i] = dedup_pairs_by_head_qname[seq_i];
        }
    }
    return qnames_to_keep;
}


const int APPROX_FACTOR = 500;
std::unordered_set<std::string> existing_caches;
std::unordered_map<std::string, google::dense_hash_map<uint64_t, read_realignment_t> > cached_virus_alignments;

uint64_t virus_read_id = 1;
std::vector<read_seq_t> virus_reads_seqs;

std::pair<int,int> remap_virus_reads_supp(region_t* virus_region, bool vregion_rc, std::vector<uint64_t>& read_seq_ids,
        region_t* host_region, bool hregion_rc, std::vector<read_realignment_t>& host_alignemnts,
        std::vector<read_realignment_t>* virus_read_realignments,
        StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter,
        window_t& host_window, window_t& virus_window) {

    std::string region_str = (vregion_rc ? "-" : "+") + virus_region->to_str();

    google::dense_hash_map<uint64_t, read_realignment_t>& cached_region_alignments = cached_virus_alignments[region_str];
    if (!existing_caches.count(region_str)) {
        cached_region_alignments.set_empty_key(0);
        existing_caches.insert(region_str);
    }

    std::vector<read_realignment_t> host_read_realignments_local, virus_read_realignments_local;
    for (int i = 0; i < read_seq_ids.size(); i++) {
        uint64_t read_seq_id = read_seq_ids[i];
        if (!cached_region_alignments.count(read_seq_id)) {
            StripedSmithWaterman::Alignment alignment;
            read_seq_t &read_seq = virus_reads_seqs[read_seq_id];
            std::string &read_seq_str = vregion_rc ? read_seq.seq_rc : read_seq.seq;
            aligner.Align(read_seq_str.c_str(), chrs[contig_id2name[virus_region->contig_id]].first + virus_region->start,
                          virus_region->len(), filter, &alignment, 0);

            if (accept_alignment(alignment, read_seq_str)) {
                std::string cigar_str = alignment_cigar_to_bam_cigar(alignment.cigar);
                std::pair<int, const uint32_t *> cigar = cigar_str_to_array(cigar_str);
                cached_region_alignments[read_seq_id] = read_realignment_t(read_seq.read, alignment.ref_begin, alignment.ref_end,
                                                                         cigar.first, cigar.second, alignment.sw_score, vregion_rc);
            } else {
                // cache failure to align
                cached_region_alignments[read_seq_id] = read_realignment_t(read_seq.read);
            }
        }

        read_realignment_t& rr = cached_region_alignments[read_seq_id];
        if (!rr.accepted()) continue;

        host_read_realignments_local.push_back(host_alignemnts[i]);
        virus_read_realignments_local.push_back(rr);
    }

    auto remove_rr = [&virus_region] (read_realignment_t& rr, window_t accept_window, bool rc) {
        int MIN_QUAL = 10;
        bool remove = rr.offset_start < accept_window.start || rr.offset_end > accept_window.end;
        if (rr.left_clip_len() > config.max_sc_dist && rr.left_clip_avg_qual() >= MIN_QUAL) {
            remove |= abs(rr.offset_start - accept_window.start) > config.max_sc_dist;
        }
        if (rr.right_clip_len() > config.max_sc_dist && rr.right_clip_avg_qual() >= MIN_QUAL) {
            remove |= abs(rr.offset_end - accept_window.end) > config.max_sc_dist;
        }
        return remove;
    };

    google::dense_hash_map<std::string, std::vector<std::string> > deduped_qnames = dedup_reads(
            host_read_realignments_local, virus_read_realignments_local, host_region, !hregion_rc, virus_region, !vregion_rc);

    google::dense_hash_set<std::string> to_be_removed;
    to_be_removed.set_empty_key("");
    host_window = get_accept_window(host_read_realignments_local, hregion_rc, host_region->len(), deduped_qnames);
    virus_window = get_accept_window(virus_read_realignments_local, vregion_rc, virus_region->len(), deduped_qnames);
    for (read_realignment_t& rr : host_read_realignments_local) {
        if (remove_rr(rr, host_window, hregion_rc)) {
            to_be_removed.insert(bam_get_qname(rr.read));
        }
    }
    for (read_realignment_t& rr : virus_read_realignments_local) {
        if (remove_rr(rr, virus_window, vregion_rc)) {
            to_be_removed.insert(bam_get_qname(rr.read));
        }
    }

    int h_score = 0, v_score = 0;
    for (int i = 0; i < virus_read_realignments_local.size(); i++) {
        read_realignment_t& hrr = host_read_realignments_local[i];
        read_realignment_t& vrr = virus_read_realignments_local[i];
        if (to_be_removed.count(bam_get_qname(hrr.read))) continue;

        if (deduped_qnames.count(bam_get_qname(hrr.read))) {
            h_score += hrr.score;
            v_score += vrr.score;
        }
        if (virus_read_realignments != NULL) {
            virus_read_realignments->push_back(vrr);
        }
    }

    return {h_score, v_score};
}

std::pair<region_t*, bool> remap_virus_reads(region_score_t& host_region_score, std::vector<read_realignment_t>& virus_read_realignments,
        google::dense_hash_map<std::string, bam1_t*>& virus_reads_by_name, window_t& host_window, window_t& virus_window) {
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;

    std::vector<read_realignment_t>& host_read_realignments = host_region_score.fwd ?
            reads_per_region[host_region_score.region->id] : reads_per_region_rc[host_region_score.region->id];

    // retrieve virus mates
    std::vector<read_realignment_t> remapped_host_reads;
    std::vector<bam1_t*> remapped_reads_mates;
    std::vector<read_realignment_t> remapped_host_anchors;
    for (read_realignment_t& rr : host_read_realignments) {
        std::string qname = bam_get_qname(rr.read);
        // we call good clipped read a read where the anchor maps to host and the clip to virus, or viceversa
        // when this happens, we often create two copies of the read, one in host and one in anchor
        // because of this, there are often two reads in host: an unclipped one, and a good anchor (or its copy)
        // in such case, the anchor is marked as secondary, because it should not use its own copy in the virus as a mate
        if (!(rr.read->core.flag & BAM_FSECONDARY) && virus_reads_by_name.count(qname)) {
            remapped_host_reads.push_back(rr);
            remapped_reads_mates.push_back(virus_reads_by_name[qname]);
        }

        if ((host_region_score.fwd && rr.right_clip_len())
        || (!host_region_score.fwd && rr.left_clip_len())) {
            remapped_host_anchors.push_back(rr);
        }

    }

    // extract virus regions
    std::vector<region_t*> virus_regions;
    extract_regions(remapped_reads_mates, virus_regions);
    for (int i = 0; i < virus_regions.size(); i++) {
        if (is_host_region(virus_regions[i])) {
            delete virus_regions[i];
            virus_regions[i] = nullptr;
        }
    }
    virus_regions.erase(std::remove(virus_regions.begin(), virus_regions.end(), nullptr), virus_regions.end());
    if (virus_regions.empty()) return {nullptr, false};

    for (region_t* region : virus_regions) {
        region->start = (region->start/APPROX_FACTOR) * APPROX_FACTOR;
        region->end = (region->end/APPROX_FACTOR) * APPROX_FACTOR + APPROX_FACTOR;
        region->end = std::min(region->end, (int) chrs[contig_id2name[region->contig_id]].second);
    }

    std::vector<uint64_t> read_seq_ids;
    std::vector<read_realignment_t> host_alignments;
    read_seq_ids.reserve(remapped_reads_mates.size()+remapped_host_anchors.size());
    for (int i = 0; i < remapped_reads_mates.size(); i++) {
        bam1_t* read = remapped_reads_mates[i];
        if (!read->id) {
            read->id = virus_read_id++;
            std::string seq = get_sequence(read, true);
            virus_reads_seqs.push_back(read_seq_t(read, seq));
        }
        read_seq_ids.push_back(read->id);
        host_alignments.push_back(remapped_host_reads[i]);
    }
    for (read_realignment_t& anchor : remapped_host_anchors) {
        if (!anchor.read->id) {
            anchor.read->id = virus_read_id++;

            std::string qname = bam_get_qname(anchor.read);

            // create a "regular" mate out of the clipped portion
            bam1_t* rc_read = bam_dup1(anchor.read);
            if (host_region_score.fwd && anchor.right_clip_len() && !anchor.rc) {
                set_to_forward(rc_read);
                rc_read->core.flag |= BAM_FREVERSE;
            }
            else if (!host_region_score.fwd && anchor.left_clip_len() && anchor.rc) {
                set_to_reverse(rc_read);
                rc_read->core.flag &= ~BAM_FREVERSE;
            }

            std::string clip = get_sequence(rc_read, true);
            virus_reads_seqs.push_back(read_seq_t(rc_read, clip));
        }
        read_seq_ids.push_back(anchor.read->id);
        host_alignments.push_back(anchor);
    }

    region_t* best_region = NULL;
    std::pair<int, int> best_score = {-1, -1};
    bool is_best_vregion_fwd;
    for (region_t* virus_region : virus_regions) {
        window_t temp1, temp2;
        std::pair<int, int> score = remap_virus_reads_supp(virus_region, false, read_seq_ids, host_region_score.region, !host_region_score.fwd,
                host_alignments, NULL, aligner, filter, temp1, temp2);
        std::pair<int, int> rc_score = remap_virus_reads_supp(virus_region, true, read_seq_ids, host_region_score.region, !host_region_score.fwd,
                host_alignments, NULL, aligner, filter, temp1, temp2);

        std::pair<int, int> max_score = std::max(score, rc_score);
        if (best_score < max_score) {
            best_score = max_score;
            delete best_region;
            best_region = virus_region;
            is_best_vregion_fwd = (max_score == score);
        } else {
            delete virus_region;
        }
    }

    // FIXME: what if best_region is null?
    std::cerr << "BEST VIRUS REGION: " << best_region->to_str() << " " << best_score.first << std::endl;
    remap_virus_reads_supp(best_region, !is_best_vregion_fwd, read_seq_ids, host_region_score.region, !host_region_score.fwd,
            host_alignments, &virus_read_realignments, aligner, filter, host_window, virus_window);

    google::dense_hash_map<std::string, std::vector<std::string> > deduped_reads = dedup_reads(
            host_alignments, virus_read_realignments, host_region_score.region, host_region_score.fwd, best_region, is_best_vregion_fwd);

    std::cerr << "BEFORE " << host_region_score.to_str() << std::endl;
    google::dense_hash_set<std::string> remapped_virus_qnames;
    remapped_virus_qnames.set_empty_key("");
    for (read_realignment_t& rr : virus_read_realignments) {
        remapped_virus_qnames.insert(bam_get_qname(rr.read));
    }
    host_region_score.score = 0;
    for (read_realignment_t& rr : host_read_realignments) {
        if (remapped_virus_qnames.count(bam_get_qname(rr.read)) && deduped_reads.count(bam_get_qname(rr.read))) {
            host_region_score.score += rr.score;
        }
    }
    std::cerr << "AFTER " << host_region_score.to_str() << std::endl << std::endl;

    return {best_region, is_best_vregion_fwd};
}


double mean(std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double std_dev(std::vector<double>& v) {
    double m = mean(v);
    double var = 0;
    for (double d : v) {
        var += (d-m) * (d-m);
    }
    var /= v.size();
    return std::sqrt(var);
}
struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1, T2> const &pair) const
    {
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);
        return h1 ^ h2;
    }
};
int overlap(int s1, int e1, int s2, int e2) {
    return std::max(0, std::min(e1, e2) - std::max(s1, s2) + 1);
}
int calculate_coverage(std::vector<read_realignment_t> reads_realignments) {
    if (reads_realignments.empty()) return 0;

    std::sort(reads_realignments.begin(), reads_realignments.end(),
            [](const read_realignment_t& rr1, const read_realignment_t& rr2) {
        return rr1.offset_start < rr2.offset_start;
    });

    uint16_t coverage = 0;
    uint16_t curr_start = 0, curr_end = 0;
    for (int i = 0; i < reads_realignments.size(); i++) {
        read_realignment_t& rr = reads_realignments[i];
        if (rr.suspicious || rr.too_long) continue;

        if (overlap(curr_start, curr_end, rr.offset_start, rr.offset_end) > 0) { // overlaps, extend if possible
            curr_end = std::max(curr_end, rr.offset_end);
        } else {
            coverage += curr_end - curr_start + 1;
            curr_start = rr.offset_start;
            curr_end = rr.offset_end;
        }
    }
    coverage += curr_end - curr_start + 1;

    return coverage;
}

std::pair<std::vector<double>, std::vector<double> > gen_population_for_ins(std::string bam_fname, std::string workdir,
        stats_t bam_stats) {
    open_samFile_t* open_sam = open_samFile(bam_fname.c_str());

    std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
    std::vector<double> population_host, population_virus;
    bam1_t* read = bam_init1();
    char contig[1000]; int pos;

    std::vector<double> cov2_dist, cov3_dist, cov4_dist;
    while (rnd_pos_fin >> contig >> pos) {
        std::vector<double>& population = virus_names.count(contig) ? population_virus : population_host;
        if (population.size() > POP_SIZE) break;

        std::unordered_set<std::pair<int, int>, pair_hash> used_pairs;

        char region[1000];
        sprintf(region, "%s:%d-%d", contig, pos-bam_stats.max_is, pos);
        hts_itr_t* iter = sam_itr_querys(open_sam->idx, open_sam->header, region);
        while (sam_itr_next(open_sam->file, iter, read) >= 0) {
            if (!is_valid(read) || bam_is_rev(read)) continue;
            int mendpos = get_mate_endpos(read);
            if (read->core.pos < pos-bam_stats.max_is || mendpos > pos+bam_stats.max_is) continue;

            if (!bam_is_rev(read) && bam_is_mrev(read) && read->core.isize > 0 && read->core.isize >= bam_stats.min_is
                && read->core.isize <= 2*bam_stats.max_is) {

                if (used_pairs.count({read->core.pos, mendpos})) continue; // we already sampled a duplicate pair

                int start = read->core.pos+30, end = mendpos-30;
                if (start >= end) continue;

                if (start <= pos && pos <= end) {
                    used_pairs.insert({read->core.pos, mendpos});
                    population.push_back(pos-read->core.pos);
                    population.push_back(mendpos-pos);
                }
            }
        }
        sam_itr_destroy(iter);
    }
    close_samFile(open_sam);
    bam_destroy1(read);

    return {population_host, population_virus};
}

google::dense_hash_set<std::string> mismatch_filter(region_t* host_region, std::vector<read_realignment_t>& hrrs) {

    std::sort(hrrs.begin(), hrrs.end(),
            [](const read_realignment_t& rr1, const read_realignment_t& rr2) {
        return std::tie(rr1.offset_start, rr1.offset_end) < std::tie(rr2.offset_start, rr2.offset_end);
    });

    std::vector<std::string> read_sequences;
    for (read_realignment_t& hrr : hrrs) {
        std::string temp_sequence = get_sequence(hrr.read, true);
        if (hrr.rc) get_rc(temp_sequence);

        char read_sequence[10000];
        int pos_dst = 0, pos_src = 0;
        for (int i = 0; i < hrr.cigar_len; i++) {
            char op = bam_cigar_opchr(hrr.cigar[i]);
            if (op == 'M') {
                int len = bam_cigar_oplen(hrr.cigar[i]);
                strncpy(read_sequence + pos_dst, temp_sequence.c_str() + pos_src, len);
                pos_dst += len;
                pos_src += len;
            }   else if (op == 'D') {
                int len = bam_cigar_oplen(hrr.cigar[i]);
                memset(read_sequence+pos_dst, 'D', len);
                pos_dst += len;
            } else if (op == 'I') {
                int len = bam_cigar_oplen(hrr.cigar[i]);
                memset(read_sequence+pos_dst, 'I', len);
                pos_src += len;
                pos_dst += len;
            } else if (op == 'S') {
                pos_src += bam_cigar_oplen(hrr.cigar[i]);
            }
        }
        read_sequence[pos_dst] = '\0';
        read_sequences.push_back(read_sequence);

        std::cerr << bam_get_qname(hrr.read) << " " << cigar_array_to_str(hrr.cigar_len, hrr.cigar) << " " << read_sequence << std::endl;
    }


    // find homopolymers
    // TODO: if slow, implement a linear version
    int count[256];
    std::vector<std::pair<int, int> > homopolymers;
    for (std::string seq : read_sequences) {
        std::pair<int, int> curr_homopolymer = {-1, -1};
        int best_score = 0;
        for (int i = 0; i < seq.length(); i++) {
            count['A'] = count['C'] = count['G'] = count['T'] = 0;
            for (int j = i; j < seq.length(); j++) {
                if (seq[j] == 'A' || seq[j] == 'C' || seq[j] == 'G' || seq[j] == 'T') {
                    count['A']--;
                    count['C']--;
                    count['G']--;
                    count['T']--;
                    count[seq[j]] += 2;

                    int curr_score = max(count['A'], count['C'], count['G'], count['T']);
                    if (curr_score < 10 || curr_score < 0.8*(j-i+1)) continue;

                    if (curr_score > best_score || (curr_score == best_score && j-i < curr_homopolymer.second-curr_homopolymer.first)) {
                        best_score = curr_score;
                        curr_homopolymer = {i,j};
                    }
                }
            }
        }

        homopolymers.push_back(curr_homopolymer);
    }

    for (int i = 0; i < homopolymers.size(); i++) {
        std::pair<int, int> hp = homopolymers[i];
        if (hp.first != -1) {
            std::cerr << "HOMO " << bam_get_qname(hrrs[i].read) << " " << read_sequences[i] << " "
                << read_sequences[i].substr(hp.first, hp.second-hp.first+1) << std::endl;
        }
    }

    bool is_indel[255];
    std::fill(is_indel, is_indel+255, false);
    is_indel['D'] = is_indel['I'] = true;

    // mismatch filter
    int s = 0, e = 0;
    std::vector<int> read_offsets, read_errors;
    char prev_consensus = ' ';
    for (int pos = 0; pos < host_region->len(); pos++) {
        while (s < hrrs.size() && hrrs[s].offset_end < pos) s++;
        while (e < hrrs.size() && hrrs[e].offset_start < pos) {
            read_offsets.push_back(0);
            read_errors.push_back(0);
            e++;
        }

        int a = 0, c = 0, g = 0, t = 0, d = 0, ins = 0, tot = 0;
        for (int i = s; i < e; i++) {
            if (read_offsets[i] >= read_sequences[i].length()) continue;
            char base = read_sequences[i][read_offsets[i]];
            if (base == 'A') a++;
            else if (base == 'C') c++;
            else if (base == 'G') g++;
            else if (base == 'T') t++;
            else if (base == 'D') d++;
            else if (base == 'I') ins++;
            tot++;
        }

        char consensus = 'N';
        if (a > max(c,g,t,d)) consensus = 'A';
        else if (c > max(a,g,t,d,ins)) consensus = 'C';
        else if (g > max(a,c,t,d,ins)) consensus = 'G';
        else if (t > max(a,c,g,d,ins)) consensus = 'T';
        else if (d > max(a,c,g,t,ins)) consensus = 'D';
        else if (ins > max(a,c,g,t,d)) consensus = 'I';

        for (int i = s; i < e; i++) {
            if (read_offsets[i] >= read_sequences[i].length()) continue;

            // erros in homopolymer are fine - see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4719071/
            if (read_offsets[i] >= homopolymers[i].first && read_offsets[i] <= homopolymers[i].second) {
                read_offsets[i]++;
                continue;
            }

            char base = read_sequences[i][read_offsets[i]];
            if (base != consensus) {
                read_errors[i]++;
                if ((is_indel[base] && read_offsets[i] > 0 && read_sequences[i][read_offsets[i]] != base)
                ||  (is_indel[consensus] && prev_consensus != consensus)) { // indel error is penalized more than just a mismatch
                    read_errors[i]++;
                }
            }
            read_offsets[i]++;
        }

        if (consensus == 'I') pos--;

        prev_consensus = consensus;
    }

    google::dense_hash_set<std::string> qnames_to_remove;
    qnames_to_remove.set_empty_key("");
    for (int i = 0; i < hrrs.size(); i++) {
        if (read_errors[i] > (hrrs[i].offset_end-hrrs[i].offset_start+1)*0.06) {
            qnames_to_remove.insert(bam_get_qname(hrrs[i].read));

            std::cerr << "EPURATING " << bam_get_qname(hrrs[i].read) << " " << cigar_array_to_str(hrrs[i].cigar_len, hrrs[i].cigar);
            std::cerr << " from ID=" << virus_integration_id;
            std::cerr << " because it has " << read_errors[i] << " errors." << std::endl;
        }
    }

    return qnames_to_remove;
}

int fix_tid(open_samFile_t* sam_file, bam1_t* read) {
    std::string target_name = sam_file->header->target_name[read->core.tid];
    return contig_id2tid[contig_name2id[target_name]];
}

int main(int argc, char* argv[]) {

    nucl_bm['A'] = nucl_bm['a'] = 0;
    nucl_bm['C'] = nucl_bm['c'] = 1;
    nucl_bm['G'] = nucl_bm['g'] = 2;
    nucl_bm['T'] = nucl_bm['t'] = 3;
    nucl_bm['N'] = nucl_bm['n'] = 0;

    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';

    std::string host_reference_fname  = argv[1];
    std::string virus_reference_fname  = argv[2];
    std::string workdir = argv[3];
    std::vector<std::string> workspaces, bam_fnames;
    for (int i = 4; i < argc; i+=2) {
        bam_fnames.push_back(argv[i]);
        workspaces.push_back(argv[i+1]);
    }

    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }


    FILE* fasta = fopen(host_reference_fname.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(fasta));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        chrs[std::string(seq->name.s)] = {new char[seq->seq.l+1], seq->seq.l};
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
    }
    kseq_destroy(seq);
    fclose(fasta);


    // read virus list
    fasta = fopen(virus_reference_fname.c_str(), "r");
    seq = kseq_init(fileno(fasta));
    while ((l = kseq_read(seq)) >= 0) {
        virus_names[seq->name.s] = seq->seq.l;
        chrs[std::string(seq->name.s)] = {new char[seq->seq.l + seq->seq.l/2 + 1], seq->seq.l + seq->seq.l/2};
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
        strncpy(chrs[std::string(seq->name.s)].first+seq->seq.l, seq->seq.s, seq->seq.l/2);
    }
    kseq_destroy(seq);
    fclose(fasta);

    config = parse_config(workdir + "/config.txt");
    global_stats.min_is = INT32_MAX;
    global_stats.max_is = 0;
    for (std::string workspace : workspaces) {
        stats_t stats = parse_stats(workspace + "/stats.txt");
        bam_stats.push_back(stats);
        global_stats.min_is = std::min(global_stats.min_is, stats.min_is);
        global_stats.max_is = std::max(global_stats.max_is, stats.max_is);
    }

    open_samFile_t* host_and_virus_file = open_samFile((workspaces[0] + "/retained-pairs-remapped.sorted.bam").data()); // just for the header
    for (int i = 0; i < host_and_virus_file->header->n_targets; i++) {
        int contig_id = contig_name2id[host_and_virus_file->header->target_name[i]];
        contig_id2tid[contig_id] = i;
        contig_tid2id[i] = contig_id;
    }

    for (int i = 0; i < workspaces.size(); i++) {
        std::string workspace = workspaces[i];
        std::pair<std::vector<double>, std::vector<double> > populations = gen_population_for_ins(bam_fnames[i], workdir, bam_stats[i]);
        host_maxIS = std::max(host_maxIS, (int) (mean(populations.first) + 3*std_dev(populations.first)));
    }


    /* === READ CHIMERIC PAIRS/READS === */

    /* == Pouring reads which are not good clips into host_reads and virus_reads == */

    bam1_t* read = bam_init1();

    // extracting paired host-virus reads
    std::vector<bam1_t*> host_reads, virus_reads;
    for (int i = 0; i < workspaces.size(); i++) {
        std::string workspace = workspaces[i];

        std::unordered_set<std::string> host_qnames;
        open_samFile_t* host_file = open_samFile((workspace + "/host-side.sorted.bam").data(), true);
        while (sam_read1(host_file->file, host_file->header, read) >= 0) {
            bam_aux_update_int(read, "BF", i);
            host_reads.push_back(bam_dup1(read));
            host_qnames.insert(bam_get_qname(read));
        }
        close_samFile(host_file);

        std::unordered_set<std::string> virus_qnames;
        open_samFile_t* virus_file = open_samFile((workspace + "/virus-side.sorted.bam").data(), true);
        while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
            bam_aux_update_int(read, "BF", i);
            virus_reads.push_back(bam_dup1(read));
            virus_qnames.insert(bam_get_qname(read));
        }
        close_samFile(virus_file);

        /* == */

        std::unordered_map<std::string, std::vector<std::pair<bam1_t*, bool> > > good_clipped_pairs;
        // pairs where both are good clipped reads. Bool = true indicates a host read

        /* == Reading host anchors == */

        // pos of good clips for first (resp. second) reads
        std::unordered_map<std::string, std::pair<int, int> > good_clips_pos[2];

        open_samFile_t* host_clips_file = open_samFile((workspace + "/host-clips.sorted.bam").data(), true);
        while (sam_read1(host_clips_file->file, host_clips_file->header, read) >= 0) {
            std::string qname = bam_get_qname(read);
            bool first = qname[qname.length()-1] == '1';
            qname = qname.substr(0, qname.length()-4);
            good_clips_pos[first ? 0 : 1][qname] = {fix_tid(host_clips_file, read), read->core.pos};
        }
        close_samFile(host_clips_file);

        open_samFile_t* virus_clips_file = open_samFile((workspace + "/virus-clips.sorted.bam").data(), true);
        while (sam_read1(virus_clips_file->file, virus_clips_file->header, read) >= 0) {
            std::string qname = bam_get_qname(read);
            bool first = qname[qname.length()-1] == '1';
            qname = qname.substr(0, qname.length()-4);
            good_clips_pos[first ? 0 : 1][qname] = {fix_tid(virus_clips_file, read), read->core.pos};
        }
        close_samFile(virus_clips_file);

        open_samFile_t* host_anchors_file = open_samFile((workspace + "/host-anchors.sorted.bam").data(), true);
        while (sam_read1(host_anchors_file->file, host_anchors_file->header, read) >= 0) {
            bam_aux_update_int(read, "BF", i);
            if (host_qnames.count(bam_get_qname(read))) {
                bam1_t* virus_copy = bam_dup1(read);
                auto vc_pos = good_clips_pos[read->core.flag & BAM_FREAD1 ? 0 : 1][bam_get_qname(read)];
                virus_copy->core.tid = vc_pos.first;
                virus_copy->core.pos = vc_pos.second;
                virus_reads.push_back(virus_copy);
                bam1_t* host_copy = bam_dup1(read);
                host_copy->core.flag |= BAM_FSECONDARY;
                host_reads.push_back(host_copy);
            } else if (virus_qnames.count(bam_get_qname(read))) {
                host_reads.push_back(bam_dup1(read));
            } else {
                good_clipped_pairs[bam_get_qname(read)].push_back({bam_dup1(read), true});
            }
        }
        close_samFile(host_anchors_file);

        /* == */

        /* == Reading virus anchors == */

        open_samFile_t* virus_anchors_file = open_samFile((workspace + "/virus-anchors.sorted.bam").data(), true);
        while (sam_read1(virus_anchors_file->file, virus_anchors_file->header, read) >= 0) {
            bam_aux_update_int(read, "BF", i);
            if (host_qnames.count(bam_get_qname(read))) {
                virus_reads.push_back(bam_dup1(read));
                bam1_t* rc_copy = rc_read(bam_dup1(read));
                rc_copy->core.flag |= BAM_FSECONDARY;
                host_reads.push_back(rc_copy);
            } else if (virus_qnames.count(bam_get_qname(read))) {
                bam1_t* copy = bam_dup1(read);
                auto hc_pos = good_clips_pos[read->core.flag & BAM_FREAD1 ? 0 : 1][bam_get_qname(read)];
                copy->core.tid = hc_pos.first;
                copy->core.pos = hc_pos.second;
                host_reads.push_back(copy);
            } else {
                good_clipped_pairs[bam_get_qname(read)].push_back({bam_dup1(read), false});
            }
        }
        close_samFile(virus_anchors_file);

        /* == */

        for (auto& e : good_clipped_pairs) {
            if (e.second.size() == 2) {
                bam1_t *r1 = e.second[0].first, *r2 = e.second[1].first;
                bool r1_is_host = e.second[0].second, r2_is_host = e.second[1].second;
                char r1_dir = bam_aux2A(bam_aux_get(r1, "CD")), r2_dir = bam_aux2A(bam_aux_get(r2, "CD"));

                // r1 will be the host read, r2 the virus read and the secondary host read
                auto point_twd_bp = [](bam1_t *r, char dir) {
                    return (bam_is_rev(r) && dir == 'L') || (!bam_is_rev(r) && dir == 'R');
                };

                if (r1_is_host != r2_is_host) {
                    // if one is from host and one from virus, then it's easy
                    host_reads.push_back(r1);
                    virus_reads.push_back(r2);
                    if (point_twd_bp(r2, r2_dir)) {
                        // virus read is pointing towards bp, so it must be rc in human
                        bam1_t *rc_copy = rc_read(bam_dup1(r2));
                        rc_copy->core.flag |= BAM_FSECONDARY;
                        host_reads.push_back(rc_copy);
                    } else {
                        // virus read pointing away from bp, no need to rc it
                        bam1_t *copy = bam_dup1(r2);
                        copy->core.flag |= BAM_FSECONDARY;
                        host_reads.push_back(copy);
                    }
                } else if (r1_is_host) { // both reads on host
                    if (!point_twd_bp(r1, r1_dir)) {
                        std::swap(r1, r2);
                    }
                    host_reads.push_back(r1);

                    bam1_t* r2_copy = bam_dup1(r2);
                    rc_read(r2_copy);
                    r2_copy->core.flag |= BAM_FSECONDARY;
                    host_reads.push_back(r2_copy);

                    auto vc_pos = good_clips_pos[r2->core.flag & BAM_FREAD1 ? 0 : 1][bam_get_qname(r2)];
                    r2->core.tid = vc_pos.first;
                    r2->core.pos = vc_pos.second;
                    virus_reads.push_back(r2);
                } else { // both reads on virus
                    // r1 must be the mate pointing away from bp
                    if (point_twd_bp(r1, r1_dir)) {
                        std::swap(r1, r2);
                    }

                    auto hc_pos = good_clips_pos[r1->core.flag & BAM_FREAD1 ? 0 : 1][bam_get_qname(r1)];
                    r1->core.tid = hc_pos.first;
                    r1->core.pos = hc_pos.second;
                    host_reads.push_back(r1);

                    bam1_t *r2_copy = bam_dup1(r2);
                    rc_read(r2_copy);
                    r2_copy->core.flag |= BAM_FSECONDARY;
                    host_reads.push_back(r2_copy);

                    virus_reads.push_back(r2);
                }
            } else { // destroy unitilised reads
                for (std::pair<bam1_t*, bool> p : e.second) {
                    bam_destroy1(p.first);
                }
            }
        }
    }

    bam_destroy1(read);


    std::unordered_set<std::string> host_qnames;
    google::dense_hash_map<std::string, bam1_t*> virus_reads_by_name;
    virus_reads_by_name.set_empty_key("");
    for (bam1_t* read : host_reads) {
        host_qnames.insert(bam_get_qname(read));
    }
    for (int i = 0; i < virus_reads.size(); i++) {
        std::string qname = bam_get_qname(virus_reads[i]);
        if (!host_qnames.count(qname)) {
            bam_destroy1(virus_reads[i]);
            virus_reads[i] = nullptr;
        } else {
            virus_reads_by_name[qname] = virus_reads[i];
        }
    }
    virus_reads.erase(std::remove(virus_reads.begin(), virus_reads.end(), nullptr), virus_reads.end());

    std::cerr << "HOST READS: " << host_reads.size() << std::endl;
    std::cerr << "VIRUS READS: " << virus_reads.size() << std::endl;

    /* ====== */


    /* === DEDUP READS === */

    dedup_reads(host_reads, virus_reads);

    /* ====== */


    /* === EXTRACT REGIONS === */

    std::string host_regions_path = workdir + "/host-regions";
    std::vector<region_t*> host_regions;
    std::ifstream host_regions_fin(host_regions_path);

    if (host_regions_fin.good()) {
        int region_id, start, end;
        std::string contig_name;
        while (host_regions_fin >> region_id >> contig_name >> start >> end) {
            int contig_id = contig_name2id[contig_name];
            region_t* region = new region_t(contig_id, contig_id2tid[contig_id], start, end);
            region->id = region_id;
            host_regions.push_back(region);
        }
    } else {
        extract_regions(host_reads, host_regions);

        // filter virus regions
        for (int i = 0; i < host_regions.size(); i++) {
            if (is_virus_region(host_regions[i])) {
                delete host_regions[i];
                host_regions[i] = nullptr;
            }
        }
        host_regions.erase(std::remove(host_regions.begin(), host_regions.end(), nullptr), host_regions.end());

        std::sort(host_regions.begin(), host_regions.end(), [](const region_t* reg1, const region_t* reg2) {
            if (reg1->contig_id == reg2->contig_id) return reg1->start < reg2->start;
            return reg1->contig_id < reg2->contig_id;
        });
        for (int i = 0; i < host_regions.size(); i++) {
            host_regions[i]->id = i;
        }

        std::ofstream host_regions_file(host_regions_path);
        for (region_t* region : host_regions) {
            host_regions_file << region->id << " " << contig_id2name[region->contig_id] << " " <<
            region->start << " " << region->end << std::endl;
        }
        host_regions_file.close();
    }
    host_regions_fin.close();

    std::cerr << "HOST REGIONS: " << host_regions.size() << std::endl;

    /* ====== */


    /* === COMPUTE READS-REGIONS ALIGNMENTS === */

    reads_per_region = new std::vector<read_realignment_t>[host_regions.size()];
    reads_per_region_rc = new std::vector<read_realignment_t>[host_regions.size()];
    cstring_to_cigar.set_empty_key("");

    std::string scores_file_path = workdir + "/reads.scores.bin";
    FILE* scores_file_in_bin = fopen(scores_file_path.c_str(), "rb");
    if (scores_file_in_bin) {
        std::cerr << "Reading qnames." << std::endl;
        std::vector<bam1_t*> id_to_read;
        load_qnames_indices(workdir, host_reads, id_to_read);

        std::cerr << "Reading cigars." << std::endl;
        load_cigars_indices(workdir, cid_to_cigar);

        std::cerr << "Reading mappings." << std::endl;
        char line[16];
        while (fread(line, 16, 1, scores_file_in_bin)) {
            uint32_t region_id = *((uint32_t *) (line + 0));
            uint32_t read_id = *((uint32_t *) (line + 4));
            bam1_t* read = id_to_read[read_id];
            uint32_t cigar_id = *((uint32_t *) (line + 8));
            bool is_rc = cigar_id & 0x80000000;
            std::pair<int, const uint32_t *>& cigar_v = cid_to_cigar[cigar_id & 0x7FFFFFFF];
            uint16_t offset_start = *((uint16_t *) (line + 12));
            uint16_t offset_end = offset_start + bam_cigar2rlen(cigar_v.first, cigar_v.second);
            uint16_t score = *((uint16_t *) (line + 14));
            read_realignment_t rr(read, offset_start, offset_end, cigar_v.first, cigar_v.second, score, is_rc);
            add_realignment_to_region(rr, region_id);
        }
        fclose(scores_file_in_bin);
        std::cerr << "Mappings read." << std::endl;
    } else {
        scores_file_out_bin = fopen(scores_file_path.c_str(), "wb");

        /* == INDEX REGIONS == */
        kmer_index = new std::vector<region_t *>[1 << KMER_BITS];
        kmer_index_rc = new std::vector<region_t *>[1 << KMER_BITS];
        mtx_kmers = new std::mutex[1 << KMER_BITS];

        ctpl::thread_pool thread_pool1(config.threads);
        std::vector<std::future<void> > futures;

        int regions_per_thread = host_regions.size() / config.threads;

        for (int i = 0; i < config.threads - 1; i++) {
            std::future<void> future = thread_pool1.push(index_regions, &host_regions, i * regions_per_thread,
                                                         (i + 1) * regions_per_thread);
            futures.push_back(std::move(future));
        }
        std::future<void> future = thread_pool1.push(index_regions, &host_regions,
                                                     (config.threads - 1) * regions_per_thread,
                                                     host_regions.size());
        futures.push_back(std::move(future));

        thread_pool1.stop(true);
        for (int i = 0; i < futures.size(); i++) {
            try {
                futures[i].get();
            } catch (char const *s) {
                std::cout << s << std::endl;
            }
        }
        delete[] mtx_kmers;
        /* ==== */

        /* == Associate reads to regions == */
        std::cerr << "READS: " << host_reads.size() << std::endl;

        // randomize reads to avoid imbalances in multi-threading
        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(host_reads), std::end(host_reads), rng);

        write_qnames_indices(workdir, host_reads);

        mtx_regions = new std::mutex[host_regions.size()];

        ctpl::thread_pool thread_pool2(config.threads);
        futures.clear();

        int reads_chunks = config.threads * 5;
        int reads_per_thread = host_reads.size() / reads_chunks;

        for (int i = 0; i < reads_chunks - 1; i++) {
            std::future<void> future = thread_pool2.push(associate_reads_to_regions, &host_reads,
                    i * reads_per_thread, (i + 1) * reads_per_thread);
            futures.push_back(std::move(future));
        }
        future = thread_pool2.push(associate_reads_to_regions, &host_reads,
                (reads_chunks - 1) * reads_per_thread, host_reads.size());
        futures.push_back(std::move(future));

        thread_pool2.stop(true);
        for (int i = 0; i < futures.size(); i++) {
            try {
                futures[i].get();
            } catch (char const *s) {
                std::cout << s << std::endl;
            }
        }

        delete[] mtx_regions;
        delete[] kmer_index;
        delete[] kmer_index_rc;

        // write file
        std::cerr << "WRITING FILE reads.score.bin" << std::endl;
        google::dense_hash_map<bam1_t*, uint32_t> read_ids;
        google::dense_hash_map<std::string, uint32_t> cigar_ids;
        read_ids.set_empty_key(NULL);
        cigar_ids.set_empty_key("");
        for (int i = 0; i < host_reads.size(); i++) {
            read_ids[host_reads[i]] = i;
        }
        for (int i = 0; i < host_regions.size(); i++) {
            auto& v_fwd = reads_per_region[i];
            for (read_realignment_t& rr : v_fwd) {
                char line[16]; // each line is 16 bytes
                std::string cigar = cigar_array_to_str(rr.cigar_len, rr.cigar);
                if (cigar_ids.count(cigar) == 0) cigar_ids[cigar] = cigar_ids.size();
                uint32_t cigar_id = cigar_ids[cigar];
                write_alignment_to_bytes(line, i, read_ids[rr.read], cigar_id, rr.offset_start, rr.rc, rr.score);
                fwrite(line, 1, 16, scores_file_out_bin);
            }

            auto& v_rev = reads_per_region_rc[i];
            for (read_realignment_t& rr : v_rev) {
                char line[16]; // each line is 16 bytes
                std::string cigar = cigar_array_to_str(rr.cigar_len, rr.cigar);
                if (cigar_ids.count(cigar) == 0) cigar_ids[cigar] = cigar_ids.size();
                uint32_t cigar_id = cigar_ids[cigar];
                write_alignment_to_bytes(line, i, read_ids[rr.read], cigar_id, rr.offset_start, rr.rc, rr.score);
                fwrite(line, 1, 16, scores_file_out_bin);
            }
        }

        std::ofstream cigar_map_out(workdir + "/cigars-map");
        for (auto& e : cigar_ids) {
            cigar_map_out << e.second << " " << e.first << std::endl;
        }
        cigar_map_out.close();

        // sort regions by number of matched reads, descending
        std::vector<int> sorted_regions_id, sorted_regions_id_rc;
        for (int i = 0; i < host_regions.size(); i++) {
            sorted_regions_id.push_back(i);
            sorted_regions_id_rc.push_back(i);
        }
        std::sort(sorted_regions_id.begin(), sorted_regions_id.end(), [](const int i1, const int i2) {
            return reads_per_region[i1].size() > reads_per_region[i2].size();
        });
        std::sort(sorted_regions_id_rc.begin(), sorted_regions_id_rc.end(), [](const int i1, const int i2) {
            return reads_per_region_rc[i1].size() > reads_per_region_rc[i2].size();
        });

        fclose(scores_file_out_bin);
        /* ==== */
    }

    /* ====== */


    // sort reads-region associations by score
    std::unordered_map<bam1_t*, region_t*> read_aln_count;
    google::dense_hash_set<std::string> not_uniquely_aligned;
    not_uniquely_aligned.set_empty_key("");
    for (int i = 0; i < host_regions.size(); i++) {
        auto rr_cmp = [](const read_realignment_t& ras1, const read_realignment_t& ras2) {
            return ras1.score > ras2.score;
        };
        std::sort(reads_per_region[i].begin(), reads_per_region[i].end(), rr_cmp);
        std::sort(reads_per_region_rc[i].begin(), reads_per_region_rc[i].end(), rr_cmp);

        region_t* region = host_regions[i];
        for (auto it = reads_per_region[i].begin(); it != reads_per_region[i].end(); it++) {
            if (!read_aln_count.count(it->read)) {
                read_aln_count[it->read] = region;
            } else {
                region_t* other_region = read_aln_count[it->read];
                if (!region->overlaps_with(other_region)) {
                    not_uniquely_aligned.insert(bam_get_qname(it->read));
                }
            }
        }
        for (auto it = reads_per_region_rc[i].begin(); it != reads_per_region_rc[i].end(); it++) {
            if (!read_aln_count.count(it->read)) {
                read_aln_count[it->read] = region;
            } else {
                region_t* other_region = read_aln_count[it->read];
                if (!region->overlaps_with(other_region)) {
                    not_uniquely_aligned.insert(bam_get_qname(it->read));
                }
            }
        }
    }

    // find uniquely aligned reads
    google::dense_hash_set<std::string> uniquely_aligned;
    uniquely_aligned.set_empty_key("");
    for (bam1_t* r : host_reads) {
        std::string qname = bam_get_qname(r);
        if (!not_uniquely_aligned.count(qname)) {
            uniquely_aligned.insert(qname);
        }
    }
    std::cerr << "UNIQUE READS: " << uniquely_aligned.size() << std::endl;


    /* === COMPUTE INTEGRATIONS === */

    google::dense_hash_set<std::string> already_used;
    already_used.set_empty_key("");

//    std::ofstream host_bp_seq_out(workdir + "/host_bp_seqs.fa");
    std::ofstream virus_bp_seq_out(workdir + "/virus_bp_seqs.fa");

    auto region_score_cmp = [](const region_score_t& r1, const region_score_t& r2) {
        if (r1.score == r2.score) return r1.id < r2.id;
        return r1.score < r2.score;
    };
    std::priority_queue<region_score_t, std::vector<region_score_t>, decltype(region_score_cmp)> region_scores(region_score_cmp);

    int region_score_id = 0;
    for (region_t* region : host_regions) {
        std::pair<int, int> score_and_reads = compute_region_score(region, false, already_used);
        int score = score_and_reads.first;
        int reads = score_and_reads.second;
        if (score > 0) {
            region_scores.push(region_score_t(region_score_id++, region, true, score, reads, reads));
        }

        score_and_reads = compute_region_score(region, true, already_used);
        score = score_and_reads.first;
        reads = score_and_reads.second;
        if (score > 0) {
            region_scores.push(region_score_t(region_score_id++, region, false, score, reads, reads));
        }
    }

    virus_reads_seqs.push_back(read_seq_t());


    std::vector<bam1_t*> remapped;
    while (!region_scores.empty()) {

        std::vector<read_realignment_t> kept_host_read_realignments, kept_virus_read_realignments;
        std::pair<region_t*, bool> best_virus_region = {NULL, false};
        google::dense_hash_map<std::string, std::vector<std::string> > dedup_qnames;
        while (!region_scores.empty()) {

            /* === VIRUS REALIGNMENT === */
            region_score_t best_region_score = region_scores.top();
            region_scores.pop();
            kept_host_read_realignments.clear();
            kept_virus_read_realignments.clear();

            std::cerr << "CURRENT TOP REGION: " << best_region_score.to_str() << std::endl;

            // remove already used read_and_score
            std::vector<read_realignment_t>& host_read_realignments = best_region_score.fwd ?
                    reads_per_region[best_region_score.region->id] : reads_per_region_rc[best_region_score.region->id];
            host_read_realignments.erase(
                    std::remove_if(host_read_realignments.begin(), host_read_realignments.end(), [&already_used](read_realignment_t& ras) {
                        return already_used.count(bam_get_qname(ras.read));
                    }),
                    host_read_realignments.end());
            best_region_score.reads = host_read_realignments.size();

            // remap virus reads and find best virus regions
            std::vector<read_realignment_t> virus_read_realignments;
            window_t host_window, virus_window;
            best_virus_region = remap_virus_reads(best_region_score, virus_read_realignments, virus_reads_by_name,
                    host_window, virus_window);


            if (best_virus_region.first == NULL) continue;

            region_scores.push(best_region_score);
            if (best_region_score.id != region_scores.top().id) continue;


            /* === DEDUP AFTER REALIGNMENT === */
            region_scores.pop();

            dedup_qnames = dedup_reads(host_read_realignments, virus_read_realignments, best_region_score.region,
                    best_region_score.fwd, best_virus_region.first, best_virus_region.second);


            // kept_host_read_realignments and kept_virus_read_realignments contain deduped realignments where both host
            // virus sides were remapped
            google::dense_hash_set<std::string> remapped_virus_qnames;
            remapped_virus_qnames.set_empty_key("");
            for (read_realignment_t& vrr : virus_read_realignments) {
                if (dedup_qnames.count(bam_get_qname(vrr.read))) {
                    remapped_virus_qnames.insert(bam_get_qname(vrr.read));
                    kept_virus_read_realignments.push_back(vrr);
                }
            }

            best_region_score.reads = 0;
            best_region_score.score = 0;
            for (read_realignment_t& hrr : host_read_realignments) {
                if (remapped_virus_qnames.count(bam_get_qname(hrr.read))) {
                    best_region_score.reads++;
                    best_region_score.score += hrr.score;
                    kept_host_read_realignments.push_back(hrr);
                }
            }
            region_scores.push(best_region_score);
            if (best_region_score.id != region_scores.top().id) continue;
            region_scores.pop();

            google::dense_hash_set<std::string> qnames_to_remove =
                    mismatch_filter(best_region_score.region, kept_host_read_realignments);
            auto to_remove = [&qnames_to_remove] (const read_realignment_t& rr) {
                return qnames_to_remove.count(bam_get_qname(rr.read));
            };

            kept_host_read_realignments.erase(
                    std::remove_if(kept_host_read_realignments.begin(), kept_host_read_realignments.end(), to_remove),
                    kept_host_read_realignments.end());
            kept_virus_read_realignments.erase(
                    std::remove_if(kept_virus_read_realignments.begin(), kept_virus_read_realignments.end(), to_remove),
                    kept_virus_read_realignments.end());

            // identify suspicious alignments - i.e. all low-complexity or clipped wrong because of garbage tail
            auto is_suspicious = [](read_realignment_t& rr, window_t accept_window) {
                if (is_low_complexity(rr.read, rr.rc, rr.left_clip_len(), rr.read->core.l_qseq-rr.right_clip_len())) { // mark low-complex as suspicious
                    return true;
                }
                if (rr.left_clip_len() > config.max_sc_dist && abs(rr.offset_start - accept_window.start) > config.max_sc_dist) { // mark wrongly clipped as suspicious
                    return true;
                }
                if (rr.right_clip_len() > config.max_sc_dist && abs(rr.offset_end - accept_window.end) > config.max_sc_dist) {
                    return true;
                }
                return false;
            };
            for (read_realignment_t& hrr : kept_host_read_realignments) {
                hrr.suspicious = is_suspicious(hrr, host_window);
            }
            for (read_realignment_t& vrr : kept_virus_read_realignments) {
                vrr.suspicious = is_suspicious(vrr, virus_window);
            }

            best_region_score.reads = 0;
            best_region_score.score = 0;
            for (read_realignment_t& hrr : kept_host_read_realignments) {
                if (!hrr.suspicious) {
                    best_region_score.score += hrr.score;
                    best_region_score.reads++;
                }
            }

            region_scores.push(best_region_score);
            if (best_region_score.id != region_scores.top().id) continue;

            std::cerr << "CONFIRMED TOP REGION: " << best_region_score.to_str() << std::endl;
            break;
        }

        if (region_scores.empty()) break;

        region_score_t best_region_score = region_scores.top();
        region_t* best_region = best_region_score.region;
        bool host_region_fwd = best_region_score.fwd;
        if (best_region_score.score <= 100) { // TODO: use read len instead
            std::cerr << "REMAINING REGIONS HAVE LOW SCORE." << std::endl;
            break;
        }

        for (read_realignment_t& rr : kept_host_read_realignments) {
            std::string qname = bam_get_qname(rr.read);
            already_used.insert(qname);
            for (std::string dup_qname : dedup_qnames[qname]) {
                already_used.insert(dup_qname);
            }
        }

        std::cerr << "HOST READS: " << kept_host_read_realignments.size() << std::endl;

        // edit virus reads and find bp
        edit_remapped_reads(best_virus_region.first, kept_virus_read_realignments, !best_virus_region.second);

        int v_min_pos = INT32_MAX, v_max_pos = 0;
        std::vector<double> virus_score_ratios;
        for (read_realignment_t& rr : kept_virus_read_realignments) {
            virus_score_ratios.push_back(rr.score/double(rr.offset_end-rr.offset_start+1));
            v_min_pos = std::min(v_min_pos, rr.read->core.pos);
            v_max_pos = std::max(v_max_pos, bam_endpos(rr.read));
        }
        int_breakpoint_t virus_bp(best_virus_region.first->contig_id, v_min_pos, v_max_pos, best_virus_region.second);

        // edit host reads and find host bp
        edit_remapped_reads(best_region, kept_host_read_realignments, !host_region_fwd);

        int reads_w_dups = 0, unique_reads_w_dups = 0;
        std::vector<double> host_score_ratios;
        for (read_realignment_t& rr : kept_host_read_realignments) {
            host_score_ratios.push_back(rr.score/double(rr.offset_end-rr.offset_start+1));
            reads_w_dups += dedup_qnames[bam_get_qname(rr.read)].size();
            if (uniquely_aligned.count(bam_get_qname(rr.read))) {
                unique_reads_w_dups += dedup_qnames[bam_get_qname(rr.read)].size();
            }
        }

        int score = 0, h_min_pos = INT32_MAX, h_max_pos = 0;
        for (read_realignment_t& rr : kept_host_read_realignments) {
            h_min_pos = std::min(h_min_pos, rr.read->core.pos);
            h_max_pos = std::max(h_max_pos, bam_endpos(rr.read));
            score += rr.score;
        }
        int_breakpoint_t host_bp(best_region->contig_id, h_min_pos, h_max_pos, host_region_fwd);


        // remove pairs having IS too long
        google::dense_hash_map<std::string, int> h_dists = get_distances(kept_host_read_realignments, host_bp, *best_region);
        google::dense_hash_map<std::string, int> v_dists = get_distances(kept_virus_read_realignments, virus_bp, *(best_virus_region.first));
        auto is_too_long = [&h_dists, &v_dists](const read_realignment_t& rr) {
            std::string qname = bam_get_qname(rr.read);
            return h_dists[qname] + v_dists[qname] > global_stats.max_is;
        };
        for (read_realignment_t& rr : kept_host_read_realignments) {
            if (is_too_long(rr)) {
                rr.too_long = true;
                bam_aux_update_str(rr.read, "CT", 2, "TL");
            }
            if (rr.suspicious) {
                bam_aux_update_str(rr.read, "CT", 2, "SS");
            }
        }
        for (read_realignment_t& rr : kept_virus_read_realignments) {
            if (is_too_long(rr)) {
                rr.too_long = true;
                bam_aux_update_str(rr.read, "CT", 2, "TL");
            }
            if (rr.suspicious) {
                bam_aux_update_str(rr.read, "CT", 2, "SS");
            }
        }

        best_region_score.good_reads = 0;
        best_region_score.score = 0;
        google::dense_hash_set<std::string> used_qnames;
        used_qnames.set_empty_key("");
        for (read_realignment_t& hrr : kept_host_read_realignments) {
            if (!hrr.suspicious && !hrr.too_long) {
                best_region_score.score += hrr.score;
                if (!used_qnames.count(bam_get_qname(hrr.read))) {
                    best_region_score.good_reads++;
                    used_qnames.insert(bam_get_qname(hrr.read));
                }
            }
        }

        double host_coverage = calculate_coverage(kept_host_read_realignments)/double(global_stats.max_is);
        double virus_coverage = calculate_coverage(kept_virus_read_realignments)/double(global_stats.max_is);

        int split_reads = 0;
        for (read_realignment_t& rr : kept_host_read_realignments) {
            if (host_bp.fwd && rr.right_clip_len() >= 30
            && abs((host_bp.max_pos-best_region_score.region->start)-rr.offset_end) < config.max_sc_dist) {
                split_reads++;
            } else if (!host_bp.fwd && rr.left_clip_len() >= 30
            && abs((host_bp.min_pos-best_region_score.region->start)-rr.offset_start) < config.max_sc_dist) {
                split_reads++;
            }
        }

        char out_region[10000];
        virus_bp_seq_out << ">" << virus_integration_id << "_" << (virus_bp.fwd ? "R" : "L") << std::endl;
        strncpy(out_region, chrs[contig_id2name[virus_bp.contig_id]].first+virus_bp.min_pos, virus_bp.max_pos-virus_bp.min_pos);
        out_region[virus_bp.max_pos-virus_bp.min_pos] = '\0';
        virus_bp_seq_out << out_region << std::endl;

        virus_integration_t v_int(best_region->contig_id, host_bp, virus_bp, best_region_score.reads,
                best_region_score.good_reads, split_reads,
                reads_w_dups, unique_reads_w_dups, best_region_score.score,
                mean(host_score_ratios), mean(virus_score_ratios), host_coverage, virus_coverage);
        std::cout << v_int.to_str() << std::endl;

        samFile* writer = open_bam_writer(workdir+"/readsx", std::to_string(v_int.id)+".bam", host_and_virus_file->header);
        for (read_realignment_t& rr : kept_virus_read_realignments) {
            // we duplicated viruses to deal with breakpoints across the border - scale back those alignment to original virus
            int virus_len = virus_names[contig_id2name[virus_bp.contig_id]];
            if (rr.read->core.pos >= virus_len) { // read past virus_len - simply scale back pos
                rr.read->core.pos %= virus_len;
            } else if (bam_endpos(rr.read) > virus_len) { // read split between start and end of virus
                // TODO
            }

            int ok = sam_write1(writer, host_and_virus_file->header, rr.read);
            if (ok < 0) throw "Failed to write to " + std::string(writer->fn);
        }
        for (read_realignment_t& rr : kept_host_read_realignments) {
            int ok = sam_write1(writer, host_and_virus_file->header, rr.read);
            if (ok < 0) throw "Failed to write to " + std::string(writer->fn);
        }
        sam_close(writer);

        std::cerr << "INTEGRATION ID: " << v_int.id << std::endl;
        std::cerr << "ALREADY USED: " << already_used.size() << std::endl << std::endl;
    }


    close_samFile(host_and_virus_file);
    for (bam1_t* r : host_reads) bam_destroy1(r);
    for (bam1_t* r : virus_reads) bam_destroy1(r);
    for (auto& e : chrs) {
        delete[] e.second.first;
    }
    for (region_t* region : host_regions) delete region;
    for (auto& e : cstring_to_cigar) {
        delete[] e.second.second;
    }
    for (auto& e : cid_to_cigar) {
        delete[] e.second;
    }
    delete[] reads_per_region;
    delete[] reads_per_region_rc;

//    host_bp_seq_out.close();
    virus_bp_seq_out.clear();

    /* ==== */

}
