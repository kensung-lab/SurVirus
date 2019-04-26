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
#include <cstring>
#include <algorithm>
#include <random>
#include <chrono>
#include <sys/time.h>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "cluster.h"
#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"
#include "alglib/statistics.h"

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

typedef unsigned long long ull;

int KMER_LEN = 13;
int KMER_BITS = KMER_LEN * 2;
ull KMER_MASK = (1ll << KMER_BITS)-1;

int POP_SIZE = 10000;

config_t config;
stats_t stats;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;
std::unordered_map<int, int> contig_id2tid, contig_tid2id;
std::unordered_map<std::string, std::pair<char*, size_t> > chrs;

std::unordered_map<std::string, uint32_t> cigar_ids;

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
        sprintf(buffer, "%s:%c%d", contig_id2name[contig_id].c_str(), fwd ? '+' : '-', bp());
        return buffer;
    }
};
struct virus_integration_t {

    int id;
    int_breakpoint_t host_bp, virus_bp;
    int reads, score;
    double p_value, mean;
    double h_score_ratio, v_score_ratio;

    virus_integration_t(int contig_id, int_breakpoint_t& host_bp, int_breakpoint_t& virus_bp, int reads,
            int score, double p_value, double mean, double h_score_ratio, double v_score_ratio) :
    id(virus_integration_id++), host_bp(host_bp), virus_bp(virus_bp), reads(reads), score(score), p_value(p_value), mean(mean),
            h_score_ratio(h_score_ratio), v_score_ratio(v_score_ratio) {}

    std::string to_str() {
        char buffer[4096];
        sprintf(buffer, "ID=%d %s %s %d %d %lf %lf %lf %lf", id, host_bp.to_str().c_str(), virus_bp.to_str().c_str(), reads,
                score, p_value, mean, h_score_ratio, v_score_ratio);
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
};

struct read_realignment_t {
    bam1_t* read;
    uint16_t offset_start, offset_end;
    uint8_t cigar_len;
    const uint32_t* cigar;
    uint16_t score;
    bool rc;

    read_realignment_t() : read(NULL), offset_start(0), offset_end(0), cigar_len(0), cigar(NULL), score(0), rc(false) {}
    read_realignment_t(bam1_t* read) : read(read), offset_start(0), offset_end(0), cigar_len(0), cigar(NULL), score(0), rc(false) {};
    read_realignment_t(bam1_t* read, uint16_t offset_start, uint16_t offset_end, uint8_t cigar_len, const uint32_t* cigar,
            uint16_t score, bool rc)
    : read(read), offset_start(offset_start), offset_end(offset_end), cigar_len(cigar_len), cigar(cigar), score(score), rc(rc) {}

    bool accepted() { return cigar != NULL; }

    bool is_left_clipped() {
        return bam_cigar_opchr(cigar[0]) == 'S';
    }
    bool is_right_clipped() {
        return bam_cigar_opchr(cigar[cigar_len-1]) == 'S';
    }
};

struct region_score_t {
    int id;
    region_t* region;
    bool fwd;
    int score, reads;

    region_score_t(int id, region_t* region, bool fwd, int score, int reads)
    : id(id), region(region), fwd(fwd), score(score), reads(reads) {};

    std::string to_str() {
        return std::to_string(id) + " " + region->to_str() + " READS=" + std::to_string(reads) + " SCORE=" + std::to_string(score);
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


char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
};

std::unordered_set<std::string> virus_names;
bool is_host_region(const region_t *region) {
    return virus_names.count(contig_id2name[region->contig_id]) == 0;
};
bool is_virus_region(const region_t *region) {
    return virus_names.count(contig_id2name[region->contig_id]);
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
    max_dists.push_back(stats.max_is);

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
                if (a1 != a2 && anchor_t::can_merge(*a1, *a2, stats.max_is) &&
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
                if (anchor_t::can_merge(*new_anchor, *map_it->second, stats.max_is)) {
                    pq.push(aa_distance_t(anchor_t::distance(*new_anchor, *map_it->second), new_anchor,
                                          map_it->second));
                }
            }
            anchors_map.insert({new_anchor->start, new_anchor});
            anchors_map.insert({new_anchor->end, new_anchor});
        }

        anchors.erase(std::remove_if(anchors.begin(), anchors.end(), [] (const anchor_t* a) {return a->dead;}), anchors.end());
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
                                           std::max(0, a->start - stats.max_is + a->len()),
                                           std::min(a->end + stats.max_is - a->len(),
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


bool accept_alignment(StripedSmithWaterman::Alignment &alignment, std::string &query) {
    return !is_poly_ACGT(query.c_str()+alignment.query_begin, alignment.query_end-alignment.query_begin+1)
    && alignment.sw_score >= 30; //read->core.l_qseq;
};


void write_alignment_to_bytes(char bytes[16], uint32_t region_id, uint32_t read_id, std::string& cigar, uint16_t ref_begin,
        bool is_rc, uint16_t score) {

    mtx.lock();
    if (cigar_ids.count(cigar) == 0) cigar_ids[cigar] = cigar_ids.size();
    uint32_t cigar_id = cigar_ids[cigar];
    mtx.unlock();
    cigar_id = (is_rc ? 0x80000000 : 0) | (cigar_id & 0x7FFFFFFF); // forcing is_rc inside the cigar_id to make 16 bytes

    memcpy(bytes, &region_id, 4);
    memcpy(bytes+4, &read_id, 4);
    memcpy(bytes+8, &cigar_id, 4);
    memcpy(bytes+12, &ref_begin, 2);
    memcpy(bytes+14, &score, 2);
}


void add_realignment_to_region(read_realignment_t& rr, int region_id) {
    if (rr.read->core.flag & BAM_FSECONDARY) { // good clip
        if (rr.is_left_clipped()) {
            reads_per_region_rc[region_id].push_back(rr);
        }
        if (rr.is_right_clipped()) {
            reads_per_region[region_id].push_back(rr);
        }
    } else { // the read points toward the bp
        std::vector<read_realignment_t>& reads_per_region_v = rr.rc ? reads_per_region_rc[region_id] : reads_per_region[region_id];
        reads_per_region_v.push_back(rr);
    }
}

inline void remap_read(std::string seq, bam1_t* read, int read_id, region_t* region, bool is_rc, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
                  region->end - region->start, filter, &alignment, 0);

    if (!accept_alignment(alignment, seq)) return;

    std::string cigar = alignment_cigar_to_bam_cigar(alignment.cigar);
    auto cigar_v = cigar_str_to_array(cigar);

    mtx_regions[region->id].lock();
    read_realignment_t rr(read, alignment.ref_begin, alignment.ref_end, cigar_v.first, cigar_v.second, alignment.sw_score, is_rc);
    add_realignment_to_region(rr, region->id);
    mtx_regions[region->id].unlock();

    char line[16]; // each line is 16 bytes
    write_alignment_to_bytes(line, region->id, read_id, cigar, alignment.ref_begin, is_rc, alignment.sw_score);

    mtx.lock();
    fwrite(line, 1, 16, scores_file_out_bin);
    mtx.unlock();


    remappings++;
    if (remappings % 1000000 == 0) {
        std::cerr << remappings << " remappings done." << std::endl;
    }
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

    StripedSmithWaterman::Aligner aligner(1, 2, 4, 1, false);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    for (int r = start+1; r <= end; r++) {
        bam1_t* read = (*host_reads)[r-1];

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

                        remap_read(seq, read, r-1, region, false, aligner, filter);
                    } else if (abs(last_read_for_region[region->id]) != r) {
                        last_read_for_region[region->id] = r;
                        last_pos[region->id] = i;
                    }
                }

                for (region_t* region : kmer_index_rc[kmer]) {
                    if (last_read_for_region_rc[region->id] == r && i - last_pos_rc[region->id] >= KMER_LEN) {
                        last_read_for_region_rc[region->id] = -r;

                        remap_read(seq_rc, read, r-1, region, true, aligner, filter);
                    } else if (abs(last_read_for_region_rc[region->id]) != r) {
                        last_read_for_region_rc[region->id] = r;
                        last_pos_rc[region->id] = i;
                    }
                }
            }
        }
    }

    delete[] last_read_for_region;
    delete[] last_pos;
    delete[] last_read_for_region_rc;
    delete[] last_pos_rc;
}


std::pair<int, int> compute_region_score(region_t* region, bool rc, std::unordered_set<bam1_t*>& already_used) {
    int score = 0, reads = 0;
    std::vector<read_realignment_t>& read_realignments = rc ? reads_per_region_rc[region->id] : reads_per_region[region->id];
    for (read_realignment_t& rr : read_realignments) {
        if (already_used.count(rr.read) == 0) {
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

//        if (shared_qnames.count(bam_get_qname(read))) {
//            if (rc) set_to_reverse(read);
//            else set_to_forward(read);
//        } else {
            if (rr.rc) set_to_reverse(read);
            else set_to_forward(read);
//        }

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
void load_cigars_indices(std::string& workspace, std::unordered_map<uint32_t, std::pair<int, const uint32_t*> >& id_to_cigar) {
    std::ifstream cigars_map_in(workspace + "/cigars-map");
    int i; std::string cigar_str;
    while (cigars_map_in >> i >> cigar_str) {
        id_to_cigar[i] = cigar_str_to_array(cigar_str);
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
std::pair<int, int> get_accept_window(std::vector<read_realignment_t>& realigned_reads, bool rc, int region_len,
        google::dense_hash_set<std::string>& deduped_qnames) {

    int clip_positions[100000];
    std::fill(clip_positions, clip_positions+region_len, 0);
    for (read_realignment_t& rr : realigned_reads) {
        if (!deduped_qnames.count(bam_get_qname(rr.read))) continue;
        if (rc && rr.is_left_clipped()) {
            clip_positions[rr.offset_start]++;
        }
        if (!rc && rr.is_right_clipped()) {
            clip_positions[rr.offset_end]++;
        }
    }

    // if it is FWD, we want to get the last maximum element
    int bp_pos = rc ? get_first_max_element(clip_positions, region_len) : get_last_max_element(clip_positions, region_len);

    if (clip_positions[bp_pos] < 2) { // cannot get bp from clips, find point of max score density
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
            for (int i = region_len-1-stats.max_is; i >= 0; i--) {
                score_windows[i] -= score_points[i+stats.max_is];
            }
        } else {
            for (int i = stats.max_is; i < region_len; i++) {
                score_windows[i] -= score_points[i-stats.max_is];
            }
        }

        bp_pos = rc ? get_last_max_element(score_windows, region_len) : get_first_max_element(score_windows, region_len);
    }

    int accept_window_start, accept_window_end;
    if (rc) {
        accept_window_start = bp_pos;
        accept_window_end = accept_window_start + stats.max_is;
    } else {
        accept_window_end = bp_pos;
        accept_window_start = accept_window_end - stats.max_is;
    }
    return {accept_window_start, accept_window_end};
}


void dedup_reads(std::vector<bam1_t*>& host_reads, std::vector<bam1_t*>& virus_reads) {

    google::dense_hash_map<std::string, bam1_t*> virus_reads_by_name;
    virus_reads_by_name.set_empty_key("");
    for (bam1_t* read : virus_reads) {
        virus_reads_by_name[bam_get_qname(read)] = read;
    }

    // DEDUP BY SEQUENCE

    // dedup pairs
    google::dense_hash_map<std::string, std::vector<std::pair<bam1_t*, bam1_t*> > > pairs_by_seq;
    pairs_by_seq.set_empty_key("");
    for (bam1_t* hread : host_reads) {
        if (virus_reads_by_name.count(bam_get_qname(hread))) {
            bam1_t* vread = virus_reads_by_name[bam_get_qname(hread)];
            std::string joint_string = get_sequence(hread, true) + "-" + get_sequence(vread, true);
            pairs_by_seq[joint_string].push_back({hread, vread});
        }
    }
    for (auto& e : pairs_by_seq) {
        std::sort(e.second.begin(), e.second.end(), [](const std::pair<bam1_t*, bam1_t*> p1, const std::pair<bam1_t*, bam1_t*> p2) {
            return strcmp(bam_get_qname(p1.first), bam_get_qname(p2.first)) < 0;
        });
    }

    // dedup host/virus anchors
    google::dense_hash_map<std::string, std::vector<bam1_t*> > host_anchors_by_seq;
    host_anchors_by_seq.set_empty_key("");
    for (bam1_t* hread : host_reads) {
        if (!virus_reads_by_name.count(bam_get_qname(hread))) {
            std::string seq = get_sequence(hread, true);
            host_anchors_by_seq[seq].push_back(hread);
        }
    }
    for (auto& e : host_anchors_by_seq) {
        std::sort(e.second.begin(), e.second.end(), [](const bam1_t* r1, const bam1_t* r2) {
            return strcmp(bam_get_qname(r1), bam_get_qname(r2)) < 0;
        });
    }

    // clear reads and insert only one for sequence
    host_reads.clear();
    virus_reads.clear();
    for (auto& e : pairs_by_seq) {
        host_reads.push_back(e.second[0].first);
        virus_reads.push_back(e.second[0].second);
    }
    for (auto& e : host_anchors_by_seq) {
        host_reads.push_back(e.second[0]);
    }


    std::cerr << "AFTER DEDUP BY SEQ" << std::endl;
    std::cerr << "HOST READS " << host_reads.size() << std::endl;
    std::cerr << "VIRUS READS " << virus_reads.size() << std::endl;


    // DEDUP BY POSITION AND CIGAR

    // dedup pairs
    std::vector<std::pair<bam1_t*, bam1_t*> > paired_reads;
    std::vector<bam1_t*> host_anchors;
    for (bam1_t* hread : host_reads) {
        if (virus_reads_by_name.count(bam_get_qname(hread))) {
            bam1_t* vread = virus_reads_by_name[bam_get_qname(hread)];
            paired_reads.push_back({hread, vread});
        } else {
            host_anchors.push_back(hread);
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

    host_reads.clear();
    virus_reads.clear();

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
            std::vector<std::pair<bam1_t*, bam1_t*> >& pairs_by_curr_seq = pairs_by_seq[curr_seq];
            std::vector<std::pair<bam1_t*, bam1_t*> >& pairs_by_last_valid_seq = pairs_by_seq[last_valid_seq];
            pairs_by_last_valid_seq.insert(pairs_by_last_valid_seq.end(), pairs_by_curr_seq.begin(), pairs_by_curr_seq.end());
        } else {
            host_reads.push_back(h1);
            virus_reads.push_back(v1);
            last_valid_i = i;
        }
    }
    host_reads.insert(host_reads.end(), host_anchors.begin(), host_anchors.end());

    std::cerr << "AFTER DEDUP BY POS" << std::endl;
    std::cerr << "HOST READS " << host_reads.size() << std::endl;
    std::cerr << "VIRUS READS " << virus_reads.size() << std::endl;
}

google::dense_hash_set<std::string> dedup_reads(std::vector<read_realignment_t>& host_read_realignments,
                                                std::vector<read_realignment_t> virus_read_realignments) {
    if (host_read_realignments.empty() || virus_read_realignments.empty()) {
        google::dense_hash_set<std::string> empty_set;
        empty_set.set_empty_key("");
        return empty_set;
    }

    std::unordered_map<std::string, read_realignment_t> host_reads_by_name;
    for (read_realignment_t& rr : host_read_realignments) {
        host_reads_by_name[bam_get_qname(rr.read)] = rr;
    }

    std::sort(virus_read_realignments.begin(), virus_read_realignments.end(),
              [] (const read_realignment_t& rr1, const read_realignment_t& rr2) {
                  auto rr1_tie = std::tie(rr1.offset_start, rr1.offset_end, rr1.cigar_len);
                  auto rr2_tie = std::tie(rr2.offset_start, rr2.offset_end, rr2.cigar_len);
                  if (rr1_tie != rr2_tie) return rr1_tie < rr2_tie;

                  // compare cigar itself
                  for (int i = 0; i < rr1.cigar_len; i++) {
                      if (rr1.cigar[i] != rr2.cigar[i]) return rr1.cigar[i] < rr2.cigar[i];
                  }
                  return false;
              });

    auto cigar_eq = [] (const read_realignment_t& rr1, const read_realignment_t& rr2) {
        if (rr1.cigar_len != rr2.cigar_len) return false;
        for (int i = 0; i < rr1.cigar_len; i++) if (rr1.cigar[i] != rr2.cigar[i]) return false;
        return true;
    };


    google::dense_hash_set<std::string> qnames_to_keep;
    qnames_to_keep.set_empty_key("");
    for (int i = 0; i < virus_read_realignments.size(); i++) {
        read_realignment_t& vrr_i = virus_read_realignments[i];
        read_realignment_t& hrr_i = host_reads_by_name[bam_get_qname(vrr_i.read)];

        bool has_prev_dup = false;
        for (int j = i-1; j >= 0; j--) {
            read_realignment_t& vrr_j = virus_read_realignments[j];
            if (vrr_i.offset_start != vrr_j.offset_start || !cigar_eq(vrr_i, vrr_j)) break;

            read_realignment_t& hrr_j = host_reads_by_name[bam_get_qname(vrr_j.read)];
            if (hrr_i.offset_start == hrr_j.offset_start && cigar_eq(hrr_i, hrr_j)) {
                has_prev_dup = true;
                break;
            }
        }

        if (!has_prev_dup) {
            qnames_to_keep.insert(bam_get_qname(vrr_i.read));
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
        StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {

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

        read_realignment_t& ras = cached_region_alignments[read_seq_id];
        if (!ras.accepted()) continue;

        host_read_realignments_local.push_back(host_alignemnts[i]);
        virus_read_realignments_local.push_back(ras);
    }


    google::dense_hash_set<std::string> deduped_qnames = dedup_reads(host_read_realignments_local, virus_read_realignments_local);

    google::dense_hash_set<std::string> outside_window;
    outside_window.set_empty_key("");
    std::pair<int,int> host_accept_window = get_accept_window(host_read_realignments_local, hregion_rc, host_region->len(), deduped_qnames);
    std::pair<int,int> virus_accept_window = get_accept_window(virus_read_realignments_local, vregion_rc, virus_region->len(), deduped_qnames);
    for (read_realignment_t& rr : host_read_realignments_local) {
        if (rr.offset_start < host_accept_window.first || rr.offset_end > host_accept_window.second) {
            outside_window.insert(bam_get_qname(rr.read));
        }
    }
    for (read_realignment_t& rr : virus_read_realignments_local) {
        if (rr.offset_start < virus_accept_window.first || rr.offset_end > virus_accept_window.second) {
            outside_window.insert(bam_get_qname(rr.read));
        }
    }

    int h_score = 0, v_score = 0;
    for (int i = 0; i < virus_read_realignments_local.size(); i++) {
        read_realignment_t& hrr = host_read_realignments_local[i];
        read_realignment_t& vrr = virus_read_realignments_local[i];
        if (outside_window.count(bam_get_qname(hrr.read))) continue;

        if (hrr.offset_start >= host_accept_window.first && hrr.offset_end <= host_accept_window.second
         && vrr.offset_start >= virus_accept_window.first && vrr.offset_end <= virus_accept_window.second) {
            if (deduped_qnames.count(bam_get_qname(hrr.read))) {
                h_score += hrr.score;
                v_score += vrr.score;
            }
            if (virus_read_realignments != NULL) {
                virus_read_realignments->push_back(vrr);
            }
        }
    }
    return {v_score, h_score};
}

std::pair<region_t*, bool> remap_virus_reads(region_score_t& host_region_score, std::vector<read_realignment_t>& virus_read_realignments,
        google::dense_hash_map<std::string, bam1_t*>& virus_reads_by_name, std::unordered_set<bam1_t*>& already_used) {
    StripedSmithWaterman::Aligner aligner(1, 2, 4, 1, false);
    StripedSmithWaterman::Filter filter;

    std::vector<read_realignment_t>& host_read_realignments = host_region_score.fwd ?
            reads_per_region[host_region_score.region->id] : reads_per_region_rc[host_region_score.region->id];
    // remove already used read_and_score
    host_read_realignments.erase(
            std::remove_if(host_read_realignments.begin(), host_read_realignments.end(), [&already_used](read_realignment_t& ras) {
                return already_used.count(ras.read);
            }),
            host_read_realignments.end());
    host_region_score.reads = host_read_realignments.size();

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

        if ((host_region_score.fwd && rr.is_right_clipped())
        || (!host_region_score.fwd && rr.is_left_clipped())) {
            remapped_host_anchors.push_back(rr);
        }
    }

    // extract virus regions
    std::vector<region_t*> virus_regions;
    extract_regions(remapped_reads_mates, virus_regions);
    virus_regions.erase(
            std::remove_if(virus_regions.begin(), virus_regions.end(), is_host_region), virus_regions.end()
    );
    for (region_t* region : virus_regions) {
        region->start = (region->start/APPROX_FACTOR) * APPROX_FACTOR;
        region->end = (region->end/APPROX_FACTOR) * APPROX_FACTOR + APPROX_FACTOR;
        region->end = std::min(region->end, (int) chrs[contig_id2name[region->contig_id]].second);
    }

    if (virus_regions.empty()) return {nullptr, false};

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
            if (host_region_score.fwd && anchor.is_right_clipped() && !anchor.rc) {
                set_to_forward(rc_read);
                rc_read->core.flag |= BAM_FREVERSE;
            }
            else if (!host_region_score.fwd && anchor.is_left_clipped() && anchor.rc) {
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
        std::pair<int, int> score = remap_virus_reads_supp(virus_region, false, read_seq_ids, host_region_score.region, !host_region_score.fwd,
                host_alignments, NULL, aligner, filter);
        std::pair<int, int> rc_score = remap_virus_reads_supp(virus_region, true, read_seq_ids, host_region_score.region, !host_region_score.fwd,
                host_alignments, NULL, aligner, filter);

        std::pair<int, int> max_score = std::max(score, rc_score);
        if (best_score < max_score) {
            best_score = max_score;
            best_region = virus_region;
            is_best_vregion_fwd = (max_score == score);
        }
    }

    // FIXME: what if best_region is null?
    std::cerr << "BEST VIRUS REGION: " << best_region->to_str() << " " << best_score.first << std::endl;
    remap_virus_reads_supp(best_region, !is_best_vregion_fwd, read_seq_ids, host_region_score.region, !host_region_score.fwd,
            host_alignments, &virus_read_realignments, aligner, filter);

    google::dense_hash_set<std::string> deduped_reads = dedup_reads(host_alignments, virus_read_realignments);

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
std::vector<double> gen_population_for_ins(std::string bam_fname, std::string workdir) {
    open_samFile_t* open_sam = open_samFile(bam_fname.c_str());

    std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
    std::vector<double> population;
    bam1_t* read = bam_init1();
    while (population.size() < POP_SIZE) {
        char contig[1000]; int pos;
        rnd_pos_fin >> contig >> pos;

        std::unordered_set<std::pair<int, int>, pair_hash> used_pairs;

        char region[1000];
        sprintf(region, "%s:%d-%d", contig, pos-2*stats.max_is, pos);
        hts_itr_t* iter = sam_itr_querys(open_sam->idx, open_sam->header, region);
        while (sam_itr_next(open_sam->file, iter, read) >= 0) {
            if (!is_valid(read, false) || bam_is_rev(read)) continue;

            if (is_inward(read, stats.min_is) && read->core.isize > 0 && read->core.isize >= stats.min_is
            && read->core.isize <= 2*stats.max_is) {
                if (used_pairs.count({read->core.pos, read->core.mpos})) continue; // we already sampled a duplicate pair

                int start = read->core.pos+read->core.l_qseq/2, end = get_mate_endpos(read)-read->core.l_qseq/2;
                if (start >= end) continue;

                if (start <= pos && pos <= end) {
                    used_pairs.insert({read->core.pos, read->core.mpos});
                    population.push_back((double) read->core.isize);
                }
            }
        }
        sam_itr_destroy(iter);
    }
    close_samFile(open_sam);
    bam_destroy1(read);

    return population;
}

int dist_from_bp(const read_realignment_t& rr, int_breakpoint_t& bp) {
    int read_endpos = bp.fwd ? rr.read->core.pos : bam_endpos(rr.read);
    return abs(read_endpos - bp.bp());
};
google::dense_hash_map<std::string, int> get_distances(std::vector<read_realignment_t>& realignments, int_breakpoint_t& bp) {
    google::dense_hash_map<std::string, int> dists;
    dists.set_empty_key("");
    for (read_realignment_t& rr : realignments) {
        std::string qname = bam_get_qname(rr.read);
        int dist = dist_from_bp(rr, bp);
        if (!dists.count(qname) || dists[qname] < dist) {
            dists[qname] = dist;
        }
    }
    return dists;
}
std::vector<std::pair<double, std::string> > get_dists_for_KStest(std::vector<read_realignment_t>& host_deduped_realignments,
                                         int_breakpoint_t& host_bp,
                                         std::vector<read_realignment_t>& virus_deduped_realignments,
                                         int_breakpoint_t& virus_bp) {

    // in case there are two reads with the same qname (i.e. a mate and an anchor/clip), choose the farthest from the bp
    google::dense_hash_map<std::string, int> h_dists = get_distances(host_deduped_realignments, host_bp);
    google::dense_hash_map<std::string, int> v_dists = get_distances(virus_deduped_realignments, virus_bp);

    std::vector<std::pair<double, std::string> > dists;
    for (auto e : h_dists) {
        dists.push_back({e.second + v_dists[e.first], e.first});
    }
    return dists;
}


int main(int argc, char* argv[]) {

    nucl_bm['A'] = nucl_bm['a'] = 0;
    nucl_bm['C'] = nucl_bm['c'] = 1;
    nucl_bm['G'] = nucl_bm['g'] = 2;
    nucl_bm['T'] = nucl_bm['t'] = 3;
    nucl_bm['N'] = nucl_bm['n'] = 0;

    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';

    std::string bam_fname = argv[1];
    std::string host_reference_fname  = argv[2];
    std::string virus_reference_fname  = argv[3];
    std::string workdir = argv[4];
    std::string workspace = argv[5];
    int forced_region_score_id = -1;
    if (argc > 6) {
        forced_region_score_id = std::stoi(argv[6]);
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
        virus_names.insert(seq->name.s);
        chrs[std::string(seq->name.s)] = {new char[seq->seq.l+1], seq->seq.l};
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
    }
    kseq_destroy(seq);
    fclose(fasta);

    config = parse_config(workdir + "/config.txt");
    stats = parse_stats(workspace + "/stats.txt");

    open_samFile_t* host_and_virus_file = open_samFile((workspace + "/retained-pairs-remapped.sorted.bam").data()); // just for the header
    for (int i = 0; i < host_and_virus_file->header->n_targets; i++) {
        int contig_id = contig_name2id[host_and_virus_file->header->target_name[i]];
        contig_id2tid[contig_id] = i;
        contig_tid2id[i] = contig_id;
    }


    /* === READ CHIMERIC PAIRS/READS === */

    /* == Pouring reads which are not good clips into host_reads and virus_reads == */

    bam1_t* read = bam_init1();

    // extracting paired host-virus reads
    std::vector<bam1_t*> host_reads;
    std::unordered_set<std::string> host_qnames;
    open_samFile_t* host_file = open_samFile((workspace + "/host-side.sorted.bam").data(), true);
    while (sam_read1(host_file->file, host_file->header, read) >= 0) {
        host_reads.push_back(bam_dup1(read));
        host_qnames.insert(bam_get_qname(read));
    }

    std::vector<bam1_t*> virus_reads;
    std::unordered_set<std::string> virus_qnames;
    open_samFile_t* virus_file = open_samFile((workspace + "/virus-side.sorted.bam").data(), true);
    while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
        virus_reads.push_back(bam_dup1(read));
        virus_qnames.insert(bam_get_qname(read));
    }
    close_samFile(virus_file);

    /* == */

    std::unordered_map<std::string, std::vector<std::pair<bam1_t*, bool> > > good_clipped_pairs;
    // pairs where both are good clipped reads. Bool = true indicates a host read

    /* == Reading host anchors == */

    open_samFile_t* host_anchors_file = open_samFile((workspace + "/host-anchors.sorted.bam").data(), true);
    while (sam_read1(host_anchors_file->file, host_anchors_file->header, read) >= 0) {
        if (host_qnames.count(bam_get_qname(read))) {
            virus_reads.push_back(bam_dup1(read));
            bam1_t* copy = bam_dup1(read);
            copy->core.flag |= BAM_FSECONDARY;
            host_reads.push_back(copy);
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
        if (host_qnames.count(bam_get_qname(read))) {
            virus_reads.push_back(bam_dup1(read));
            bam1_t* rc_copy = rc_read(bam_dup1(read));
            rc_copy->core.flag |= BAM_FSECONDARY;
            host_reads.push_back(rc_copy);
        } else if (virus_qnames.count(bam_get_qname(read))) {
            host_reads.push_back(bam_dup1(read));
        } else {
            good_clipped_pairs[bam_get_qname(read)].push_back({bam_dup1(read), false});
        }
    }
    close_samFile(virus_anchors_file);

    /* == */

    std::cerr << "HOST READS: " << host_reads.size() << std::endl;

    for (auto& e : good_clipped_pairs) {
        if (e.second.size() == 2) {
            bam1_t* r1 = e.second[0].first,* r2 = e.second[1].first;
            bool r1_is_host = e.second[0].second, r2_is_host = e.second[1].second;
            char r1_dir = bam_aux2A(bam_aux_get(r1, "CD")), r2_dir = bam_aux2A(bam_aux_get(r2, "CD"));

            // r1 will be the host read, r2 the virus read and the secondary host read

            if (r1_is_host != r2_is_host) {
                // if one is from host and one from virus, then it's easy
                host_reads.push_back(r1);
                virus_reads.push_back(r2);
                if ((bam_is_rev(r2) && r2_dir == 'L') || (!bam_is_rev(r2) && r2_dir == 'R')) {
                    // virus read is pointing towards bp, so it must be rc in human
                    bam1_t* rc_copy = rc_read(bam_dup1(r2));
                    rc_copy->core.flag |= BAM_FSECONDARY;
                    host_reads.push_back(rc_copy);
                } else {
                    // virus read pointing away from bp, no need to rc it
                    bam1_t* copy = bam_dup1(r2);
                    copy->core.flag |= BAM_FSECONDARY;
                    host_reads.push_back(copy);
                }
            } else if (r1_is_host) { // both reads on host
                host_reads.push_back(r1);
                bam1_t* r2_copy = bam_dup1(r2);
                r2_copy->core.flag |= BAM_FSECONDARY;
                host_reads.push_back(r2);

                r2_copy = bam_dup1(r2);
                if ((bam_is_rev(r2) && r2_dir == 'L') || (!bam_is_rev(r2) && r2_dir == 'R')) {
                    // r2 is pointing towards bp, so it must be rc in virus
                    rc_read(r2_copy);
                }
                virus_reads.push_back(r2_copy);
            } else { // both reads on virus
                if ((bam_is_rev(r1) && r1_dir == 'L') || (!bam_is_rev(r1) && r1_dir == 'R')) {
                    // r2 is pointing towards bp, so it must be rc in virus
                    rc_read(r1);
                }
                host_reads.push_back(r1);

                virus_reads.push_back(r2);
                bam1_t* r2_copy = bam_dup1(r2);
                if ((bam_is_rev(r2) && r2_dir == 'L') || (!bam_is_rev(r2) && r2_dir == 'R')) {
                    // r2 is pointing towards bp, so it must be rc in virus
                    rc_read(r2_copy);
                }
                r2_copy->core.flag |= BAM_FSECONDARY;
                host_reads.push_back(r2_copy);
            }
        }
    }

    google::dense_hash_map<std::string, bam1_t*> virus_reads_by_name;
    virus_reads_by_name.set_empty_key("");
    for (bam1_t* read : virus_reads) {
        virus_reads_by_name[bam_get_qname(read)] = read;
    }

    std::cerr << "HOST READS: " << host_reads.size() << std::endl;

    /* ====== */


    /* === DEDUP READS === */

    dedup_reads(host_reads, virus_reads);

    /* ====== */


    /* === EXTRACT REGIONS === */

    std::string host_regions_path = workspace + "/host-regions";
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
        host_regions.erase(
                std::remove_if(host_regions.begin(), host_regions.end(), is_virus_region),
                host_regions.end());

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

    std::string scores_file_path = workspace + "/reads.scores.bin";
    FILE* scores_file_in_bin = fopen(scores_file_path.c_str(), "rb");
    if (scores_file_in_bin) {
        std::cerr << "Reading qnames." << std::endl;
        std::vector<bam1_t*> id_to_read;
//        std::unordered_map<uint32_t, bam1_t*> id_to_read;
        load_qnames_indices(workspace, host_reads, id_to_read);

        std::cerr << "Reading cigars." << std::endl;
//        std::vector<std::pair<int, const uint32_t*> > id_to_cigar; TODO: if use vector version (faster) google perftools crashes
        std::unordered_map<uint32_t, std::pair<int, const uint32_t*> > id_to_cigar;
        load_cigars_indices(workspace, id_to_cigar);

        std::cerr << "Reading mappings." << std::endl;
        char line[16];
        while (fread(line, 16, 1, scores_file_in_bin)) {
            uint32_t region_id = *((uint32_t *) (line + 0));
            uint32_t read_id = *((uint32_t *) (line + 4));
            bam1_t* read = id_to_read[read_id];
            uint32_t cigar_id = *((uint32_t *) (line + 8));
            bool is_rc = cigar_id & 0x80000000;
            std::pair<int, const uint32_t *>& cigar_v = id_to_cigar[cigar_id & 0x7FFFFFFF];
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
        /* ==== */

        /* == Associate reads to regions == */
        std::cerr << "READS: " << host_reads.size() << std::endl;

        // randomize reads to avoid imbalances in multi-threading
        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(host_reads), std::end(host_reads), rng);

        write_qnames_indices(workspace, host_reads);

        mtx_regions = new std::mutex[host_regions.size()];

        ctpl::thread_pool thread_pool2(config.threads);
        futures.clear();

        int reads_chunks = config.threads * 5;
        int reads_per_thread = host_reads.size() / reads_chunks;

        for (int i = 0; i < reads_chunks - 1; i++) {
            std::future<void> future = thread_pool2.push(associate_reads_to_regions, &host_reads,// &region_trees,
                    i * reads_per_thread, (i + 1) * reads_per_thread);
            futures.push_back(std::move(future));
        }
        future = thread_pool2.push(associate_reads_to_regions, &host_reads,// &region_trees,
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

        std::ofstream cigar_map_out(workspace + "/cigars-map");
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
    for (int i = 0; i < host_regions.size(); i++) {
        auto rr_cmp = [](const read_realignment_t& ras1, const read_realignment_t& ras2) {
            return ras1.score > ras2.score;
        };
        std::sort(reads_per_region[i].begin(), reads_per_region[i].end(), rr_cmp);
        std::sort(reads_per_region_rc[i].begin(), reads_per_region_rc[i].end(), rr_cmp);
    }


    /* === GENERATE POPULATION FOR KS-TEST === */

    std::ofstream pop_outf(workdir + "/population.txt");
    std::vector<double> population = gen_population_for_ins(bam_fname, workdir);
    for (double d : population) {
        pop_outf << d << std::endl;
    }
    pop_outf.close();

    alglib::real_1d_array pop_v;
    pop_v.setcontent(population.size(), population.data());
    std::cerr << "POPULATION MEAN: " << mean(population) << std::endl;

    /* ==== */


    /* === COMPUTE INTEGRATIONS === */

    std::unordered_set<bam1_t*> already_used;

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
        if (score > 0) region_scores.push(region_score_t(region_score_id++, region, true, score, reads));

        score_and_reads = compute_region_score(region, true, already_used);
        score = score_and_reads.first;
        reads = score_and_reads.second;
        if (score > 0) region_scores.push(region_score_t(region_score_id++, region, false, score, reads));
    }


    virus_reads_seqs.push_back(read_seq_t());

    std::vector<bam1_t*> remapped;
    while (!region_scores.empty()) {

        std::vector<read_realignment_t> virus_read_realignments;
        std::pair<region_t*, bool> best_virus_region;
        google::dense_hash_set<std::string> dedup_qnames;
        while (!region_scores.empty()) {
            region_score_t best_region_score = region_scores.top();
            region_scores.pop();

            if (forced_region_score_id >= 0 && best_region_score.id != forced_region_score_id) continue;
            std::cerr << "CURRENT TOP REGION: " << best_region_score.to_str() << std::endl;

            // remap virus reads and find best virus regions
            virus_read_realignments.clear();
            best_virus_region = remap_virus_reads(best_region_score, virus_read_realignments, virus_reads_by_name, already_used);

            if (best_virus_region.first == NULL) continue;

            region_scores.push(best_region_score);
            if (best_region_score.id == region_scores.top().id) {
                region_t* best_region = best_region_score.region;
                bool host_region_fwd = best_region_score.fwd;
                region_scores.pop();

                std::vector<read_realignment_t>& host_read_realignments = host_region_fwd ?
                        reads_per_region[best_region->id] : reads_per_region_rc[best_region->id];
                dedup_qnames = dedup_reads(host_read_realignments, virus_read_realignments);

                best_region_score.reads = 0;
                best_region_score.score = 0;
                for (read_realignment_t& rr : host_read_realignments) {
                    if (dedup_qnames.count(bam_get_qname(rr.read))) {
                        best_region_score.reads++;
                        best_region_score.score += rr.score;
                    }
                }

                region_scores.push(best_region_score);
                if (best_region_score.id == region_scores.top().id) {
                    std::cerr << "CONFIRMED TOP REGION: " << best_region_score.to_str() << std::endl;
                    break;
                }
            }
        }

        if (region_scores.empty()) break;

        region_score_t best_region_score = region_scores.top();
        region_t* best_region = best_region_score.region;
        bool host_region_fwd = best_region_score.fwd;
        if (best_region_score.score == 0) {
            std::cerr << "REMAINING REGIONS HAVE SCORE 0" << std::endl;
            break;
        }


        // remove host reads s.t. the virus read was not remapped to best virus region
        std::vector<read_realignment_t>& host_read_realignments = host_region_fwd ? reads_per_region[best_region->id] : reads_per_region_rc[best_region->id];
        std::unordered_set<std::string> remapped_virus_qnames;
        for (read_realignment_t& rr : virus_read_realignments) {
            remapped_virus_qnames.insert(bam_get_qname(rr.read));
        }
        host_read_realignments.erase(
                std::remove_if(host_read_realignments.begin(), host_read_realignments.end(), [&remapped_virus_qnames](read_realignment_t& ras) {
                    return !remapped_virus_qnames.count(bam_get_qname(ras.read));
                }),
                host_read_realignments.end());

        for (read_realignment_t& rr : host_read_realignments) {
            already_used.insert(rr.read);
        }

        std::cerr << "HOST READS: " << host_read_realignments.size() << std::endl;


        samFile* writer = open_bam_writer(workspace+"/readsx", std::to_string(virus_integration_id)+".bam", host_and_virus_file->header);

        std::vector<double> host_score_ratios, virus_score_ratios;

        // edit virus reads and find host bp
        edit_remapped_reads(best_virus_region.first, virus_read_realignments, !best_virus_region.second);

        std::vector<read_realignment_t> virus_deduped_realignments;
        for (read_realignment_t& rr : virus_read_realignments) {
            if (dedup_qnames.count(bam_get_qname(rr.read))) {
                virus_deduped_realignments.push_back(rr);
                virus_score_ratios.push_back(rr.score/double(rr.offset_end-rr.offset_start+1));
            }
        }
        for (read_realignment_t& rr : virus_deduped_realignments) {
            int ok = sam_write1(writer, host_and_virus_file->header, rr.read);
            if (ok < 0) throw "Failed to write to " + std::string(writer->fn);
        }
        int v_min_pos = INT32_MAX, v_max_pos = 0;
        for (read_realignment_t& rr : virus_deduped_realignments) {
            v_min_pos = std::min(v_min_pos, rr.read->core.pos);
            v_max_pos = std::max(v_max_pos, bam_endpos(rr.read));
        }
        int_breakpoint_t virus_bp(best_virus_region.first->contig_id, v_min_pos, v_max_pos, best_virus_region.second);

        // edit host reads and find host bp
        edit_remapped_reads(best_region, host_read_realignments, !host_region_fwd);

        std::vector<read_realignment_t> host_deduped_realignments;
        for (read_realignment_t& rr : host_read_realignments) {
            if (dedup_qnames.count(bam_get_qname(rr.read))) {
                host_deduped_realignments.push_back(rr);
                host_score_ratios.push_back(rr.score/double(rr.offset_end-rr.offset_start+1));
            }
        }
        for (read_realignment_t& rr : host_deduped_realignments) {
            int ok = sam_write1(writer, host_and_virus_file->header, rr.read);
            if (ok < 0) throw "Failed to write to " + std::string(writer->fn);
        }
        int score = 0, h_min_pos = INT32_MAX, h_max_pos = 0;
        for (read_realignment_t& rr : host_deduped_realignments) {
            h_min_pos = std::min(h_min_pos, rr.read->core.pos);
            h_max_pos = std::max(h_max_pos, bam_endpos(rr.read));
            score += rr.score;
        }
        int_breakpoint_t host_bp(best_region->contig_id, h_min_pos, h_max_pos, host_region_fwd);

        sam_close(writer);

        std::vector<std::pair<double, std::string> > sample_w_qname = get_dists_for_KStest(host_deduped_realignments, host_bp,
                                                          virus_deduped_realignments, virus_bp);

        std::ofstream dist_file(workspace + "/readsx/" + std::to_string(virus_integration_id) + ".dist");
        std::vector<double> sample;
        for (auto p : sample_w_qname) {
            if (!virus_reads_by_name.count(p.second)) continue; // keep only reads where one mate maps to host and one to virus
            dist_file << p.second << " " << p.first << std::endl;
            sample.push_back(p.first);
        }
        dist_file.close();

        double p_value = -1.0;
        if (sample.size() >= 5) {
//            if (sample.size() > 100) {
//                std::random_shuffle(sample.begin(), sample.end());
//                sample.resize(100);
//            }
//            std::vector<double> pop2 = population;
//            if (population.size() > 100) {
//                std::random_shuffle(population.begin(), population.end());
//                pop2.resize(100);
//            }
//            p_value = ks_test(pop2, sample);
            alglib::real_1d_array sample_v;
            sample_v.setcontent(sample.size(), sample.data());
            for (double d : population) {
            }

            double p1, p2;
            alglib::mannwhitneyutest(pop_v, pop_v.length(), sample_v, sample_v.length(), p_value, p1, p2);
        }

        virus_integration_t v_int(best_region->contig_id, host_bp, virus_bp,
                                  best_region_score.reads, best_region_score.score, p_value, mean(sample),
                                  mean(host_score_ratios), mean(virus_score_ratios));
        std::cout << v_int.to_str() << std::endl;

        std::cerr << "INTEGRATION ID: " << v_int.id << std::endl;
        std::cerr << "ALREADY USED: " << already_used.size() << std::endl << std::endl;

        region_scores.pop();
        std::cerr << "REMAINING REGIONS: " << region_scores.size() << std::endl;
    }

    /* ==== */

    close_samFile(host_file);
    close_samFile(host_and_virus_file);
}
