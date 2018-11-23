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

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "cluster.h"
#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"
#include "libs/IntervalTree.h"

typedef unsigned long long ull;

int KMER_LEN = 11;
int KMER_BITS = KMER_LEN * 2;
ull KMER_MASK = (1ll << KMER_BITS)-1;

config_t config;
stats_t stats;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;
std::unordered_map<int, int> contig_id2tid;
std::unordered_map<std::string, std::pair<char*, size_t> > chrs;

const char fwd = 'F', rev = 'R';

ull nucl_bm[256] = { 0 };
char bm_nucl[4] = { 'A', 'C', 'G', 'T' };
char nucl2chr[16];

const int SKIP_READ = -1;

int virus_integration_id = 0;

struct virus_integration_t {
    int id;
    int contig_id;
    int pos;
    char dir;
    int reads, score;

    virus_integration_t(int contig_id, int pos, char dir, int reads, int score) :
    id(virus_integration_id++), contig_id(contig_id), pos(pos), dir(dir), reads(reads), score(score) {}

    std::string print() {
        std::stringstream ss;
        ss << "ID=" << id << " " << contig_id2name[contig_id] << " " << pos << " " << pos
           << " " << reads << " " << (dir == 'R' ? reads : 0) << " "
           << (dir == 'R' ? 0 : reads) << " " << score;
        return ss.str();
    }
};


struct region_t {
    int id;
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    int score = 0;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : id(-1), contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}

    std::string to_str() {
        return contig_id2name[contig_id] + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    int len() {
        return end-start+1;
    }
};

struct realign_info_t {
    uint16_t offset;
    uint8_t cigar_len;
    const uint32_t* cigar;
    uint16_t score;

    realign_info_t(uint16_t offset, uint8_t cigar_len, const uint32_t* cigar, uint16_t score) : offset(offset), cigar_len(cigar_len), cigar(cigar), score(score) {}
};
struct read_and_score_t {
    bam1_t* read;
    realign_info_t realign_info;

    read_and_score_t(bam1_t* read, uint16_t offset, uint8_t cigar_len, const uint32_t* cigar, uint16_t score)
    : read(read), realign_info(realign_info_t(offset, cigar_len, cigar, score)) {}
};

struct region_score_t {
    region_t* region;
    bool fwd;
    int score, reads;

    region_score_t(region_t* region, bool fwd, int score, int reads) : region(region), fwd(fwd), score(score), reads(reads) {};
};


char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
};

int find(int* parents, int i) {
    int root = i;
    while (root != parents[root]) {
        root = parents[root];
    }
    while (i != root) {
        int newp = parents[i];
        parents[i] = root;
        i = newp;
    }
    return root;
}
void merge(int* parents, int* sizes, int x, int y) {
    int i = find(parents, x);
    int j = find(parents, y);
    if (i == j) return;

    if (sizes[i] < sizes[j]) {
        parents[i] = j;
        sizes[j] += sizes[i];
    } else {
        parents[j] = i;
        sizes[i] += sizes[j];
    }
}

void remove_cluster_from_mm(std::multimap<int, cluster_t*>& mm, cluster_t* c, int pos) {
    auto bounds = mm.equal_range(pos);
    for (auto it = bounds.first; it != bounds.second; it++) {
        if (it->second == c) {
            mm.erase(it);
            break;
        }
    }
}
void remove_cluster_from_mm(std::multimap<int, cluster_t*>& mm, cluster_t* c) {
    remove_cluster_from_mm(mm, c, c->a1.start);
    remove_cluster_from_mm(mm, c, c->a1.end);
}

void extract_regions(std::vector<bam1_t *>& reads, bam_hdr_t* header, std::vector<region_t*>& regions,
        bool process_xa = true) {

    std::vector<cluster_t*> clusters;
    std::multimap<int, cluster_t*> clusters_map;
    for (bam1_t* read : reads) {
        std::string tname = header->target_name[read->core.tid];
        anchor_t a(bam_is_rev(read) ? 'L' : 'R', contig_name2id[tname], read->core.pos, bam_endpos(read), 0);
        cluster_t* c = new cluster_t(a, a, DISC_TYPES.DC, 1);
        c->id = clusters.size();
        clusters.push_back(c);
        clusters_map.insert(std::make_pair(c->a1.start, c));
        clusters_map.insert(std::make_pair(c->a1.end, c));

        uint8_t* xa = bam_aux_get(read, "XA");
        if (xa != NULL && process_xa) {
            std::string xa_s = bam_aux2Z(xa);
            size_t pos = 0, prev = 0;
            while ((pos = xa_s.find(';', pos+1)) != std::string::npos) {
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

                anchor_t a(xa_dir == '-' ? 'L' : 'R', contig_name2id[xa_tname], xa_pos, xa_pos+qlen, 0);
                cluster_t* c = new cluster_t(a, a, DISC_TYPES.DC, 1);
                c->id = clusters.size();
                clusters.push_back(c);
                clusters_map.insert(std::make_pair(c->a1.start, c));
                clusters_map.insert(std::make_pair(c->a1.end, c));
            };
        }
    }

    sort(clusters.begin(), clusters.end(), [](const cluster_t* c1, const cluster_t* c2) {
        return c1->a1 < c2->a1;
    });

    // merge first similar clusters (essentially downsample reads)
    int prev;
    do {
        prev = clusters.size();
        for (int i = 0; i < clusters.size()-1; i++) {
            cluster_t* c1 = clusters[i], * c2 = clusters[i+1];
            if (c1 != NULL && std::abs(c1->a1.start-c2->a1.start) <= 10 && std::abs(c1->a1.end-c2->a1.end) <= 10) {
                cluster_t* new_cluster = cluster_t::merge(c1, c2);
                new_cluster->id = std::min(c1->id, c2->id);
                assert(std::min(c1->id, c2->id) >= 0);
                clusters[i] = new_cluster;
                clusters[i+1] = NULL;
            }
        }

        clusters.erase(std::remove(clusters.begin(), clusters.end(), (cluster_t*) NULL), clusters.end());
    } while (prev != clusters.size());

//    for (cluster_t* c : clip_clusters) {
//        c->id = -1;
//        clusters.push_back(c);
//        clusters_map.insert(std::make_pair(c->a1.start, c));
//        clusters_map.insert(std::make_pair(c->a1.end, c));
//    }

    std::vector<int> max_dists;
    for (int i = 0; i < 10; i++) max_dists.push_back(i);
    for (int i = 10; i < 100; i+=10) max_dists.push_back(i);
//    for (int i = 100; i < config.max_is; i+=100) max_dists.push_back(i);
    max_dists.push_back(stats.max_is);

    for (int max_dist : max_dists) {
        clusters_map.clear();
        for (cluster_t* c : clusters) {
            if (c->dead) continue;
            clusters_map.insert(std::make_pair(c->a1.start, c));
            clusters_map.insert(std::make_pair(c->a1.end, c));
        }

        std::priority_queue<cc_distance_t> pq;
        for (cluster_t* c1 : clusters) {
            if (c1->dead) continue;
            auto end = clusters_map.upper_bound(c1->a1.end+max_dist);
            for (auto map_it = clusters_map.lower_bound(c1->a1.start); map_it != end; map_it++) {
                cluster_t* c2 = map_it->second;
                if (c1 != c2 && cluster_t::can_merge(c1, c2, stats.max_is) &&
                    (c1->a1.start <= c2->a1.start)) {
                    pq.push(cc_distance_t(cluster_t::distance(c1, c2), c1, c2));
                }
            }
        }

        while (!pq.empty()) {
            cc_distance_t ccd = pq.top();
            pq.pop();

            if (ccd.c1->dead || ccd.c2->dead) continue;

            cluster_t* new_cluster = cluster_t::merge(ccd.c1, ccd.c2);
            new_cluster->id = std::max(ccd.c1->id, ccd.c2->id); // clip clusters have id -1
            clusters.push_back(new_cluster);

            ccd.c1->dead = true;
            remove_cluster_from_mm(clusters_map, ccd.c1);
            ccd.c2->dead = true;
            remove_cluster_from_mm(clusters_map, ccd.c2);

            auto end = clusters_map.upper_bound(new_cluster->a1.end + max_dist);
            for (auto map_it = clusters_map.lower_bound(new_cluster->a1.start - max_dist);
                 map_it != end; map_it++) {
                if (cluster_t::can_merge(new_cluster, map_it->second, stats.max_is)) {
                    pq.push(cc_distance_t(cluster_t::distance(new_cluster, map_it->second), new_cluster,
                                          map_it->second));
                }
            }
            clusters_map.insert(std::make_pair(new_cluster->a1.start, new_cluster));
            clusters_map.insert(std::make_pair(new_cluster->a1.end, new_cluster));
        }
    }


    for (cluster_t* c : clusters) {
        if (c->dead) continue;
        regions.push_back(new region_t(c->a1.contig_id, contig_id2tid[c->a1.contig_id],
                                             std::max(0, c->a1.start-stats.max_is+c->a1.size()),
                                             std::min(c->a1.end+stats.max_is-c->a1.size(),
                                             (int) chrs[contig_id2name[c->a1.contig_id]].second)));
    }

//    /*for (std::vector<bam1_t*>& cluster : clusters) {
//        std::sort(cluster.begin(), cluster.end(), [](const bam1_t* r1, const bam1_t* r2) {return r1->core.pos < r2->core.pos;});
//        bam1_t* fr = *(cluster.begin());
//        bam1_t* lr = *(cluster.rbegin());
//        int contig_id = contig_name2id[reference_file->header->target_name[fr->core.tid]];
//        regions.push_back(region_t(contig_id, fr->core.tid, fr->core.pos-config.max_is, lr->core.pos+config.max_is));
//    }
//
//
//    // remove clusters of size 0
//    read_clusters.erase(std::remove_if(read_clusters.begin(), read_clusters.end(),
//                                       [](std::vector<bam1_t*> v) {return v.empty();}), read_clusters.end());*/

//    for (auto& v : read_clusters) {
//        std::cout << header->target_name[v[0]->core.tid] << " " << v[0]->core.pos << " "
//                  << (bam_is_rev(v[0]) ? '-' : '+') << std::endl;
//        for (bam1_t* r : v) {
//            std::cout << bam_get_qname(r) << std::endl;
//        }
//        std::cout << std::endl;
//    }

//    for (cluster_t* c : clusters) delete c;
//
//    delete[] parents;
//    delete[] sizes;
}

int compute_score_supp(region_t& region, std::vector<bam1_t *>& reads, bool do_rc, std::vector<int>* scores,
                       std::vector<int> *offsets, std::vector<std::string> *cigars,
                        StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter) {

    std::vector<StripedSmithWaterman::Alignment> alignments;
    int reads_positions[1000000];
    std::fill(reads_positions, reads_positions+region.len(), 0);

    for (bam1_t* r : reads) {
        std::string s = get_sequence(r);
        StripedSmithWaterman::Alignment alignment;
        int mask_len = s.length() / 2;
        if (mask_len < 15) mask_len = 15;

        if ((do_rc && !bam_is_rev(r)) || (!do_rc && bam_is_rev(r))) {
            get_rc(s);
        }
        aligner.Align(s.c_str(), chrs[contig_id2name[region.contig_id]].first + region.start, region.end - region.start,
                      filter, &alignment, mask_len);

        alignments.push_back(alignment);
    }

    for (StripedSmithWaterman::Alignment& alignment : alignments) {
        reads_positions[alignment.ref_begin] += alignment.sw_score;
    }

    for (int i = 1; i < region.len(); i++) {
        reads_positions[i] += reads_positions[i-1];
    }

    int subregion_end = 0, max_score = 0;
    for (int i = stats.max_is+1; i < region.len(); i++) {
        if (reads_positions[i]-reads_positions[i-stats.max_is-1] > max_score) {
            max_score = reads_positions[i] - reads_positions[i-stats.max_is-1];
            subregion_end = i;
        }
    }

    int score = 0;
    for (int i = 0; i < alignments.size(); i++) {
//        accepted &= !is_poly_ACGT(s.c_str()+alignment.query_begin, alignment.query_end-alignment.query_begin+1);
        bool accepted = alignments[i].sw_score >= reads[i]->core.l_qseq
                        && alignments[i].ref_begin >= subregion_end-stats.max_is
                        && alignments[i].ref_begin <= subregion_end;

        if (accepted) {
            score += alignments[i].sw_score;
        }

        if (scores != NULL) {
            scores->push_back(alignments[i].sw_score);
        }
        if (offsets != NULL) {
            if (accepted) {
                offsets->push_back(alignments[i].ref_begin);
            } else {
                offsets->push_back(SKIP_READ);
            }
        }
        if (cigars != NULL) {
            // wrapper that returns M in case of = or X
            std::stringstream ss;
            char op = ' '; int len = 0;
            for (uint32_t c : alignments[i].cigar) {
                if (op != _cigar_int_to_op(c)) {
                    if (op != ' ') ss << len << op;
                    op = _cigar_int_to_op(c);
                    len = cigar_int_to_len(c);
                } else {
                    len += cigar_int_to_len(c);
                }
            }
            ss << len << op;
            cigars->push_back(ss.str());
        }
    }

    region.score = score;
    return score;
}

void compute_score(region_t& region, std::vector<bam1_t*>& reads, std::vector<int>* scores, std::vector<int>* offsets,
                   std::vector<std::string>* cigars, bool& is_rc,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    int score = compute_score_supp(region, reads, false, NULL, NULL, NULL, aligner, filter);
    int rc_score = compute_score_supp(region, reads, true, NULL, NULL, NULL, aligner, filter);
    region.score = std::max(score, rc_score);
    if (score >= rc_score) {
        is_rc = false;
        if (offsets != NULL) {
            compute_score_supp(region, reads, false, scores, offsets, cigars, aligner, filter);
        }
    } else {
        is_rc = true;
        if (offsets != NULL) {
            compute_score_supp(region, reads, true, scores, offsets, cigars, aligner, filter);
        }
    }
}


region_t* remap(std::vector<region_t*>& regions, std::vector<bam1_t*>& reads, std::vector<std::pair<bam1_t*, int> >& remapped_now,
           std::vector<std::pair<bam1_t*, int> >& not_remapped) {

    StripedSmithWaterman::Aligner aligner(2, 2, 3, 1, false);
    StripedSmithWaterman::Filter filter, filter_w_cigar;
    filter_w_cigar.report_cigar = true;

    bool is_rc;
    for (region_t* region : regions) {
        compute_score(*region, reads, NULL, NULL, NULL, is_rc, aligner, filter);
    }
    std::sort(regions.begin(), regions.end(), [](const region_t* r1, const region_t* r2) {return r1->score > r2->score;});

    region_t* best_region = regions[0];

    std::vector<int> offsets, scores;
    std::vector<std::string> cigars;
    compute_score(*best_region, reads, &scores, &offsets, &cigars, is_rc, aligner, filter_w_cigar);

    int forward = 0, forward_pos = -1, reverse = 0, reverse_pos = INT32_MAX;
    for (int i = 0; i < reads.size(); i++) {
        bam1_t* read = reads[i];
        if (offsets[i] == SKIP_READ) {
            not_remapped.push_back(std::make_pair(read, scores[i]));
        } else {
            read->core.tid = best_region->original_bam_id;
            read->core.pos = best_region->start + offsets[i];
            if (is_rc && !bam_is_rev(read)) {
                read->core.flag |= BAM_FREVERSE; //sets flag to true
            } else if (!is_rc && bam_is_rev(read)) {
                read->core.flag &= ~BAM_FREVERSE; //sets flag to false
            }

            if (bam_is_rev(read)) {
                reverse_pos = std::min(reverse_pos, read->core.pos);
                reverse++;
            } else {
                forward_pos = std::max(forward_pos, bam_endpos(read));
                forward++;
            }

            remapped_now.push_back(std::make_pair(read, scores[i]));
        }
    }

    return best_region;
}

// TODO: unify with isolate_relevant_pairs
std::vector<region_t*>* kmer_index, * kmer_index_rc;
std::vector<read_and_score_t>* reads_per_region, * reads_per_region_rc;
std::mutex* mtx_kmers, * mtx_regions;
//std::unordered_set<bam1_t*>* reads_by_region;

//inline ull mask(ull kmer, int seg_n) {
//    ull first_n_segs_mask = (1ll << SEG_BITS*seg_n)-1;
//    ull first_n_segs = kmer & first_n_segs_mask;
//    ull remaining_segs_shifted = (kmer >> SEG_BITS) & ~first_n_segs_mask;
//    return (first_n_segs | remaining_segs_shifted) & MASKED_KMER_MASK;
//}

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

void index_regions(int id, std::vector<region_t*>* reference_regions, int start, int end) {
    mtx.lock();
    std::cerr << "Thread " << id << ": " << "indexing regions " << start << " to " << end << " " << std::endl;
    mtx.unlock();

    char region_str[100000];
    for (int r = start; r < end; r++) {
        region_t* region = (*reference_regions)[r];
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
std::ofstream scores_file_out;

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


bool accept(StripedSmithWaterman::Alignment& alignment, std::string& query) {
    return !is_poly_ACGT(query.c_str()+alignment.query_begin, alignment.query_end-alignment.query_begin+1)
    && alignment.sw_score >= 30; //read->core.l_qseq;
};

inline void remap_read(std::string seq, bam1_t* read, region_t* region, bool rc, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter, int mask_len) {
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
                  region->end - region->start, filter, &alignment, mask_len);

    if (!accept(alignment, seq)) return;

    std::vector<read_and_score_t>& reads_per_region_v = rc ? reads_per_region_rc[region->id] : reads_per_region[region->id];

    std::string cigar = alignment_cigar_to_bam_cigar(alignment.cigar);
    auto cigar_v = cigar_str_to_array(cigar);

    mtx_regions[region->id].lock();
    reads_per_region_v.push_back(read_and_score_t(read, alignment.ref_begin, cigar_v.first, cigar_v.second, alignment.sw_score));
    mtx_regions[region->id].unlock();
    mtx.lock();
    scores_file_out << region->id << " " << bam_get_qname(read) << " " << cigar << " " << alignment.ref_begin << " "
    << (rc ? rev : fwd) << " " << alignment.sw_score << std::endl;
    mtx.unlock();

    remappings++;
    if (remappings % 1000000 == 0) {
        std::cerr << remappings << " remappings done." << std::endl;
    }
}

void associate_reads_to_regions(int id, std::vector<bam1_t*>* reference_reads, std::unordered_map<int, IntervalTree<region_t*>* >* region_trees, int start, int end) {
    mtx.lock();
    std::cerr << "Thread " << id << ": " << "considering reads " << start << " to " << end << " " << std::endl;
    mtx.unlock();

    const int MAX_REGIONS = 1000000;
    int* last_read_for_region = new int[MAX_REGIONS];
    int* last_pos = new int[MAX_REGIONS];
    int* last_read_for_region_rc = new int[MAX_REGIONS];
    int* last_pos_rc = new int[MAX_REGIONS];
    std::fill(last_read_for_region, last_read_for_region+MAX_REGIONS, end);
    std::fill(last_read_for_region_rc, last_read_for_region_rc+MAX_REGIONS, end);

    StripedSmithWaterman::Aligner aligner(2, 2, 3, 1, false);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    for (int r = start; r < end; r++) {
        bam1_t* read = (*reference_reads)[r];
        std::string seq = get_sequence(read, true);
        std::string seq_rc = seq;
        get_rc(seq_rc);

        int mask_len = seq.length() / 2;
        if (mask_len < 15) mask_len = 15;

        if (read->core.qual >= 50) {
            std::vector<Interval<region_t*> > regions = (*region_trees)[read->core.tid]->findOverlapping(read->core.pos, bam_endpos(read));
            for (Interval<region_t*> region_i : regions) {
                region_t* region = region_i.value;
                remap_read(seq, read, region, bam_is_rev(read), aligner, filter, mask_len);
            }
            continue;
        }

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

                        remap_read(seq, read, region, false, aligner, filter, mask_len);

//                        remappings++;
//
//                        if (remappings % 1000000 == 0) {
//                            std::cerr << remappings << " remappings done." << std::endl;
//                        }
//
////                        aligner.Align(seq.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
////                                      region->end - region->start, filter, &alignment, mask_len);
//
//                        if (!accept(alignment)) continue;
//
//                        mtx_regions[region->id].lock();
//                        reads_per_region[region->id].push_back(read_and_score_t(read, alignment.sw_score));
//                        mtx_regions[region->id].unlock();
////                        mtx.lock();
////                        scores_file_out << region->id << " " << alignment.cigar_string << " " << bam_get_qname(read)
////                        << " " << fwd << " " << alignment.sw_score << std::endl;
//                        mtx.unlock();
                    } else if (abs(last_read_for_region[region->id]) != r) {
                        last_read_for_region[region->id] = r;
                        last_pos[region->id] = i;
                    }
                }

                for (region_t* region : kmer_index_rc[kmer]) {
                    if (last_read_for_region_rc[region->id] == r && i - last_pos_rc[region->id] >= KMER_LEN) {
                        last_read_for_region_rc[region->id] = -r;

                        remap_read(seq_rc, read, region, true, aligner, filter, mask_len);

//                        remappings++;
//
//                        if (remappings % 1000000 == 0) {
//                            std::cerr << remappings << " remappings done." << std::endl;
//                        }
//
////                        aligner.Align(seq_rc.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
////                                      region->end - region->start, filter, &alignment, mask_len);
//
//                        if (!accept(alignment)) continue;
//
//                        mtx_regions[region->id].lock();
//                        reads_per_region_rc[region->id].push_back(read_and_score_t(read, alignment.sw_score));
//                        mtx_regions[region->id].unlock();
//                        mtx.lock();
//                        scores_file_out << region->id << " " << region->to_str() << " " << bam_get_qname(read) << " "
//                        << rev << " " << alignment.sw_score << std::endl;
//                        mtx.unlock();
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


std::string print(ull kmer, int len) {
    char s[KMER_LEN];
    s[len] = '\0';
    while (len > 0) {
        len--;
        s[len] = bm_nucl[kmer%4];
        kmer /= 4;
    }
    return s;
}


std::pair<int, int> compute_region_score(const int region_id, bool rc, std::unordered_set<bam1_t*>& already_used) {
    int score = 0, reads = 0;
    std::vector<read_and_score_t>& read_scores = rc ? reads_per_region_rc[region_id] : reads_per_region[region_id];
    for (read_and_score_t read_and_score : read_scores) {
        if (already_used.count(read_and_score.read) == 0) {
            score += read_and_score.realign_info.score;
            reads++;
        }
    }
    return std::make_pair(score, reads);
};

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

void edit_remapped_reads(region_t* region, bool rc, std::unordered_set<bam1_t*>& already_used) {
    std::vector<read_and_score_t>& read_scores = rc ? reads_per_region_rc[region->id] : reads_per_region[region->id];
    for (read_and_score_t read_and_score : read_scores) {
        if (already_used.count(read_and_score.read) == 0) {
            bam1_t* read = read_and_score.read;
            read->core.tid = region->original_bam_id;
            read->core.pos = region->start + read_and_score.realign_info.offset;
            if ((rc && !bam_is_rev(read)) || (!rc && bam_is_rev(read))) {
                std::string seq = get_sequence(read);
                get_rc(seq);
                uint8_t* s = bam_get_seq(read);
                for (int i = 0; i < read->core.l_qseq; ++i){
                    bam1_seq_seti(s, i, seq_nt16_table[seq[i]]);
                }
            }
            if (rc) read->core.flag |= BAM_FREVERSE;
            else read->core.flag &= ~BAM_FREVERSE;
            int l_data = read->l_data - read->core.n_cigar + read_and_score.realign_info.cigar_len*4;
            uint32_t m_data = l_data;
            kroundup32(m_data);
            uint8_t* data = new uint8_t[m_data];
            memset(data, 0, m_data);
            memcpy(data, read->data, read->core.l_qname);
            memcpy(data+read->core.l_qname, read_and_score.realign_info.cigar, read_and_score.realign_info.cigar_len*4);
            memcpy(data+read->core.l_qname+read_and_score.realign_info.cigar_len*4, bam_get_seq(read),
                    read->l_data-read->core.n_cigar*4-read->core.l_qname);
            read->l_data = l_data;
            read->m_data = m_data;
            read->data = data;
            read->core.n_cigar = read_and_score.realign_info.cigar_len;
        }
    }
}


int main(int argc, char* argv[]) {

    nucl_bm['A'] = nucl_bm['a'] = 0;
    nucl_bm['C'] = nucl_bm['c'] = 1;
    nucl_bm['G'] = nucl_bm['g'] = 2;
    nucl_bm['T'] = nucl_bm['t'] = 3;
    nucl_bm['N'] = nucl_bm['n'] = 0;

    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';

    std::string reference_fname  = argv[1];
    std::string virus_reference_fname  = argv[2];
    std::string workdir = argv[3];
    std::string workspace = argv[4];

    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }

    FILE* fastaf = fopen(reference_fname.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(fastaf));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        chrs[std::string(seq->name.s)] = std::make_pair(new char[seq->seq.l+1], seq->seq.l);
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
    }
    kseq_destroy(seq);

    // read virus list
    std::unordered_set<std::string> virus_names;
    FILE* virus_ref_fasta = fopen(virus_reference_fname.c_str(), "r");
    seq = kseq_init(fileno(virus_ref_fasta));
    while ((l = kseq_read(seq)) >= 0) {
        virus_names.insert(seq->name.s);
    }
    kseq_destroy(seq);

    config = parse_config(workdir + "/config.txt");
    stats = parse_stats(workspace + "/stats.txt");

    open_samFile_t* virus_file = open_samFile((workspace + "/virus-side.sorted.bam").data(), true);
    open_samFile_t* reference_file = open_samFile((workspace + "/reference-side.sorted.bam").data(), true);
    for (int i = 0; i < reference_file->header->n_targets; i++) {
        int contig_id = contig_name2id[reference_file->header->target_name[i]];
        contig_id2tid[contig_id] = i;
    }

    std::unordered_set<std::string> virus_qnames;
    bam1_t* read = bam_init1();
    while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        virus_qnames.insert(qname);
    }
    close_samFile(virus_file);


    // in order to dedup reads inserted twice (perhaps once because of discordant pair and once because of clips)
    std::unordered_set<std::string> ref_hashes, virus_hashes;
    auto read_hash = [](bam1_t* read) {
        return std::string(bam_get_qname(read)) + "-" + get_sequence(read);
    };

    // extracting paired reference-virus reads
    std::vector<bam1_t*> reference_reads;
    std::unordered_set<std::string> reference_qnames;
    while (sam_read1(reference_file->file, reference_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (virus_qnames.count(qname) > 0) {
            std::string h = read_hash(read);
            if (ref_hashes.count(h) == 0) {
                reference_reads.push_back(bam_dup1(read));
                ref_hashes.insert(h);
                reference_qnames.insert(qname);
            }
        }
    }

    virus_file = open_samFile((workspace + "/virus-side.sorted.bam").data(), true);
    std::vector<bam1_t*> virus_reads;
    while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (reference_qnames.count(qname) > 0) {
            virus_reads.push_back(bam_dup1(read));
        }
    }
//    close_samFile(virus_file);


    // extracting good clips, clips that map to virus (resp. reference) if anchor was on reference (resp. virus)
    std::unordered_set<std::string> good_clips;

    open_samFile_t* virus_clips_file = open_samFile((workspace + "/virus-clips.sorted.bam").data(), true);
    while (sam_read1(virus_clips_file->file, virus_clips_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (virus_names.count(virus_clips_file->header->target_name[read->core.tid]) == 0) { // clip was mapped to reference
            good_clips.insert(qname);
            std::string h = read_hash(read);
            if (ref_hashes.count(h) == 0) {
                reference_reads.push_back(bam_dup1(read));
                ref_hashes.insert(h);
            }
        }
    }
    open_samFile_t* reference_clips_file = open_samFile((workspace + "/reference-clips.sorted.bam").data(), true);
    while (sam_read1(reference_clips_file->file, reference_clips_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (virus_names.count(reference_clips_file->header->target_name[read->core.tid]) > 0) { // clip was mapped to virus
            good_clips.insert(qname);
            virus_reads.push_back(bam_dup1(read));
        }
    }
    close_samFile(virus_clips_file);
    close_samFile(reference_clips_file);


    // extracting the good anchors, anchors of good clips
    open_samFile_t* virus_anchors_file = open_samFile((workspace + "/virus-anchors.sorted.bam").data(), true);
    while (sam_read1(virus_anchors_file->file, virus_anchors_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (good_clips.count(qname)) {
            virus_reads.push_back(bam_dup1(read));
        }
    }
    open_samFile_t* reference_anchors_file = open_samFile((workspace + "/reference-anchors.sorted.bam").data(), true);
    while (sam_read1(reference_anchors_file->file, reference_anchors_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (good_clips.count(qname)) {
            std::string h = read_hash(read);
            if (ref_hashes.count(h) == 0) {
                reference_reads.push_back(bam_dup1(read));
                ref_hashes.insert(h);
            }
        }
    }
    close_samFile(virus_anchors_file);
    close_samFile(reference_anchors_file);

    std::unordered_map<std::string, bam1_t*> virus_reads_by_name;
    for (bam1_t* read : virus_reads) {
        std::string qname = bam_get_qname(read);
        virus_reads_by_name[qname] = read;
    }


    std::vector<region_t*> reference_regions; //, virus_regions;
    extract_regions(reference_reads, reference_file->header, reference_regions);

    // filter virus regions
    auto is_virus_region = [&virus_names](const region_t* region) {return virus_names.count(contig_id2name[region->contig_id]); };
    reference_regions.erase(
            std::remove_if(reference_regions.begin(), reference_regions.end(), is_virus_region),
            reference_regions.end());

    std::sort(reference_regions.begin(), reference_regions.end(), [](const region_t* reg1, const region_t* reg2) {
        if (reg1->contig_id == reg2->contig_id) return reg1->start < reg2->start;
        return reg1->contig_id < reg2->contig_id;
    });
    for (int i = 0; i < reference_regions.size(); i++) {
        reference_regions[i]->id = i;
    }
    std::cerr << "REFERENCE REGIONS: " << reference_regions.size() << std::endl;

    std::string reference_regions_path = workspace + "/reference-regions";
    std::ofstream reference_regions_file(reference_regions_path);
    for (region_t* region : reference_regions) {
        reference_regions_file << region->id << " " << region->to_str() << std::endl;
    }
    reference_regions_file.close();

//    extract_regions(virus_reads, virus_file->header, virus_regions);

//    for (std::vector<bam1_t*>& cluster : clusters) {
//        std::sort(cluster.begin(), cluster.end(), [](const bam1_t* r1, const bam1_t* r2) {return r1->core.pos < r2->core.pos;});
//        bam1_t* fr = *(cluster.begin());
//        bam1_t* lr = *(cluster.rbegin());
//        int contig_id = contig_name2id[reference_file->header->target_name[fr->core.tid]];
//        regions.push_back(region_t(contig_id, fr->core.tid, fr->core.pos-config.max_is, lr->core.pos+config.max_is));
//    }

    KMER_LEN = 13;
    KMER_BITS = KMER_LEN * 2;
    KMER_MASK = (1ll << KMER_BITS)-1;

    reads_per_region = new std::vector<read_and_score_t>[reference_regions.size()];
    reads_per_region_rc = new std::vector<read_and_score_t>[reference_regions.size()];

    std::string scores_file_path = workspace + "/reads.scores";
    std::ifstream scores_file_in(scores_file_path);
    if (scores_file_in.good()) {
        std::unordered_map<std::string, bam1_t*> reads_map;
        for (bam1_t* read : reference_reads) {
            reads_map[std::string(bam_get_qname(read))] = read;
        }

        int region_id;
        uint16_t score;
        std::string bam_name, cigar;
        uint16_t offset;
        char strand;
        while (scores_file_in >> region_id >> bam_name >> cigar >> offset >> strand >> score) {
            bam1_t* read = reads_map[bam_name];
            auto cigar_v = cigar_str_to_array(cigar);
            if (strand == fwd) {
                reads_per_region[region_id].push_back(read_and_score_t(read, offset, cigar_v.first, cigar_v.second, score));
            } else {
                reads_per_region_rc[region_id].push_back(read_and_score_t(read, offset, cigar_v.first, cigar_v.second, score));
            }
        }
    } else {
        scores_file_out.open(scores_file_path);

        // Index regions
        kmer_index = new std::vector<region_t *>[1 << KMER_BITS];
        kmer_index_rc = new std::vector<region_t *>[1 << KMER_BITS];
        mtx_kmers = new std::mutex[1 << KMER_BITS];

        ctpl::thread_pool thread_pool1(config.threads);
        std::vector<std::future<void> > futures;

        int regions_per_thread = reference_regions.size() / config.threads;

        for (int i = 0; i < config.threads - 1; i++) {
            std::future<void> future = thread_pool1.push(index_regions, &reference_regions, i * regions_per_thread,
                                                         (i + 1) * regions_per_thread);
            futures.push_back(std::move(future));
        }
        std::future<void> future = thread_pool1.push(index_regions, &reference_regions,
                                                     (config.threads - 1) * regions_per_thread,
                                                     reference_regions.size());
        futures.push_back(std::move(future));

        thread_pool1.stop(true);
        for (int i = 0; i < futures.size(); i++) {
            try {
                futures[i].get();
            } catch (char const *s) {
                std::cout << s << std::endl;
            }
        }

        // Associate reads to regions
        std::cerr << "READS: " << reference_reads.size() << std::endl;

        // index reads by sequence
        std::unordered_map<std::string, std::vector<bam1_t*> > reads_by_seq;
        std::vector<std::string> read_seqs;
        for (bam1_t* r : reference_reads) {
            std::string seq = get_sequence(r);
            if (reads_by_seq[seq].empty()) {
                read_seqs.push_back(seq);
            }
            reads_by_seq[seq].push_back(r);
        }

        // randomize reads to avoid imbalances
        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(reference_reads), std::end(reference_reads), rng);

        // creating an interval tree of regions
        std::unordered_map<int, IntervalTree<region_t*>* > region_trees;
        std::unordered_map<int, std::vector<Interval<region_t*> > > intervals;
        for (region_t* region : reference_regions) {
            intervals[region->original_bam_id].push_back(Interval<region_t*>(region->start, region->end, region));
        }
        for (int i = 0; i < reference_file->header->n_targets; i++) {
            region_trees[i] = new IntervalTree<region_t*>(intervals[i]);
        }

        mtx_regions = new std::mutex[reference_regions.size()];

        ctpl::thread_pool thread_pool2(config.threads);
        futures.clear();

        int reads_chunks = config.threads * 5;
        int reads_per_thread = reference_reads.size() / reads_chunks;

        for (int i = 0; i < reads_chunks - 1; i++) {
            std::future<void> future = thread_pool2.push(associate_reads_to_regions, &reference_reads, &region_trees,
                    i * reads_per_thread, (i + 1) * reads_per_thread);
            futures.push_back(std::move(future));
        }
        future = thread_pool2.push(associate_reads_to_regions, &reference_reads, &region_trees,
                (reads_chunks - 1) * reads_per_thread, reference_reads.size());
        futures.push_back(std::move(future));

        thread_pool2.stop(true);
        for (int i = 0; i < futures.size(); i++) {
            try {
                futures[i].get();
            } catch (char const *s) {
                std::cout << s << std::endl;
            }
        }

        scores_file_out.close();
    }


    std::ofstream result_out(workspace + "/result.txt");

    std::unordered_set<bam1_t*> already_used;

    auto region_score_cmp = [](const region_score_t& r1, const region_score_t& r2) {return r1.score < r2.score;};
    std::priority_queue<region_score_t, std::vector<region_score_t>, decltype(region_score_cmp)> region_scores(region_score_cmp);

    for (int i = 0; i < reference_regions.size(); i++) {
        region_t* region = reference_regions[i];
        std::pair<int, int> score_and_reads = compute_region_score(region->id, false, already_used);

        int score = score_and_reads.first;
        int reads = score_and_reads.second;
        if (score > 0) region_scores.push(region_score_t(region, true, score, reads));

        score_and_reads = compute_region_score(region->id, true, already_used);
        score = score_and_reads.first;
        reads = score_and_reads.second;
        if (score > 0) region_scores.push(region_score_t(region, false, score, reads));
    }

//    while (!region_scores.empty()) {
//        region_score_t rst = region_scores.top();
//        std::cout << rst.region->to_str() << " " << (rst.fwd ? "FWD " : "REV ")
//                  << rst.score << " " << rst.reads << std::endl;
//        region_scores.pop();
//    }

    int id = 1;
    int prev = 0;
//    std::vector<bam1_t*> remaining = reference_reads;
    std::vector<bam1_t*> remapped;
    bool has_been_remapped = true;
    do {
        std::vector<std::pair<bam1_t*, int> > not_remapped, remapped_now;

//        bool is_rc;
//        for (region_t& region : reference_regions) {
//            compute_score(region, remaining, NULL, NULL, is_rc, aligner, filter);
//        }
//        std::sort(reference_regions.begin(), reference_regions.end(), [](const region_t& r1, const region_t& r2) {return r1.score > r2.score;});
//
//        region_t* best_reference_region = &reference_regions[0];
//
//        has_been_remapped = false;
//        std::vector<int> offsets;
//        std::vector<std::string> cigars;
//        compute_score(*best_reference_region, remaining, &offsets, &cigars, is_rc, aligner, filter_w_cigar);
//
//        int forward = 0, forward_pos = -1, reverse = 0, reverse_pos = INT32_MAX;
//        for (int i = 0; i < remaining.size(); i++) {
//            bam1_t* read = remaining[i];
//            if (offsets[i] == SKIP_READ) {
//                not_remapped.push_back(read);
//            } else {
//                read->core.tid = best_reference_region->original_bam_id;
//                read->core.pos = best_reference_region->start + offsets[i];
//                if (is_rc && !bam_is_rev(read)) {
//                    read->core.flag |= BAM_FREVERSE; //sets flag to true
//                } else if (!is_rc && bam_is_rev(read)) {
//                    read->core.flag &= ~BAM_FREVERSE; //sets flag to false
//                }
//
//                if (bam_is_rev(read)) {
//                    reverse_pos = std::min(reverse_pos, read->core.pos);
//                    reverse++;
//                } else {
//                    forward_pos = std::max(forward_pos, bam_endpos(read));
//                    forward++;
//                }
//
//                remapped_now.push_back(read);
//                has_been_remapped = true;
//            }
//        }

        if (region_scores.empty()) break;

        while (true) {
            region_score_t top_region = region_scores.top();

            region_t* best_region = top_region.region;
            bool fwd = top_region.fwd;
            auto score_and_reads = compute_region_score(best_region->id, !fwd, already_used);
            int new_score = score_and_reads.first;
            int new_reads = score_and_reads.second;
            int score = top_region.score;

            if (score == new_score) break;

            region_scores.pop();
            top_region.score = new_score;
            top_region.reads = new_reads;
            region_scores.push(top_region);
        }

        region_t* best_region = region_scores.top().region;
        bool fwd = region_scores.top().fwd;
        int score = region_scores.top().score;
        int reads = region_scores.top().reads;
        if (score == 0) break;

        edit_remapped_reads(best_region, !fwd, already_used);
        samFile* writer = open_bam_writer(workspace+"/readsx", std::to_string(virus_integration_id)+".bam", reference_file->header);
        std::vector<read_and_score_t>& read_scores = fwd ? reads_per_region[best_region->id] : reads_per_region_rc[best_region->id];
        if (read_scores.size() > 1) {
//            std::sort(read_scores.begin(), read_scores.end(), [](const read_and_score_t& r1, const read_and_score_t& r2) {
//                if (r1.read->core.tid != r2.read->core.tid) return r1.read->core.tid < r2.read->core.tid;
//                return r1.read->core.pos < r2.read->core.pos;
//            });
            for (read_and_score_t read_and_score : read_scores) {
                int ok = sam_write1(writer, reference_file->header, read_and_score.read);
                if (ok < 0) throw "Failed to write to " + std::string(writer->fn);

                std::string qname = bam_get_qname(read_and_score.read);
                ok = sam_write1(writer, reference_file->header, virus_reads_by_name[qname]);
                if (ok < 0) throw "Failed to write to " + std::string(writer->fn);

            }
            sam_close(writer);
        }

        std::vector<region_score_t> top10;
        for (int i = 0; i < 10 && !region_scores.empty(); i++) {
            top10.push_back(region_scores.top());
            region_scores.pop();
            std::cout << top10[i].region->to_str() << " " << (top10[i].fwd ? fwd : rev) << " "
            << top10[i].score << std::endl;
        }
        for (int i = 0; i < top10.size(); i++) {
            region_scores.push(top10[i]);
        }

        {
            virus_integration_t v_int(best_region->contig_id, fwd ? best_region->end : best_region->start,
                                      fwd ? 'R' : 'L', reads, score);
            std::cout << v_int.print() << std::endl;
        }

        std::vector<read_and_score_t>& used_reads = fwd ? reads_per_region[best_region->id] : reads_per_region_rc[best_region->id];
        for (read_and_score_t read_and_score : used_reads) {
            already_used.insert(read_and_score.read);
        }

        std::cout << "ALREADY USED: " << already_used.size() << std::endl << std::endl;

        continue;

//        region_t* best_reference_region = remap(reference_regions, remaining, remapped_now, not_remapped);
//
//        int forward = 0, forward_pos = -1, reverse = 0, reverse_pos = INT32_MAX;
//        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
//            bam1_t* read = read_and_score.first;
//            if (bam_is_rev(read)) {
//                reverse_pos = std::min(reverse_pos, read->core.pos);
//                reverse++;
//            } else {
//                forward_pos = std::max(forward_pos, bam_endpos(read));
//                forward++;
//            }
//        }
//        if (forward_pos == -1) {
//            forward_pos = reverse_pos;
//        } else if (reverse_pos == INT32_MAX) {
//            reverse_pos = forward_pos;
//        }
//
//
//        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
//            remapped.push_back(read_and_score.first);
//        }
//
//        // extract virus reads
//        std::unordered_set<std::string> remapped_qnames;
//        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
//            remapped_qnames.insert(bam_get_qname(read_and_score.first));
//        }
//        std::vector<bam1_t*> remapped_now_virus_reads;
//        for (bam1_t* r : virus_reads) {
//            if (remapped_qnames.count(bam_get_qname(r))) {
//                remapped_now_virus_reads.push_back(r);
//            }
//        }
////        for (region_t& region : virus_regions) {
////            compute_score(region, remapped_now_virus_reads, NULL, NULL, is_rc, aligner, filter);
////        }
//
//        std::vector<std::pair<bam1_t*, int> > v1, v2;
////        region_t* best_virus_region = remap(virus_regions, remapped_now_virus_reads, v1, v2);
//
//        int tot_score = 0;
//        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
//            tot_score += read_and_score.second;
//        }
//
//        remaining.clear();
//        for (std::pair<bam1_t*, int> read_and_score : not_remapped) {
//            remaining.push_back(read_and_score.first);
//        }
//
//
//        auto distinct_qnames = [](std::vector<std::pair<bam1_t*, int> >& reads_and_scores) {
//            std::set<std::string> s;
//            for (std::pair<bam1_t*, int> read_and_score : reads_and_scores) {
//                std::string qname = bam_get_qname(read_and_score.first);
//                s.insert(qname);
//            }
//            return s.size();
//        };
//
//        char dir = forward > 0 ? 'R' : 'L';
//        virus_integration_t v_int(best_reference_region->contig_id, dir == 'R' ? forward_pos : reverse_pos, dir, distinct_qnames(remapped_now), 0);
//
//        std::stringstream ss;
//        ss << contig_id2name[v_int.contig_id] << " " << v_int.pos << " " << v_int.pos
//           << " " << v_int.reads << " " << forward << " " << reverse << " " << tot_score << std::endl;
//        std::string int_s = ss.str();
//
//        result_out << int_s << std::endl;
//
//        std::cout << "ID=" << id++ << " " << int_s << std::endl;
////        std::cout << "VIRUS " << best_virus_region->to_str() << std::endl;
//        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
//            bam1_t* r = read_and_score.first;
//            std::cout << bam_get_qname(r) << "\t" << get_sequence(r) << "\t" << r->core.pos << "\t" << read_and_score.second << std::endl;
//        }
//        for (region_t* region : reference_regions) {
//            std::cout << region->to_str() << " " << region->score << std::endl;
//        }
//        for (bam1_t* r : remapped_now_virus_reads) {
//            std::cout << bam_get_qname(r) << " " << get_sequence(r, true) << std::endl;
//        }
//
//        std::cout << std::endl;
//        prev = remapped.size();
//
//        auto it = std::find_if(reference_regions.begin(), reference_regions.end(), [best_reference_region](region_t* elem) {
//            return elem == best_reference_region;});
//        reference_regions.erase(it);
//
//        has_been_remapped = !remapped_now.empty();
    } while (has_been_remapped); // && !remaining.empty());


    samFile* reference_remapped_file = open_bam_writer(workspace, "reference-remapped.bam", reference_file->header);
    for (bam1_t* read : remapped) {
        int ok = sam_write1(reference_remapped_file, reference_file->header, read);
        if (ok < 0) throw "Failed to write to " + std::string(reference_remapped_file->fn);
    }

    sam_close(reference_remapped_file);

    close_samFile(reference_file);
    result_out.close();
}

