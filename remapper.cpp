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
#include <sys/time.h>

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

config_t config;
stats_t stats;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;
std::unordered_map<int, int> contig_id2tid, contig_tid2id;
std::unordered_map<std::string, std::pair<char*, size_t> > chrs;

std::unordered_map<std::string, uint32_t> cigar_ids;

const char fwd_char = 'F', rev_char = 'R';

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

    virus_integration_t(int contig_id, int_breakpoint_t& host_bp, int_breakpoint_t& virus_bp, int reads, int score) :
    id(virus_integration_id++), host_bp(host_bp), virus_bp(virus_bp), reads(reads), score(score) {}

    std::string to_str() {
        char buffer[4096];
        sprintf(buffer, "ID=%d %s %s %d %d", id, host_bp.to_str().c_str(), virus_bp.to_str().c_str(), reads, score);
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

struct realign_info_t {
    uint16_t offset_start, offset_end;
    uint8_t cigar_len;
    const uint32_t* cigar;
    uint16_t score;

    realign_info_t() : offset_start(0), offset_end(0), cigar_len(0), cigar(NULL), score(0) {};
    realign_info_t(uint16_t offset_start, uint16_t offset_end, uint8_t cigar_len, const uint32_t* cigar, uint16_t score) :
        offset_start(offset_start), offset_end(offset_end), cigar_len(cigar_len), cigar(cigar), score(score) {}
};
struct read_and_score_t {
    bam1_t* read;
    realign_info_t realign_info;

    read_and_score_t() : read(NULL) {}
    read_and_score_t(bam1_t* read) : read(read) {};
    read_and_score_t(bam1_t* read, uint16_t offset_start, uint16_t offset_end, uint8_t cigar_len, const uint32_t* cigar, uint16_t score)
    : read(read), realign_info(offset_start, offset_end, cigar_len, cigar, score) {}

    bool accepted() { return realign_info.cigar != NULL; }
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

void extract_regions(std::vector<bam1_t *>& reads, std::vector<region_t*>& regions, bool process_xa = true) {

    if (reads.empty()) return;

    std::vector<cluster_t*> clusters;
    for (bam1_t* read : reads) {
        anchor_t a(bam_is_rev(read) ? 'L' : 'R', contig_tid2id[read->core.tid], read->core.pos, bam_endpos(read), 0);
        cluster_t* c = new cluster_t(a, a, DISC_TYPES.DC, 1);
        c->id = clusters.size();
        clusters.push_back(c);

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
    max_dists.push_back(stats.max_is);

    for (int max_dist : max_dists) {
        std::multimap<int, cluster_t*> clusters_map;
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
            new_cluster->id = std::max(ccd.c1->id, ccd.c2->id); // seq clusters have id -1
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
        if (!c->dead) {
            regions.push_back(new region_t(c->a1.contig_id, contig_id2tid[c->a1.contig_id],
                                           std::max(0, c->a1.start - stats.max_is + c->a1.size()),
                                           std::min(c->a1.end + stats.max_is - c->a1.size(),
                                                    (int) chrs[contig_id2name[c->a1.contig_id]].second)));
        };
        delete c;
    }
}

// TODO: unify with isolate_relevant_pairs
std::vector<region_t*>* kmer_index, * kmer_index_rc;
std::vector<read_and_score_t>* reads_per_region, * reads_per_region_rc;
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


inline void remap_read(std::string seq, bam1_t* read, int read_id, region_t* region, bool is_rc, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
                  region->end - region->start, filter, &alignment, 0);

    if (!accept_alignment(alignment, seq)) return;

    std::vector<read_and_score_t>& reads_per_region_v = is_rc ? reads_per_region_rc[region->id] : reads_per_region[region->id];

    std::string cigar = alignment_cigar_to_bam_cigar(alignment.cigar);
    auto cigar_v = cigar_str_to_array(cigar);

    mtx_regions[region->id].lock();
    reads_per_region_v.push_back(read_and_score_t(read, alignment.ref_begin, alignment.ref_end, cigar_v.first, cigar_v.second, alignment.sw_score));
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

void associate_reads_to_regions(int id, std::vector<bam1_t*>* reference_reads, int start, int end) {
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
        bam1_t* read = (*reference_reads)[r-1];

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
    std::vector<read_and_score_t>& read_scores = rc ? reads_per_region_rc[region->id] : reads_per_region[region->id];
    for (read_and_score_t& read_and_score : read_scores) {
        if (already_used.count(read_and_score.read) == 0) {
            score += read_and_score.realign_info.score;
            reads++;
        }
    }
    return std::make_pair(score, reads);
};

void del_aux(bam1_t* read, const char* tag) {
    uint8_t* b = bam_aux_get(read, tag);
    if (b == NULL) return;
    bam_aux_del(read, b);
}

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

void edit_remapped_reads(region_t* region, std::vector<read_and_score_t>& read_scores, bool rc) {
    for (read_and_score_t read_and_score : read_scores) {
        bam1_t* read = read_and_score.read;

        read->core.tid = region->original_bam_id;
        read->core.pos = region->start + read_and_score.realign_info.offset_start;
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

        del_aux(read, "AS");
        del_aux(read, "XS");
        del_aux(read, "XA");
        del_aux(read, "SA");
        del_aux(read, "NM");
        del_aux(read, "MD");

        int l_aux = bam_get_l_aux(read);
        int l_data = read->core.l_qname + 4*read_and_score.realign_info.cigar_len + (read->core.l_qseq+1)/2
                + read->core.l_qseq + l_aux;
        uint32_t m_data = l_data;
        kroundup32(m_data);
        uint8_t* data = new uint8_t[m_data];
        memset(data, 0, m_data);

        uint8_t* mov_data = data;
        memcpy(mov_data, (uint8_t*) bam_get_qname(read), read->core.l_qname);
        mov_data += read->core.l_qname;
        memcpy(mov_data, read_and_score.realign_info.cigar, 4*read_and_score.realign_info.cigar_len);
        mov_data += 4*read_and_score.realign_info.cigar_len;
        memcpy(mov_data, (uint8_t*) bam_get_seq(read), (read->core.l_qseq+1)/2);
        mov_data += (read->core.l_qseq+1)/2;
        memcpy(mov_data, (uint8_t*) bam_get_qual(read), read->core.l_qseq);
        mov_data += read->core.l_qseq;
        memcpy(mov_data, (uint8_t*) bam_get_aux(read), l_aux);

        read->l_data = l_data;
        read->m_data = m_data;
        read->data = data;
        read->core.n_cigar = read_and_score.realign_info.cigar_len;
    }
}


void write_qnames_indices(std::string& workspace, std::vector<bam1_t *>& reads) {
    std::ofstream qnames_map_out(workspace + "/qnames-map");
    for (int i = 0; i < reads.size(); i++) {
        qnames_map_out << i << " " << bam_get_qname(reads[i]) << " " << get_cigar_code(reads[i]) << std::endl;
    }
    qnames_map_out.close();
}
void load_qnames_indices(std::string& workspace, std::vector<bam1_t*> reads, std::vector<bam1_t*>& id_to_read) {
    std::ifstream qnames_map_in(workspace + "/qnames-map");

    std::unordered_map<std::string, std::vector<int> > qname_to_id;
    int i; std::string qname, cigar;
    while (qnames_map_in >> i >> qname >> cigar) {
        qname_to_id[qname + " " + cigar].push_back(i);
    }
    qnames_map_in.close();

    id_to_read.resize(i+1);

    for (bam1_t* r : reads) {
        std::string qname = bam_get_qname(r), cigar = get_cigar_code(r);
        for (int i : qname_to_id[qname + " " + cigar]) {
            id_to_read[i] = r;
        }
    }
}
void load_cigars_indices(std::string& workspace, std::vector<std::pair<int, const uint32_t*> >& id_to_cigar) {
    std::ifstream cigars_map_in(workspace + "/cigars-map");
    int i; std::string cigar_str;

    std::vector<std::pair<int, std::string> > temp;
    while (cigars_map_in >> i >> cigar_str) {
        temp.push_back(std::make_pair(i, cigar_str));
    }

    id_to_cigar.resize(temp.size()+1);
    for (std::pair<int, std::string>& p : temp) {
        id_to_cigar[p.first] = cigar_str_to_array(p.second);
    }
}
void load_qnames_indices(std::string& workspace, std::vector<bam1_t*> reads, std::unordered_map<uint32_t , bam1_t*>& id_to_read) {
    std::ifstream qnames_map_in(workspace + "/qnames-map");

    std::unordered_map<std::string, std::vector<int> > qname_to_id;
    int i; std::string qname, cigar;
    while (qnames_map_in >> i >> qname >> cigar) {
        qname_to_id[qname + " " + cigar].push_back(i);
    }
    qnames_map_in.close();

    for (bam1_t* r : reads) {
        std::string qname = bam_get_qname(r), cigar = get_cigar_code(r);
        for (int i : qname_to_id[qname + " " + cigar]) {
            id_to_read[i] = r;
        }
    }
}
void load_cigars_indices(std::string& workspace, std::unordered_map<uint32_t, std::pair<int, const uint32_t*> >& id_to_cigar) {
    std::ifstream cigars_map_in(workspace + "/cigars-map");
    int i; std::string cigar_str;
    while (cigars_map_in >> i >> cigar_str) {
        id_to_cigar[i] = cigar_str_to_array(cigar_str);
    }
}


const int APPROX_FACTOR = 500;
std::unordered_set<std::string> existing_caches;
std::unordered_map<std::string, google::dense_hash_map<uint64_t, read_and_score_t> > cached_virus_alignments;

uint64_t virus_read_id = 1;
std::vector<read_seq_t> virus_reads_seqs;

int remap_virus_reads_supp(region_t* region, std::vector<uint64_t>& read_seq_ids, bool rc,
        std::vector<read_and_score_t>* virus_read_scores,
        StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {

    std::string region_str = (rc ? "-" : "+") + region->to_str();

    // FIXME: reads appearing as both read and anchor should have two different cachings, one as reads and one as anchor
    google::dense_hash_map<uint64_t, read_and_score_t>& cached_region_alignments = cached_virus_alignments[region_str];
    if (!existing_caches.count(region_str)) {
        cached_region_alignments.set_empty_key(0);
        existing_caches.insert(region_str);
    }

    std::vector<read_and_score_t> virus_read_scores_local;
    for (uint64_t read_seq_id : read_seq_ids) {
        if (!cached_region_alignments.count(read_seq_id)) {
            StripedSmithWaterman::Alignment alignment;
            read_seq_t &read_seq = virus_reads_seqs[read_seq_id];
            std::string &read_seq_str = rc ? read_seq.seq_rc : read_seq.seq;
            aligner.Align(read_seq_str.c_str(), chrs[contig_id2name[region->contig_id]].first + region->start,
                          region->len(), filter, &alignment, 0);

            if (accept_alignment(alignment, read_seq_str)) {
                std::string cigar_str = alignment_cigar_to_bam_cigar(alignment.cigar);
                std::pair<int, const uint32_t *> cigar = cigar_str_to_array(cigar_str);
                cached_region_alignments[read_seq_id] = read_and_score_t(read_seq.read, alignment.ref_begin, alignment.ref_end,
                                                                         cigar.first, cigar.second, alignment.sw_score);
            } else {
                // cache failure to align
                cached_region_alignments[read_seq_id] = read_and_score_t(read_seq.read);
            }
        }

        read_and_score_t& ras = cached_region_alignments[read_seq_id];
        if (!ras.accepted()) continue;

        virus_read_scores_local.push_back(ras);
    }

    // find clip position
    int clip_positions[1000000];
    std::fill(clip_positions, clip_positions+region->len(), 0);
    for (read_and_score_t& ras : virus_read_scores_local) {
        if (rc && cigar_int_to_op(ras.realign_info.cigar[0]) == 'S') {
            clip_positions[ras.realign_info.offset_start]++;
        }
        if (!rc &&  cigar_int_to_op(ras.realign_info.cigar[ras.realign_info.cigar_len-1]) == 'S') {
            clip_positions[ras.realign_info.offset_end]++;
        }
    }
    int clip_pos = std::distance(clip_positions, std::max_element(clip_positions, clip_positions+region->len()));

    int accept_window_start, accept_window_end;
    if (rc) {
        accept_window_start = clip_pos;
        accept_window_end = accept_window_start + stats.max_is;
    } else {
        accept_window_end = clip_pos;
        accept_window_start = accept_window_end - stats.max_is;
    }

    int score = 0;
    for (read_and_score_t& ras : virus_read_scores_local) {
        if (ras.realign_info.offset_start >= accept_window_start && ras.realign_info.offset_end <= accept_window_end) {
            score += ras.realign_info.score;
            if (virus_read_scores != NULL) {
                virus_read_scores->push_back(ras);
            }
        }
    }
    return score;
}

std::pair<region_t*, bool> remap_virus_reads(region_score_t& host_region_score, std::vector<read_and_score_t>& virus_read_scores,
        google::dense_hash_map<std::string, bam1_t*>& virus_reads_by_name, std::unordered_set<bam1_t*>& already_used) {
    StripedSmithWaterman::Aligner aligner(1, 2, 4, 1, false);
    StripedSmithWaterman::Filter filter;

    std::vector<read_and_score_t>& host_read_scores = host_region_score.fwd ?
            reads_per_region[host_region_score.region->id] : reads_per_region_rc[host_region_score.region->id];
    // remove already used read_and_score
    host_read_scores.erase(
            std::remove_if(host_read_scores.begin(), host_read_scores.end(), [&already_used](read_and_score_t& ras) {
                return already_used.count(ras.read);
            }),
            host_read_scores.end());
    host_region_score.reads = host_read_scores.size();

    // retrieve virus mates
    std::vector<bam1_t*> remapped_reads_mates;
    std::vector<bam1_t*> remapped_host_anchors;
    for (read_and_score_t& read_score : host_read_scores) {
        std::string qname = bam_get_qname(read_score.read);
        if (virus_reads_by_name.count(qname)) { // TODO: it would be useful to use successful host clips for region finding
            remapped_reads_mates.push_back(virus_reads_by_name[qname]);
        }

        if ((host_region_score.fwd && is_right_clipped(read_score.read))
        || (!host_region_score.fwd && is_left_clipped(read_score.read))) {
            remapped_host_anchors.push_back(read_score.read);
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

    if (virus_regions.empty()) return std::make_pair(nullptr, false);

    std::vector<uint64_t> read_seq_ids;
    read_seq_ids.reserve(remapped_reads_mates.size()+remapped_host_anchors.size());
    for (bam1_t* read : remapped_reads_mates) {
        if (!read->id) {
            read->id = virus_read_id++;
            std::string seq = get_sequence(read, true);
            virus_reads_seqs.push_back(read_seq_t(read, seq));
        }
        read_seq_ids.push_back(read->id);
    }
    for (bam1_t* anchor : remapped_host_anchors) {
        if (!anchor->id) {
            anchor->id = virus_read_id++;
            std::string clip;
//            if (is_right_clipped(anchor)) {
//                clip = get_sequence(anchor).substr(anchor->core.l_qseq-get_right_clip_len(anchor), get_right_clip_len(anchor));
//                get_rc(clip);
//            } else {
//                clip = get_sequence(anchor).substr(0, get_left_clip_len(anchor));
//            }
            clip = get_sequence(anchor, true);
            virus_reads_seqs.push_back(read_seq_t(anchor, clip));
        }
        read_seq_ids.push_back(anchor->id);
    }

    region_t* best_region = NULL;
    int best_score = -1;
    bool is_best_fwd;
    for (region_t* region : virus_regions) {
        int score = remap_virus_reads_supp(region, read_seq_ids, false, NULL, aligner, filter);
        int rc_score = remap_virus_reads_supp(region, read_seq_ids, true, NULL, aligner, filter);

        int max_score = std::max(score, rc_score);
        if (best_score < max_score) {
            best_score = max_score;
            best_region = region;
            is_best_fwd = (max_score == score);
        }
    }

    // FIXME: what if best_region is null?
    std::cerr << "BEST VIRUS REGION: " << best_region->to_str() << " " << best_score << std::endl;
    remap_virus_reads_supp(best_region, read_seq_ids, !is_best_fwd, &virus_read_scores, aligner, filter);

    google::dense_hash_set<std::string> remapped_virus_qnames;
    remapped_virus_qnames.set_empty_key("");
    for (read_and_score_t& read_and_score : virus_read_scores) {
        remapped_virus_qnames.insert(bam_get_qname(read_and_score.read));
    }

    std::cerr << "BEFORE " << host_region_score.to_str() << std::endl;
    host_region_score.score = 0;
    for (read_and_score_t& ras : host_read_scores) {
        if (remapped_virus_qnames.count(bam_get_qname(ras.read))) {
            host_region_score.score += ras.realign_info.score;
        }
    }
    std::cerr << "AFTER " << host_region_score.to_str() << std::endl << std::endl;

    return std::make_pair(best_region, is_best_fwd);
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
    std::string workspace = argv[4];

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
        chrs[std::string(seq->name.s)] = std::make_pair(new char[seq->seq.l+1], seq->seq.l);
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
    }
    kseq_destroy(seq);
    fclose(fasta);

    // read virus list
    fasta = fopen(virus_reference_fname.c_str(), "r");
    seq = kseq_init(fileno(fasta));
    while ((l = kseq_read(seq)) >= 0) {
        virus_names.insert(seq->name.s);
        chrs[std::string(seq->name.s)] = std::make_pair(new char[seq->seq.l+1], seq->seq.l);
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
    }
    kseq_destroy(seq);
    fclose(fasta);

    config = parse_config(workdir + "/config.txt");
    stats = parse_stats(workspace + "/stats.txt");

    open_samFile_t* virus_file = open_samFile((workspace + "/virus-side.sorted.bam").data(), true);
    open_samFile_t* host_file = open_samFile((workspace + "/host-side.sorted.bam").data(), true);
    open_samFile_t* host_and_virus_file = open_samFile((workspace + "/retained-pairs-remapped.sorted.bam").data()); // just for the header
    for (int i = 0; i < host_file->header->n_targets; i++) {
        int contig_id = contig_name2id[host_file->header->target_name[i]];
        contig_id2tid[contig_id] = i;
        contig_tid2id[i] = contig_id;
    }


    /* === READ CHIMERIC PAIRS/READS === */

    /* == Pouring host-reads into host_reads == */

    std::unordered_set<std::string> virus_qnames;
    bam1_t* read = bam_init1();
    while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        virus_qnames.insert(qname);
    }
    close_samFile(virus_file);

    // in order to dedup reads inserted twice (perhaps once because of discordant pair and once because of clips)
    std::unordered_set<std::string> host_hashes;
    auto read_hash = [](bam1_t* read) {
        return std::string(bam_get_qname(read)) + "_" + (is_first_in_pair(read) ? "1" : "2");
    };

    // extracting paired reference-virus reads
    std::vector<bam1_t*> host_reads;
    std::unordered_set<std::string> host_qnames;
    while (sam_read1(host_file->file, host_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (virus_qnames.count(qname)) {
            std::string h = read_hash(read);
            if (!host_hashes.count(h)) {
                host_reads.push_back(bam_dup1(read));
                host_hashes.insert(h);
                host_qnames.insert(qname);
            }
        }
    }

    virus_file = open_samFile((workspace + "/virus-side.sorted.bam").data(), true);
    std::vector<bam1_t*> virus_reads;
    while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (host_qnames.count(qname)) {
            virus_reads.push_back(bam_dup1(read));
        }
    }
    close_samFile(virus_file);

    /* == */

    /* == Pour host anchors having their seq remapped to virus into host_reads == */

    std::unordered_set<std::string> successful_host_clips; // host clips that map to host
    open_samFile_t* host_clips_file = open_samFile((workspace + "/host-clips.sorted.bam").data(), true);
    while (sam_read1(host_clips_file->file, host_clips_file->header, read) >= 0) {
        if (virus_names.count(host_clips_file->header->target_name[read->core.tid]) ) { // seq was mapped to virus
            std::string qname = bam_get_qname(read);
            successful_host_clips.insert(qname);
        }
    }
    close_samFile(host_clips_file);

    open_samFile_t* host_anchors_file = open_samFile((workspace + "/host-anchors.sorted.bam").data(), true);
    while (sam_read1(host_anchors_file->file, host_anchors_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (successful_host_clips.count(qname)) {
            std::string h = read_hash(read);
            if (!host_hashes.count(h)) {
                host_reads.push_back(bam_dup1(read));
                host_hashes.insert(h);
            }
        }
    }
    close_samFile(host_anchors_file);

    /* == */

    google::dense_hash_map<std::string, bam1_t*> virus_reads_by_name;
    virus_reads_by_name.set_empty_key("");
    for (bam1_t* read : virus_reads) {
        std::string qname = bam_get_qname(read);
        virus_reads_by_name[qname] = read;
    }

    std::cerr << "HOST READS: " << host_reads.size() << std::endl;

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

        std::sort(host_regions.begin(), host_regions.end(), [](const region_t *reg1, const region_t *reg2) {
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

    reads_per_region = new std::vector<read_and_score_t>[host_regions.size()];
    reads_per_region_rc = new std::vector<read_and_score_t>[host_regions.size()];

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
            uint16_t offset = *((uint16_t *) (line + 12));
            uint16_t score = *((uint16_t *) (line + 14));
            if (!is_rc) {
                reads_per_region[region_id].emplace_back(read, offset, offset, cigar_v.first, cigar_v.second, score);
            } else {
                reads_per_region_rc[region_id].emplace_back(read, offset, offset, cigar_v.first, cigar_v.second, score);
            }
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
        auto read_score_cmp = [](const read_and_score_t& ras1, const read_and_score_t& ras2) {
            return ras1.realign_info.score > ras2.realign_info.score;
        };
        std::sort(reads_per_region[i].begin(), reads_per_region[i].end(), read_score_cmp);
        std::sort(reads_per_region_rc[i].begin(), reads_per_region_rc[i].end(), read_score_cmp);
    }


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

        int best_region_score_id = -1;
        std::vector<read_and_score_t> virus_read_scores;
        std::pair<region_t*, bool> best_virus_region;
        while (!region_scores.empty()) {
            region_score_t best_region_score = region_scores.top();
            if (best_region_score.id == best_region_score_id) {
                std::cerr << "CONFIRMED TOP REGION: " << best_region_score.to_str() << std::endl;
                break;
            }
            region_scores.pop();

            best_region_score_id = best_region_score.id;

            std::cerr << "CURRENT TOP REGION: " << best_region_score.to_str() << std::endl;

            // remap virus reads and find best virus regions
            virus_read_scores.clear();
            best_virus_region = remap_virus_reads(best_region_score, virus_read_scores, virus_reads_by_name, already_used);

            if (best_virus_region.first == NULL) continue;

            region_scores.push(best_region_score);
        }

        if (region_scores.empty()) break;

        region_t* best_region = region_scores.top().region;
        bool fwd = region_scores.top().fwd;
        if (region_scores.top().score == 0) {
            std::cerr << "REMAINING REGIONS HAVE SCORE 0" << std::endl;
            break;
        }


        // remove host reads s.t. the virus read was not remapped to best virus region
        std::vector<read_and_score_t>& host_read_scores = fwd ? reads_per_region[best_region->id] : reads_per_region_rc[best_region->id];
        std::unordered_set<std::string> remapped_virus_qnames;
        for (read_and_score_t read_and_score : virus_read_scores) {
            remapped_virus_qnames.insert(bam_get_qname(read_and_score.read));
        }
        host_read_scores.erase(
                std::remove_if(host_read_scores.begin(), host_read_scores.end(), [&remapped_virus_qnames](read_and_score_t& ras) {
                    return !remapped_virus_qnames.count(bam_get_qname(ras.read));
                }),
                host_read_scores.end());


        std::cerr << "HOST READS: " << host_read_scores.size() << std::endl;

//        std::cerr << "VIRUS READS: " << remapped_reads_mates.size() << " " << remapped_host_anchors.size() << std::endl;
//        samFile* original_writer = open_bam_writer(workspace+"/original_readsx", std::to_string(virus_integration_id)+".bam", host_and_virus_file->header);
//        edit_remapped_reads(best_region, host_read_scores, !fwd);
//        for (read_and_score_t read_and_score : host_read_scores) {
//            int ok = sam_write1(original_writer, host_and_virus_file->header, read_and_score.read);
//            if (ok < 0) throw "Failed to write to " + std::string(original_writer->fn);
//        }
//        for (read_and_score_t read_and_score : remapped_reads_mates) {
//            int ok = sam_write1(original_writer, host_and_virus_file->header, read_and_score.read);
//            if (ok < 0) throw "Failed to write to " + std::string(original_writer->fn);
//        }
//        for (read_and_score_t read_and_score : remapped_host_anchors) {
//            int ok = sam_write1(original_writer, host_and_virus_file->header, read_and_score.read);
//            if (ok < 0) throw "Failed to write to " + std::string(original_writer->fn);
//        }
//        sam_close(original_writer);

//        if (virus_regions.empty()) {
//            std::cerr << "NO VIRUS REGIONS???" << " " << remapped_reads_mates.size() << std::endl;
//            break;
//        }

        samFile* writer = open_bam_writer(workspace+"/readsx", std::to_string(virus_integration_id)+".bam", host_and_virus_file->header);

        // edit virus reads and find host bp
        edit_remapped_reads(best_virus_region.first, virus_read_scores, !best_virus_region.second);
        for (read_and_score_t read_and_score : virus_read_scores) {
            int ok = sam_write1(writer, host_and_virus_file->header, read_and_score.read);
            if (ok < 0) throw "Failed to write to " + std::string(writer->fn);
        }
        int v_min_pos = INT32_MAX, v_max_pos = 0;
        for (read_and_score_t& ras : virus_read_scores) {
            v_min_pos = std::min(v_min_pos, ras.read->core.pos);
            v_max_pos = std::max(v_max_pos, bam_endpos(ras.read));
        }
        int_breakpoint_t virus_bp(best_virus_region.first->contig_id, v_min_pos, v_max_pos, best_virus_region.second);

        // edit host reads and find host bp
        edit_remapped_reads(best_region, host_read_scores, !fwd);
        for (read_and_score_t read_and_score : host_read_scores) {
            int ok = sam_write1(writer, host_and_virus_file->header, read_and_score.read);
            if (ok < 0) throw "Failed to write to " + std::string(writer->fn);
        }
        int score = 0, h_min_pos = INT32_MAX, h_max_pos = 0;
        for (read_and_score_t& ras : host_read_scores) {
            h_min_pos = std::min(h_min_pos, ras.read->core.pos);
            h_max_pos = std::max(h_max_pos, bam_endpos(ras.read));
            score += ras.realign_info.score;
        }
        int_breakpoint_t host_bp(best_region->contig_id, h_min_pos, h_max_pos, fwd);

        sam_close(writer);

        // TODO: after virus remapping, this does not make sense as it currently is
//        std::vector<region_score_t> top10;
//        for (int i = 0; i < 10 && !region_scores.empty(); i++) {
//            top10.push_back(region_scores.top());
//            region_scores.pop();
//            std::cout << top10[i].region->to_str() << " " << (top10[i].fwd ? fwd_char : rev_char) << " "
//            << top10[i].reads << " " << top10[i].score << std::endl;
//        }
//        for (int i = 0; i < top10.size(); i++) {
//            region_scores.push(top10[i]);
//        }

        virus_integration_t v_int(best_region->contig_id, host_bp, virus_bp,
                                  host_read_scores.size(), score);
        std::cout << v_int.to_str() << std::endl;

        for (read_and_score_t read_and_score : host_read_scores) {
            already_used.insert(read_and_score.read);
        }

        std::cerr << "ALREADY USED: " << already_used.size() << std::endl << std::endl;

        region_scores.pop();
        std::cerr << "REMAINING REGIONS: " << region_scores.size() << std::endl;
    }

    /* ==== */


//    std::cerr << "VIRUS REMAPPINGS: " << virus_remappings << std::endl;
//    std::cerr << "VIRUS REGIONS REMAPPED: " << virus_regions_remapped << std::endl;


    samFile* reference_remapped_file = open_bam_writer(workspace, "reference-remapped.bam", host_file->header);
    for (bam1_t* read : remapped) {
        int ok = sam_write1(reference_remapped_file, host_file->header, read);
        if (ok < 0) throw "Failed to write to " + std::string(reference_remapped_file->fn);
    }

    sam_close(reference_remapped_file);

    close_samFile(host_file);
    close_samFile(host_and_virus_file);
}


