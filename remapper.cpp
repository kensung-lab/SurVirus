#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <htslib/sam.h>
#include <htslib/hts_endian.h>
#include <cstring>
#include <algorithm>
#include <random>

#include "sam_utils.h"
#include "utils.h"
#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

typedef unsigned long long ull;

config_t config;
int max_is;

contig_map_t contig_map;
chr_seqs_map_t chr_seqs;
std::unordered_map<int, int> contig_id2tid, contig_tid2id;

std::vector<std::pair<int, const uint32_t*> > cid_to_cigar;


int virus_integration_id = 0;

std::string print_region(region_t* region) {
    return contig_map.id2name(region->contig_id) + ":" + std::to_string(region->start) + "-" + std::to_string(region->end);
}


struct region_score_t {
    int id;
    region_t* region;
    bool rev;
    int score, reads, good_reads;

    region_score_t(int id, region_t* region, bool rev, int score, int reads, int good_reads)
    : id(id), region(region), rev(rev), score(score), reads(reads), good_reads(good_reads) {};

    std::string to_string() {
        return std::to_string(id) + " " + print_region(region) + " READS=" + std::to_string(reads) + " SCORE=" + std::to_string(score)
               + (rev ? " FWD" : " REV");
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


bool is_host_region(const region_t *region) {
    return !chr_seqs.is_virus(contig_map.id2name(region->contig_id));
}



// TODO: unify with isolate_relevant_pairs
std::vector<read_realignment_t>* reads_per_region, * reads_per_region_rc;

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
}


void edit_remapped_reads(region_t* region, std::vector<read_realignment_t>& read_realignments) {
    int region_tid = contig_id2tid[region->contig_id];
    for (read_realignment_t rr : read_realignments) {
        bam1_t* read = rr.read;
        read->core.tid = region_tid;
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
            for (int i = region_len-1-max_is; i >= 0; i--) {
                score_windows[i] -= score_points[i+max_is];
            }
        } else {
            for (int i = max_is; i < region_len; i++) {
                score_windows[i] -= score_points[i-max_is];
            }
        }

        bp_pos = rc ? get_last_max_element(score_windows, region_len) : get_first_max_element(score_windows, region_len);
    }

    int accept_window_start, accept_window_end;
    if (rc) {
        accept_window_start = bp_pos;
        accept_window_end = accept_window_start + max_is;
    } else {
        accept_window_end = bp_pos;
        accept_window_start = accept_window_end - max_is;
    }
    return window_t(accept_window_start, accept_window_end);
}


google::dense_hash_map<std::string, int> get_distances(std::vector<read_realignment_t>& realignments, breakpoint_t& bp, region_t region) {
    google::dense_hash_map<std::string, int> dists;
    dists.set_empty_key("");
    for (read_realignment_t& rr : realignments) {
        if (rr.read->core.flag & BAM_FSECONDARY) continue;
        std::string qname = bam_get_qname(rr.read);

        int read_endpos;
        if (bp.rev) { // include spurious soft-clipping as part of the distance
            read_endpos = region.start + rr.offset_end + rr.right_clip_len();
        } else {
            read_endpos = region.start + rr.offset_start - rr.left_clip_len();
        }

        int dist = abs(read_endpos - bp.pos());
        if (!dists.count(qname) || dists[qname] < dist) {
            dists[qname] = dist;
        }
    }
    return dists;
}
google::dense_hash_map<std::string, int> get_distances(std::vector<read_realignment_t>& realignments, breakpoint_t& bp) {
    return get_distances(realignments, bp, region_t(0,0,0));
}
google::dense_hash_map<std::string, std::vector<std::string> > dedup_reads(
        std::vector<read_realignment_t>& host_realignments, std::vector<read_realignment_t> virus_realignments,
        region_t* host_region, bool host_rev, region_t* virus_region, bool virus_fwd) {
    if (host_realignments.empty() || virus_realignments.empty()) {
        google::dense_hash_map<std::string, std::vector<std::string> > empty_set;
        empty_set.set_empty_key("");
        return empty_set;
    }

    bool virus_rev = !virus_fwd;

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
    breakpoint_t host_bp = breakpoint_t(contig_map.id2name(host_region->contig_id), 0, host_region->len(), host_rev);
    google::dense_hash_map<std::string, int> h_dists = get_distances(host_realignments, host_bp);
    breakpoint_t virus_bp = breakpoint_t(contig_map.id2name(virus_region->contig_id), 0, virus_region->len(), virus_rev);
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
                break;
            }
        }

        if (!has_prev_dup) {
            qnames_to_keep[qname_i] = std::vector<std::string>();
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

    std::string region_str = (vregion_rc ? "-" : "+") + print_region(virus_region);

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
            aligner.Align(read_seq_str.c_str(), chr_seqs.get_seq(contig_map.id2name(virus_region->contig_id)) + virus_region->start,
                          virus_region->len(), filter, &alignment, 0);

            if (accept_alignment(alignment, read_seq_str, config.min_sc_size)) {
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
            host_read_realignments_local, virus_read_realignments_local, host_region, hregion_rc, virus_region, !vregion_rc);

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

    std::vector<read_realignment_t>& host_read_realignments = host_region_score.rev ?
            reads_per_region_rc[host_region_score.region->id] : reads_per_region[host_region_score.region->id];

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

        if ((!host_region_score.rev && rr.right_clip_len())
        || (host_region_score.rev && rr.left_clip_len())) {
            remapped_host_anchors.push_back(rr);
        }
    }

    // extract virus regions
    std::vector<region_t*> virus_regions;
    extract_regions(remapped_reads_mates, virus_regions, max_is, contig_map, contig_tid2id, chr_seqs);
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
        region->end = std::min(region->end, (int) chr_seqs.get_len(contig_map.id2name(region->contig_id)));
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
            if (!host_region_score.rev && anchor.right_clip_len() && !anchor.rc) {
                set_to_forward(rc_read);
                rc_read->core.flag |= BAM_FREVERSE;
            }
            else if (host_region_score.rev && anchor.left_clip_len() && anchor.rc) {
                set_to_reverse(rc_read);
                rc_read->core.flag &= ~BAM_FREVERSE;
            }
            rc_read->core.flag |= BAM_FSECONDARY;

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
        std::pair<int, int> score = remap_virus_reads_supp(virus_region, false, read_seq_ids, host_region_score.region, host_region_score.rev,
                host_alignments, NULL, aligner, filter, temp1, temp2);
        std::pair<int, int> rc_score = remap_virus_reads_supp(virus_region, true, read_seq_ids, host_region_score.region, host_region_score.rev,
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
    std::cerr << "BEST VIRUS REGION: " << print_region(best_region) << " " << best_score.first << std::endl;
    remap_virus_reads_supp(best_region, !is_best_vregion_fwd, read_seq_ids, host_region_score.region, host_region_score.rev,
            host_alignments, &virus_read_realignments, aligner, filter, host_window, virus_window);

    google::dense_hash_map<std::string, std::vector<std::string> > deduped_reads = dedup_reads(
            host_alignments, virus_read_realignments, host_region_score.region, host_region_score.rev, best_region, is_best_vregion_fwd);

    std::cerr << "BEFORE " << host_region_score.to_string() << std::endl;
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
    std::cerr << "AFTER " << host_region_score.to_string() << std::endl << std::endl;

    return {best_region, is_best_vregion_fwd};
}


double mean(std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
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


int main(int argc, char* argv[]) {

    std::string host_reference_fname  = argv[1];
    std::string virus_reference_fname  = argv[2];
    std::string workdir = argv[3];
    std::vector<std::string> workspaces;
    for (int i = 4; i < argc; i++) {
        workspaces.push_back(argv[i]);
    }

    contig_map = contig_map_t(workdir);

    chr_seqs.read_fasta_into_map(host_reference_fname);
    chr_seqs.read_fasta_into_map(virus_reference_fname, true, true);

    config = parse_config(workdir + "/config.txt");
    for (std::string workspace : workspaces) {
        stats_t stats = parse_stats(workspace + "/stats.txt");
        max_is = std::max(max_is, stats.max_is);
    }

    open_samFile_t* host_and_virus_file = open_samFile((workspaces[0] + "/retained-pairs.remapped.cs.bam").data()); // just for the header
    for (int i = 0; i < host_and_virus_file->header->n_targets; i++) {
        int contig_id = contig_map.name2id(host_and_virus_file->header->target_name[i]);
        contig_id2tid[contig_id] = i;
        contig_tid2id[i] = contig_id;
    }


    /* === READ CHIMERIC PAIRS/READS === */

    // extracting paired host-virus reads
    std::vector<bam1_t*> host_reads, virus_reads;
    load_reads(workdir + "/host-side.cs.bam", host_reads, false);
    load_reads(workdir + "/virus-side.cs.bam", virus_reads, false);

    google::dense_hash_map<std::string, bam1_t*> virus_reads_by_name;
    virus_reads_by_name.set_empty_key("");
    for (int i = 0; i < virus_reads.size(); i++) {
        virus_reads_by_name[bam_get_qname(virus_reads[i])] = virus_reads[i];
    }

    std::cerr << "HOST READS: " << host_reads.size() << std::endl;
    std::cerr << "VIRUS READS: " << virus_reads.size() << std::endl;

    /* ====== */


    /* === EXTRACT REGIONS === */

    std::string host_regions_path = workdir + "/host-regions";
    std::vector<region_t*> host_regions;
    std::ifstream host_regions_fin(host_regions_path);

    int region_id, start, end;
    std::string contig_name;
    while (host_regions_fin >> region_id >> contig_name >> start >> end) {
        int contig_id = contig_map.name2id(contig_name);
        region_t* region = new region_t(contig_id, start, end);
        region->id = region_id;
        host_regions.push_back(region);
    }
    host_regions_fin.close();

    std::cerr << "HOST REGIONS: " << host_regions.size() << std::endl;

    /* ====== */


    /* === COMPUTE READS-REGIONS ALIGNMENTS === */

    reads_per_region = new std::vector<read_realignment_t>[host_regions.size()];
    reads_per_region_rc = new std::vector<read_realignment_t>[host_regions.size()];

    std::string scores_file_path = workdir + "/reads.scores.bin";
    FILE* scores_file_in_bin = fopen(scores_file_path.c_str(), "rb");
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
        add_realignment_to_region(rr, region_id, config.min_sc_size, reads_per_region, reads_per_region_rc);
    }
    fclose(scores_file_in_bin);
    std::cerr << "Mappings read." << std::endl;

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
            region_scores.push(region_score_t(region_score_id++, region, false, score, reads, reads));
        }

        score_and_reads = compute_region_score(region, true, already_used);
        score = score_and_reads.first;
        reads = score_and_reads.second;
        if (score > 0) {
            region_scores.push(region_score_t(region_score_id++, region, true, score, reads, reads));
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

            std::cerr << "CURRENT TOP REGION: " << best_region_score.to_string() << std::endl;

            // remove already used read_and_score
            std::vector<read_realignment_t>& host_read_realignments = best_region_score.rev ?
                    reads_per_region_rc[best_region_score.region->id] : reads_per_region[best_region_score.region->id];
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
                    best_region_score.rev, best_virus_region.first, best_virus_region.second);


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
                char* chr = chr_seqs.get_seq(contig_map.id2name(contig_tid2id[rr.read->core.tid]));
                bool lc_q = is_low_complexity(rr.read, rr.rc, rr.left_clip_len(), rr.read->core.l_qseq-rr.right_clip_len());
                bool lc_r = is_low_complexity(chr, false, rr.read->core.pos, bam_endpos(rr.read));
                if (lc_q || lc_r) {
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

            std::unordered_set<std::string> suspicious_virus_aln;
            for (read_realignment_t& vrr : kept_virus_read_realignments) {
                if (is_suspicious(vrr, virus_window)) {
                    vrr.suspicious = true;
                    suspicious_virus_aln.insert(bam_get_qname(vrr.read));
                }
            }
            for (read_realignment_t& hrr : kept_host_read_realignments) {
                if (is_suspicious(hrr, host_window) || suspicious_virus_aln.count(bam_get_qname(hrr.read))) {
                    hrr.suspicious = true;
                }
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

            std::cerr << "CONFIRMED TOP REGION: " << best_region_score.to_string() << std::endl;
            break;
        }

        if (region_scores.empty()) break;

        region_score_t best_region_score = region_scores.top();
        region_t* best_region = best_region_score.region;
        bool host_region_rev = best_region_score.rev;
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
        edit_remapped_reads(best_virus_region.first, kept_virus_read_realignments);

        int v_min_pos = 1000000000, v_max_pos = 0;
        std::vector<double> virus_score_ratios;
        for (read_realignment_t& rr : kept_virus_read_realignments) {
            virus_score_ratios.push_back(rr.score/double(rr.offset_end-rr.offset_start+1));
            v_min_pos = std::min(v_min_pos, (int) rr.read->core.pos);
            v_max_pos = std::max(v_max_pos, (int) bam_endpos(rr.read));
        }
        std::string virus_chr = contig_map.id2name(best_virus_region.first->contig_id);
        bool virus_rev = !best_virus_region.second;
        breakpoint_t virus_bp = breakpoint_t(virus_chr, v_min_pos, v_max_pos, virus_rev);

        // edit host reads and find host bp
        edit_remapped_reads(best_region, kept_host_read_realignments);

        int reads_w_dups = 0, unique_reads_w_dups = 0;
        std::vector<double> host_score_ratios;
        for (read_realignment_t& rr : kept_host_read_realignments) {
            host_score_ratios.push_back(rr.score/double(rr.offset_end-rr.offset_start+1));
            reads_w_dups += dedup_qnames[bam_get_qname(rr.read)].size();
            if (uniquely_aligned.count(bam_get_qname(rr.read))) {
                unique_reads_w_dups += dedup_qnames[bam_get_qname(rr.read)].size();
            }
        }

        int score = 0;
        int h_min_pos = 1000000000, h_max_pos = 0;
        for (read_realignment_t& rr : kept_host_read_realignments) {
            h_min_pos = std::min(h_min_pos, (int) rr.read->core.pos);
            h_max_pos = std::max(h_max_pos, (int) bam_endpos(rr.read));
            score += rr.score;
        }
        std::string host_chr = contig_map.id2name(best_region->contig_id);
        breakpoint_t host_bp = breakpoint_t(host_chr, h_min_pos, h_max_pos, host_region_rev);

        // remove pairs having IS too long
        google::dense_hash_map<std::string, int> h_dists = get_distances(kept_host_read_realignments, host_bp, *best_region);
        google::dense_hash_map<std::string, int> v_dists = get_distances(kept_virus_read_realignments, virus_bp, *(best_virus_region.first));
        auto is_too_long = [&h_dists, &v_dists](const read_realignment_t& rr) {
            std::string qname = bam_get_qname(rr.read);
            return h_dists[qname] + v_dists[qname] > max_is;
        };
        for (read_realignment_t& rr : kept_host_read_realignments) {
            if (is_too_long(rr)) {
                rr.too_long = true;
//                bam_aux_update_str(rr.read, "CT", 2, "TL");
            }
            if (rr.suspicious) {
//                bam_aux_update_str(rr.read, "CT", 2, "SS");
            }
        }
        for (read_realignment_t& rr : kept_virus_read_realignments) {
            if (is_too_long(rr)) {
                rr.too_long = true;
//                bam_aux_update_str(rr.read, "CT", 2, "TL");
            }
            if (rr.suspicious) {
//                bam_aux_update_str(rr.read, "CT", 2, "SS");
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

        double host_coverage = calculate_coverage(kept_host_read_realignments)/double(max_is);
        double virus_coverage = calculate_coverage(kept_virus_read_realignments)/double(max_is);

        int split_reads = 0;
        for (read_realignment_t& rr : kept_host_read_realignments) {
            if (!host_bp.rev && rr.right_clip_len() >= config.min_sc_size
                && abs((host_bp.end - best_region_score.region->start) - rr.offset_end) < config.max_sc_dist) {
                split_reads++;
            } else if (host_bp.rev && rr.left_clip_len() >= config.min_sc_size
                       && abs((host_bp.start - best_region_score.region->start) - rr.offset_start) < config.max_sc_dist) {
                split_reads++;
            }
        }

        virus_bp.start %= chr_seqs.get_original_len(virus_bp.chr);
        virus_bp.end %= chr_seqs.get_original_len(virus_bp.chr);
        call_t v_int(virus_integration_id++,host_bp, virus_bp, best_region_score.reads,
                                  best_region_score.good_reads, split_reads,
                                  reads_w_dups, unique_reads_w_dups, best_region_score.score,
                                  mean(host_score_ratios), mean(virus_score_ratios), host_coverage, virus_coverage);
        std::cout << v_int.to_string() << std::endl;

        samFile* writer = open_bam_writer(workdir+"/readsx", std::to_string(v_int.id)+".bam", host_and_virus_file->header);
        for (read_realignment_t& rr : kept_virus_read_realignments) {
            // we duplicated viruses to deal with breakpoints across the border - scale back those alignment to original virus
            int virus_len = chr_seqs.get_original_len(virus_bp.chr);
            rr.read->core.pos %= virus_len;

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

    for (region_t* region : host_regions) delete region;
    for (auto& e : cid_to_cigar) {
        delete[] e.second;
    }
    delete[] reads_per_region;
    delete[] reads_per_region_rc;

    /* ==== */

}
