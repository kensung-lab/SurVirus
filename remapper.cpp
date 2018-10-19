#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <cassert>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <curses.h>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "cluster.h"
#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"

config_t config;
stats_t stats;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;
std::unordered_map<int, int> contig_id2tid;
std::unordered_map<std::string, std::pair<char*, size_t> > chrs;

const int SKIP_READ = -1;


struct virus_integration_t {
    int contig_id;
    int pos;
    char dir;
    int reads;

    virus_integration_t(int contig_id, int pos, char dir, int reads) : contig_id(contig_id), pos(pos), dir(dir), reads(reads) {}
};


struct region_t {
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    int score = 0;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}

    std::string to_str() {
        return contig_id2name[contig_id] + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    int len() {
        return end-start+1;
    }
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

std::vector<region_t> extract_regions(std::vector<bam1_t *> &reads, bam_hdr_t* header, bool process_xa = true) {

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

//    for (cluster_t* c : clip_clusters) {
//        c->id = -1;
//        clusters.push_back(c);
//        clusters_map.insert(std::make_pair(c->a1.start, c));
//        clusters_map.insert(std::make_pair(c->a1.end, c));
//    }


    std::vector<int> max_dists;
//    for (int i = 0; i < 10; i++) max_dists.push_back(i);
//    for (int i = 10; i < 100; i+=10) max_dists.push_back(i);
//    for (int i = 100; i < config.max_is; i+=100) max_dists.push_back(i);
    max_dists.push_back(stats.max_is);

    for (int max_dist : max_dists) {
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


    std::vector<region_t> regions;

    for (cluster_t* c : clusters) {
        if (c->dead) continue;
        regions.push_back(region_t(c->a1.contig_id, contig_id2tid[c->a1.contig_id],
                                   std::max(0, c->a1.start-stats.max_is),
                                   std::min(c->a1.end+stats.max_is, (int) chrs[contig_id2name[c->a1.contig_id]].second)));
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
//
    return regions;
}

int compute_score_supp(region_t &region, std::vector<bam1_t *> &reads, bool do_rc, std::vector<int>* scores,
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
            rc(s);
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


region_t* remap(std::vector<region_t>& regions, std::vector<bam1_t*>& reads, std::vector<std::pair<bam1_t*, int> >& remapped_now,
           std::vector<std::pair<bam1_t*, int> >& not_remapped) {

    StripedSmithWaterman::Aligner aligner(2, 2, 3, 1, false);
    StripedSmithWaterman::Aligner aligner_to_base(2, 2, 3, 1, true);
    StripedSmithWaterman::Filter filter, filter_w_cigar;
    filter_w_cigar.report_cigar = true;

    bool is_rc;
    for (region_t& region : regions) {
        compute_score(region, reads, NULL, NULL, NULL, is_rc, aligner, filter);
    }
    std::sort(regions.begin(), regions.end(), [](const region_t& r1, const region_t& r2) {return r1.score > r2.score;});

    region_t* best_region = &regions[0];

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


int main(int argc, char* argv[]) {
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

//    for (bam1_t* r : reference_reads) {
//        std::cerr << bam_get_qname(r) << " " << get_sequence(r) << std::endl;
//    }
//    exit(0);



    stats.max_is *= 2; // since we are clustering reads in two directions
    std::vector<region_t> reference_regions = extract_regions(reference_reads, reference_file->header);
    std::vector<region_t> virus_regions = extract_regions(virus_reads, virus_file->header, false);


//    for (std::vector<bam1_t*>& cluster : clusters) {
//        std::sort(cluster.begin(), cluster.end(), [](const bam1_t* r1, const bam1_t* r2) {return r1->core.pos < r2->core.pos;});
//        bam1_t* fr = *(cluster.begin());
//        bam1_t* lr = *(cluster.rbegin());
//        int contig_id = contig_name2id[reference_file->header->target_name[fr->core.tid]];
//        regions.push_back(region_t(contig_id, fr->core.tid, fr->core.pos-config.max_is, lr->core.pos+config.max_is));
//    }

    std::sort(reference_regions.begin(), reference_regions.end(), [](const region_t& reg1, const region_t& reg2) {
        if (reg1.contig_id == reg2.contig_id) return reg1.start < reg2.start;
        return reg1.contig_id < reg2.contig_id;
    });
    for (region_t region : reference_regions) {
        std::cout << region.to_str() << std::endl;
    }
    std::cout << std::endl;


    std::ofstream result_out(workspace + "/result.txt");

    int id = 1;
    int prev = 0;
    std::vector<bam1_t*> remaining = reference_reads;
    std::vector<bam1_t*> remapped;
    bool has_been_remapped = false;
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

        region_t* best_reference_region = remap(reference_regions, remaining, remapped_now, not_remapped);

        int forward = 0, forward_pos = -1, reverse = 0, reverse_pos = INT32_MAX;
        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
            bam1_t* read = read_and_score.first;
            if (bam_is_rev(read)) {
                reverse_pos = std::min(reverse_pos, read->core.pos);
                reverse++;
            } else {
                forward_pos = std::max(forward_pos, bam_endpos(read));
                forward++;
            }
        }
        if (forward_pos == -1) {
            forward_pos = reverse_pos;
        } else if (reverse_pos == INT32_MAX) {
            reverse_pos = forward_pos;
        }


        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
            remapped.push_back(read_and_score.first);
        }

        // extract virus reads
        std::unordered_set<std::string> remapped_qnames;
        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
            remapped_qnames.insert(bam_get_qname(read_and_score.first));
        }
        std::vector<bam1_t*> remapped_now_virus_reads;
        for (bam1_t* r : virus_reads) {
            if (remapped_qnames.count(bam_get_qname(r))) {
                remapped_now_virus_reads.push_back(r);
            }
        }
//        for (region_t& region : virus_regions) {
//            compute_score(region, remapped_now_virus_reads, NULL, NULL, is_rc, aligner, filter);
//        }

        std::vector<std::pair<bam1_t*, int> > v1, v2;
        region_t* best_virus_region = remap(virus_regions, remapped_now_virus_reads, v1, v2);

        int tot_score = 0;
        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
            tot_score += read_and_score.second;
        }

        remaining.clear();
        for (std::pair<bam1_t*, int> read_and_score : not_remapped) {
            remaining.push_back(read_and_score.first);
        }


        auto distinct_qnames = [](std::vector<std::pair<bam1_t*, int> >& reads_and_scores) {
            std::set<std::string> s;
            for (std::pair<bam1_t*, int> read_and_score : reads_and_scores) {
                std::string qname = bam_get_qname(read_and_score.first);
                s.insert(qname);
            }
            return s.size();
        };

        char dir = forward > 0 ? 'R' : 'L';
        virus_integration_t v_int(best_reference_region->contig_id, dir == 'R' ? forward_pos : reverse_pos, dir, distinct_qnames(remapped_now));

        std::stringstream ss;
        ss << contig_id2name[v_int.contig_id] << " " << v_int.pos << " " << v_int.pos
           << " " << v_int.reads << " " << forward << " " << reverse << " " << tot_score << std::endl;
        std::string int_s = ss.str();

        result_out << int_s << std::endl;

        std::cout << "ID=" << id++ << " " << int_s << std::endl;
        std::cout << "VIRUS " << best_virus_region->to_str() << std::endl;
        for (std::pair<bam1_t*, int> read_and_score : remapped_now) {
            bam1_t* r = read_and_score.first;
            std::cout << bam_get_qname(r) << "\t" << get_sequence(r) << "\t" << r->core.pos << "\t" << read_and_score.second << std::endl;
        }
        for (region_t& region : reference_regions) {
            std::cout << region.to_str() << " " << region.score << std::endl;
        }
        for (bam1_t* r : remapped_now_virus_reads) {
            std::cout << bam_get_qname(r) << " " << get_sequence(r, true) << std::endl;
        }

        std::cout << std::endl;
        prev = remapped.size();

        auto it = std::find_if(reference_regions.begin(), reference_regions.end(), [best_reference_region](region_t& elem) {return &elem == best_reference_region;});
        reference_regions.erase(it);

        has_been_remapped = !remapped_now.empty();
    } while (has_been_remapped && !remaining.empty());


    samFile* reference_remapped_file = open_bam_writer(workspace, "reference-remapped.bam", reference_file->header);
    for (bam1_t* read : remapped) {
        int ok = sam_write1(reference_remapped_file, reference_file->header, read);
        if (ok < 0) throw "Failed to write to " + std::string(reference_remapped_file->fn);
    }

    sam_close(reference_remapped_file);

    close_samFile(reference_file);
    result_out.close();
}

