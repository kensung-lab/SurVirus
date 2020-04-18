#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <htslib/sam.h>

#include "sam_utils.h"
#include "utils.h"

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

contig_map_t contig_map;
std::unordered_map<int, int> contig_id2tid;

int fix_tid(open_samFile_t* sam_file, bam1_t* read) {
    std::string target_name = sam_file->header->target_name[read->core.tid];
    return contig_id2tid[contig_map.name2id(target_name)];
}

std::vector<std::pair<bam1_t*, bam1_t*> > compose_read_pairs(std::vector<bam1_t*>& host_reads, std::vector<bam1_t*>& virus_reads) {
    std::vector<bam1_t*> primary_reads;
    for (bam1_t* read : host_reads) {
        if (is_primary(read)) {
            primary_reads.push_back(read);
        }
    }
    for (bam1_t* read : virus_reads) {
        if (is_primary(read)) {
            primary_reads.push_back(read);
        }
    }

    std::unordered_map<std::string, std::pair<bam1_t*, bam1_t*> > read_pairs_map;
    for (bam1_t* read : primary_reads) {
        std::string qname = bam_get_qname(read);
        if (!read_pairs_map.count(qname)) read_pairs_map[qname] = {read, NULL};
        else read_pairs_map[qname].second = read;
    }

    std::vector<std::pair<bam1_t*, bam1_t*> > read_pairs;
    for (auto& e : read_pairs_map) {
        read_pairs.push_back(e.second);
    }
    return read_pairs;
}

void mark_duplicates(std::vector<bam1_t*>& host_reads, std::vector<bam1_t*>& virus_reads) {

    if (host_reads.empty() || virus_reads.empty()) return;

    std::vector<std::pair<bam1_t*, bam1_t*> > read_pairs = compose_read_pairs(host_reads, virus_reads);

    // order by host and virus pos, then cigar code, and finally by average base quality
    auto pair_cmp = [] (const std::pair<bam1_t*,bam1_t*>& p1, const std::pair<bam1_t*,bam1_t*>& p2) {
        bam1_t* h1 = p1.first, * h2 = p2.first;
        bam1_t* v1 = p1.second, * v2 = p2.second;
        std::string h1_cigar = get_cigar_code(h1), v1_cigar = get_cigar_code(v1);
        std::string h2_cigar = get_cigar_code(h2), v2_cigar = get_cigar_code(v2);
        int avg_qual1 = get_avg_qual(h1, false)+get_avg_qual(v1, false);
        int avg_qual2 = get_avg_qual(h2, false)+get_avg_qual(v2, false);
        return std::tie(h1->core.tid, h1->core.pos, h1_cigar, v1->core.tid, v1->core.pos, v1_cigar, avg_qual1) <
               std::tie(h2->core.tid, h2->core.pos, h2_cigar, v2->core.tid, v2->core.pos, v2_cigar, avg_qual2);
    };
    std::sort(read_pairs.begin(), read_pairs.end(), pair_cmp);


    // MARK DUPS BY SEQUENCE
    // "head" is the pair not being marked as duplicate
    google::dense_hash_map<std::string, std::string> qname_head_by_joint_seq, head_by_qname;
    qname_head_by_joint_seq.set_empty_key("");
    head_by_qname.set_empty_key("");
    for (std::pair<bam1_t*, bam1_t*>& read_pair : read_pairs) {
        bam1_t* read1 = read_pair.first, * read2 = read_pair.second;
        std::string seq1 = get_sequence(read1, true), seq2 = get_sequence(read2, true);
        std::string joint_seq = std::min(seq1, seq2) + "-" + std::max(seq1, seq2);
        if (qname_head_by_joint_seq.count(joint_seq)) {
            head_by_qname[bam_get_qname(read1)] = qname_head_by_joint_seq[joint_seq];
        } else {
            qname_head_by_joint_seq[joint_seq] = bam_get_qname(read1);
        }
    }

    // MARK DUPS BY POSITION AND CIGAR
    int last_valid_i = 0;
    for (int i = 1; i < read_pairs.size(); i++) {
        bam1_t* h1 = read_pairs[last_valid_i].first, * h2 = read_pairs[i].first;
        bam1_t* v1 = read_pairs[last_valid_i].second, * v2 = read_pairs[i].second;
        if (head_by_qname.count(bam_get_qname(h2))) continue;

        std::string h1_cigar = get_cigar_code(h1), v1_cigar = get_cigar_code(v1);
        std::string h2_cigar = get_cigar_code(h2), v2_cigar = get_cigar_code(v2);

        if (std::tie(h1->core.tid, h1->core.pos, h1_cigar, v1->core.tid, v1->core.pos, v1_cigar) ==
            std::tie(h2->core.tid, h2->core.pos, h2_cigar, v2->core.tid, v2->core.pos, v2_cigar)) {
            head_by_qname[bam_get_qname(h2)] = bam_get_qname(h1);
        } else {
            last_valid_i = i;
        }
    }

    for (bam1_t* h_read : host_reads) {
        std::string qname = bam_get_qname(h_read);
        if (head_by_qname.count(qname)) {
            h_read->core.flag |= BAM_FDUP;
            // this corrupts the BAM file. Am I using it incorrectly or is there a bug in the library?
//            bam_aux_update_str(h_read, "DH", qname.length(), qname.c_str());
        }
    }
    for (bam1_t* v_read : virus_reads) {
        std::string qname = bam_get_qname(v_read);
        if (head_by_qname.count(bam_get_qname(v_read))) {
            v_read->core.flag |= BAM_FDUP;
            // this corrupts the BAM file. Am I using it incorrectly or is there a bug in the library?
//            bam_aux_update_str(v_read, "DH", qname.length(), qname.c_str());
        }
    }
}

int main(int argc, char* argv[]) {

    std::string workdir = argv[1];
    std::vector<std::string> workspaces, bam_fnames;
    for (int i = 2; i < argc; i+=2) {
        bam_fnames.push_back(argv[i]);
        workspaces.push_back(argv[i+1]);
    }

    contig_map = contig_map_t(workdir);
    open_samFile_t* host_and_virus_file = open_samFile((workspaces[0] + "/retained-pairs.remapped.cs.bam").data()); // just for the header
    for (int i = 0; i < host_and_virus_file->header->n_targets; i++) {
        int contig_id = contig_map.name2id(host_and_virus_file->header->target_name[i]);
        contig_id2tid[contig_id] = i;
    }
    close_samFile(host_and_virus_file);

    bam1_t* read = bam_init1();


    // extracting paired host-virus reads
    std::vector<bam1_t*> host_reads, virus_reads;
    for (int i = 0; i < workspaces.size(); i++) {
        std::string workspace = workspaces[i];

        std::unordered_set<std::string> host_qnames;
        open_samFile_t* host_file = open_samFile((workspace + "/host-side.cs.bam").data(), true);
        while (sam_read1(host_file->file, host_file->header, read) >= 0) {
            bam_aux_update_int(read, "BF", i);
            host_qnames.insert(bam_get_qname(read));
            host_reads.push_back(bam_dup1(read));
        }
        close_samFile(host_file);

        std::unordered_set<std::string> virus_qnames;
        open_samFile_t* virus_file = open_samFile((workspace + "/virus-side.cs.bam").data(), true);
        while (sam_read1(virus_file->file, virus_file->header, read) >= 0) {
            bam_aux_update_int(read, "BF", i);
            virus_qnames.insert(bam_get_qname(read));
            virus_reads.push_back(bam_dup1(read));
        }
        close_samFile(virus_file);


        // pos of good clips for first (resp. second) reads
        std::unordered_map<std::string, std::pair<int, int> > good_clips_pos[2];
        open_samFile_t* host_clips_file = open_samFile((workspace + "/host-clips.cs.bam").data(), true);
        while (sam_read1(host_clips_file->file, host_clips_file->header, read) >= 0) {
            std::string qname = bam_get_qname(read);
            bool first = qname[qname.length()-1] == '1';
            qname = qname.substr(0, qname.length()-4);
            good_clips_pos[first ? 0 : 1][qname] = {fix_tid(host_clips_file, read), read->core.pos};
        }
        close_samFile(host_clips_file);

        open_samFile_t* virus_clips_file = open_samFile((workspace + "/virus-clips.cs.bam").data(), true);
        while (sam_read1(virus_clips_file->file, virus_clips_file->header, read) >= 0) {
            std::string qname = bam_get_qname(read);
            bool first = qname[qname.length()-1] == '1';
            qname = qname.substr(0, qname.length()-4);
            good_clips_pos[first ? 0 : 1][qname] = {fix_tid(virus_clips_file, read), read->core.pos};
        }
        close_samFile(virus_clips_file);


        // pairs where both are good clipped reads. Bool = true indicates a host read
        std::unordered_map<std::string, std::vector<std::pair<bam1_t*, bool> > > good_clipped_pairs;

        /* == Reading host anchors == */

        open_samFile_t* host_anchors_file = open_samFile((workspace + "/host-anchors.cs.bam").data(), true);
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

        open_samFile_t* virus_anchors_file = open_samFile((workspace + "/virus-anchors.cs.bam").data(), true);
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
            } else { // destroy unutilised reads
                for (std::pair<bam1_t*, bool> p : e.second) {
                    bam_destroy1(p.first);
                }
            }
        }
    }

    bam_destroy1(read);

    // remove pairs not having exactly 2 primary alignments
    std::unordered_map<std::string, int> qnames_count;
    for (bam1_t* r : host_reads) {
        if (is_primary(r)) qnames_count[bam_get_qname(r)]++;
    }
    for (bam1_t* r : virus_reads) {
        if (is_primary(r)) qnames_count[bam_get_qname(r)]++;
    }

    auto not_2_primary_aln = [&qnames_count](const bam1_t* r) { return qnames_count[bam_get_qname(r)] != 2; };
    host_reads.erase(std::remove_if(host_reads.begin(), host_reads.end(), not_2_primary_aln), host_reads.end());
    virus_reads.erase(std::remove_if(virus_reads.begin(), virus_reads.end(), not_2_primary_aln), virus_reads.end());


    mark_duplicates(host_reads, virus_reads);


    open_samFile_t* host_template_header = open_samFile((workspaces[0] + "/host-side.cs.bam").c_str(), true);
    samFile* host_writer = open_bam_writer(workdir, "host-side.cs.bam", host_template_header->header);
    open_samFile_t* virus_template_header = open_samFile((workspaces[0] + "/virus-side.cs.bam").c_str(), true);
    samFile* virus_writer = open_bam_writer(workdir, "virus-side.cs.bam", virus_template_header->header);

    for (bam1_t* read : host_reads) {
        int ok = sam_write1(host_writer, host_template_header->header, read);
        if (ok < 0) throw "Failed to write to " + std::string(host_writer->fn);
        bam_destroy1(read);
    }
    for (bam1_t* read : virus_reads) {
        int ok = sam_write1(virus_writer, virus_template_header->header, read);
        if (ok < 0) throw "Failed to write to " + std::string(virus_writer->fn);
        bam_destroy1(read);
    }

    sam_close(host_writer);
    sam_close(virus_writer);

    close_samFile(host_template_header);
    close_samFile(virus_template_header);

}