#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "config.h"
#include "libs/cptl_stl.h"

config_t config;
std::unordered_map<std::string, int> contig_name2id;
std::vector<int> tid_to_contig_id;
std::unordered_set<std::string> virus_names;

std::mutex mtx, mtx_lc, mtx_rc;

const int MAX_BUFFER_SIZE = 100;


bam1_t* add_to_queue(std::deque<bam1_t*>& q, bam1_t* o, int size_limit) {
    bam1_t* t = NULL;
    while (q.size() >= size_limit) {
        t = q.front();
        q.pop_front();
    }
    q.push_back(o);
    return t;
}

void extract(int id, std::string contig, std::string bam_fname, int target_len, std::vector<bam1_t *> *lc_anchors,
             std::vector<bam1_t *> *rc_anchors) {

    open_samFile_t* bam_file_open = open_samFile(bam_fname.data());

    samFile* bam_file = bam_file_open->file;
    hts_idx_t* idx = bam_file_open->idx;
    bam_hdr_t* header = bam_file_open->header;

    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), 1, target_len);

    mtx.lock();
    std::cout << "Extracting clips for " << region << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(idx, header, region);
    bam1_t* read = bam_init1();

    int i = 0;
    while (sam_itr_next(bam_file, iter, read) >= 0) {
        if (is_unmapped(read)) continue;

        // clipped read
        if (is_left_clipped(read) && get_left_clip_len(read) >= config.min_sc_size) {
            mtx_lc.lock();
            lc_anchors->push_back(bam_dup1(read));
            mtx_lc.unlock();
        } else if (is_right_clipped(read) && get_right_clip_len(read) >= config.min_sc_size) {
            mtx_rc.lock();
            rc_anchors->push_back(bam_dup1(read));
            mtx_rc.unlock();
        }
    }
    bam_destroy1(read);
    bam_itr_destroy(iter);

    close_samFile(bam_file_open);
}

int main(int argc, char* argv[]) {

    std::string virus_name;
    FILE* virus_ref_fasta = fopen(argv[1], "r");
    kseq_t *seq = kseq_init(fileno(virus_ref_fasta));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        virus_names.insert(seq->name.s);
    }
    kseq_destroy(seq);
    fclose(virus_ref_fasta);

    std::string workdir = argv[2];
    std::string workspace = argv[3];
    std::string bam_fname = workspace + "/retained-pairs-remapped.sorted.bam";

    config = parse_config(workdir + "/config.txt");

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_name2id[contig_name] = contig_id;
    }

    ctpl::thread_pool thread_pool(config.threads);


    open_samFile_t* bam_file = open_samFile(bam_fname.c_str(), true);
    bam_hdr_t* header = bam_file->header;
    tid_to_contig_id.resize(header->n_targets);

    std::vector<bam1_t*> lc_virus_anchors, rc_virus_anchors, lc_host_anchors, rc_host_anchors;

    std::vector<std::future<void> > futures;
    for (int i = 0; i < header->n_targets; i++) {
        bool is_virus = virus_names.count(header->target_name[i]);
        std::future<void> future = thread_pool.push(extract, header->target_name[i], bam_fname, header->target_len[i],
                                                    is_virus ? &lc_virus_anchors : &lc_host_anchors,
                                                    is_virus ? &rc_virus_anchors : &rc_host_anchors);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cout << s << std::endl;
        }
    }

    samFile* virus_anchor_writer = open_bam_writer(workspace, "virus-anchors.bam", header);
    std::ofstream virus_clip_out(workspace + "/virus-clips.fa");
    for (bam1_t* anchor : lc_virus_anchors) {
        int ok = sam_write1(virus_anchor_writer, header, anchor);
        if (ok < 0) throw "Failed to write to " + std::string(virus_anchor_writer->fn);
        virus_clip_out << ">" << bam_get_qname(anchor) << "_L_" << (anchor->core.flag & BAM_FREAD1 ? "1" : "2") << "\n";
        virus_clip_out << get_left_clip(anchor) << "\n";
        bam_destroy1(anchor);
    }
    for (bam1_t* anchor : rc_virus_anchors) {
        int ok = sam_write1(virus_anchor_writer, header, anchor);
        if (ok < 0) throw "Failed to write to " + std::string(virus_anchor_writer->fn);
        virus_clip_out << ">" << bam_get_qname(anchor) << "_R_" << (anchor->core.flag & BAM_FREAD1 ? "1" : "2") << "\n";
        virus_clip_out << get_right_clip(anchor) << "\n";
        bam_destroy1(anchor);
    }
    sam_close(virus_anchor_writer);
    virus_clip_out.close();

    samFile* host_anchor_writer = open_bam_writer(workspace, "host-anchors.bam", header);
    std::ofstream host_clip_out(workspace + "/host-clips.fa");
    for (bam1_t* anchor : lc_host_anchors) {
        int ok = sam_write1(host_anchor_writer, header, anchor);
        if (ok < 0) throw "Failed to write to " + std::string(host_anchor_writer->fn);
        host_clip_out << ">" << bam_get_qname(anchor) << "_L_" << (anchor->core.flag & BAM_FREAD1 ? "1" : "2") << "\n";
        host_clip_out << get_left_clip(anchor) << "\n";
        bam_destroy1(anchor);
    }
    for (bam1_t* anchor : rc_host_anchors) {
        int ok = sam_write1(host_anchor_writer, header, anchor);
        if (ok < 0) throw "Failed to write to " + std::string(host_anchor_writer->fn);
        host_clip_out << ">" << bam_get_qname(anchor) << "_R_" << (anchor->core.flag & BAM_FREAD1 ? "1" : "2") << "\n";
        host_clip_out << get_right_clip(anchor) << "\n";
        bam_destroy1(anchor);
    }
    sam_close(host_anchor_writer);
    host_clip_out.close();

    close_samFile(bam_file);
}
