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

std::mutex mtx;

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

void categorize(int id, std::string contig, std::string bam_fname, int target_len,
                std::vector<bam1_t*>* reads, std::vector<bam1_t*>* anchor_reads,
                std::unordered_map<std::string, std::pair<char, char> >& good_clips) {

    open_samFile_t* bam_file_open = open_samFile(bam_fname.data());

    samFile* bam_file = bam_file_open->file;
    hts_idx_t* idx = bam_file_open->idx;
    bam_hdr_t* header = bam_file_open->header;

    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), 1, target_len);

    mtx.lock();
    std::cout << "Categorizing " << region << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(idx, header, region);
    bam1_t* read = bam_init1();

    int i = 0;
    while (sam_itr_next(bam_file, iter, read) >= 0) {
        if (!is_unmapped(read)) {
            if (good_clips.count(bam_get_qname(read))) {
                std::pair<char, char> good_clip_flag = good_clips[bam_get_qname(read)];
                bool is_first_read = read->core.flag & BAM_FREAD1;
                mtx.lock();
                if ((is_first_read && good_clip_flag.first) || (!is_first_read && good_clip_flag.second)) {
                    uint8_t dir = is_first_read ? good_clip_flag.first : good_clip_flag.second;
                    bam1_t* copy = bam_dup1(read);
                    bam_aux_append(copy, "CD", 'A', 1, &dir);
                    anchor_reads->push_back(copy);
                } else {
                    reads->push_back(bam_dup1(read));
                }
                mtx.unlock();
            } else if (is_dc_pair(read) && !is_poly_ACGT(read, true)) {
                std::string target_name = contig;
                std::string mate_target_name = header->target_name[read->core.mtid];
                if (virus_names.count(target_name) != virus_names.count(mate_target_name)) {
                    mtx.lock();
                    reads->push_back(bam_dup1(read));
                    mtx.unlock();
                }
            }
        }
    }

    close_samFile(bam_file_open);
    bam_destroy1(read);
    bam_itr_destroy(iter);
}

std::string print_fq(bam1_t* r) {
    std::stringstream ss;
    ss << "@" << bam_get_qname(r) << "\n";
    std::string seq = get_sequence(r);
    if (bam_is_rev(r)) get_rc(seq);
    ss << seq << "\n";
    ss << "+" << "\n";
    for (int i = 0; i < r->core.l_qseq; i++) {
        int idx = bam_is_rev(r) ? r->core.l_qseq-i-1 : i;
        ss << (char) (bam_get_qual(r)[idx] + 33);
    }
    ss << std::endl;
    return ss.str();
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

    std::unordered_map<std::string, std::pair<char, char> > good_clips;
    bam1_t* read = bam_init1();

    open_samFile_t* host_clips_file = open_samFile((workspace + "/host-clips.sorted.bam").data(), true);
    while (sam_read1(host_clips_file->file, host_clips_file->header, read) >= 0) {
        if (virus_names.count(host_clips_file->header->target_name[read->core.tid]) ) { // seq was mapped to virus
            std::string clip_name = bam_get_qname(read);
            std::string qname = clip_name.substr(0, clip_name.length()-4);
            char dir = clip_name[clip_name.length()-3];
            if (clip_name[clip_name.length()-1] == '1') {
                good_clips[qname].first = dir;
            } else {
                good_clips[qname].second = dir;
            }
        }
    }
    close_samFile(host_clips_file);

    open_samFile_t* virus_clips_file = open_samFile((workspace + "/virus-clips.sorted.bam").data(), true);
    while (sam_read1(virus_clips_file->file, virus_clips_file->header, read) >= 0) {
        if (!virus_names.count(virus_clips_file->header->target_name[read->core.tid]) ) { // seq was mapped to host
            std::string clip_name = bam_get_qname(read);
            std::string qname = clip_name.substr(0, clip_name.length()-4);
            char dir = clip_name[clip_name.length()-3];
            if (clip_name[clip_name.length()-1] == '1') {
                good_clips[qname].first = dir;
            } else {
                good_clips[qname].second = dir;
            }
        }
    }
    close_samFile(virus_clips_file);
    bam_destroy1(read);

    std::vector<bam1_t*> virus_side_reads, host_side_reads;
    std::vector<bam1_t*> virus_anchors, host_anchors;

    std::vector<std::future<void> > futures;
    for (int i = 0; i < header->n_targets; i++) {
        bool is_virus = virus_names.count(header->target_name[i]);
        std::future<void> future = thread_pool.push(categorize, header->target_name[i], bam_fname, header->target_len[i],
                                                    is_virus ? &virus_side_reads : &host_side_reads,
                                                    is_virus ? &virus_anchors : &host_anchors, good_clips);
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

    std::ofstream virus_fq(workspace + "/virus-side.fq"), host_fq(workspace + "/host-side.fq");
    for (bam1_t* r : host_side_reads) {
        host_fq << print_fq(r);
        bam_destroy1(r);
    }
    for (bam1_t* r : virus_side_reads) {
        virus_fq << print_fq(r);
        bam_destroy1(r);
    }

    samFile* virus_anchor_writer = open_bam_writer(workspace, "virus-anchors.bam", header);
    for (bam1_t* anchor : virus_anchors) {
        int ok = sam_write1(virus_anchor_writer, header, anchor);
        if (ok < 0) throw "Failed to write to " + std::string(virus_anchor_writer->fn);
        bam_destroy1(anchor);
    }
    sam_close(virus_anchor_writer);

    samFile* host_anchor_writer = open_bam_writer(workspace, "host-anchors.bam", header);
    for (bam1_t* anchor : host_anchors) {
        int ok = sam_write1(host_anchor_writer, header, anchor);
        if (ok < 0) throw "Failed to write to " + std::string(host_anchor_writer->fn);
        bam_destroy1(anchor);
    }
    sam_close(host_anchor_writer);

    close_samFile(bam_file);
}
