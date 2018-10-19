#include <iostream>
#include <fstream>
#include <unordered_set>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "sam_utils.h"
#include "config.h"
#include "libs/cptl_stl.h"

std::unordered_set<std::string> qnames;
samFile* bam_writer;
std::mutex mtx;

void filter_by_qname(int id, char* bam_fname, int t_num) {

    open_samFile_t* bam_file = open_samFile(bam_fname, false);

    char region[1000];
    sprintf(region, "%s:%d-%d", bam_file->header->target_name[t_num], 1, bam_file->header->target_len[t_num]);

    mtx.lock();
    std::cerr << "Filtering " << bam_file->header->target_name[t_num] << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region);

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (qnames.count(std::string(bam_get_qname(read)))) {
            mtx.lock();
            int ok = sam_write1(bam_writer, bam_file->header, read);
            mtx.unlock();
            if (ok < 0) throw "Failed to write to " + std::string(bam_writer->fn);
        }
    }
    bam_destroy1(read);

    sam_itr_destroy(iter);

    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];

    // load qnames
    std::ifstream qnames_fin(workspace + "/qnames-to-keep");
    std::string qname;
    while (qnames_fin >> qname) {
        qnames.insert(qname);
    }

    config_t config = parse_config(workdir + "/config.txt");

    ctpl::thread_pool thread_pool(config.threads);

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());
    bam_writer = open_bam_writer(workspace, "retained-pairs.bam", bam_file->header);

    std::vector<std::future<void> > futures;
    for (int t_num = 0; t_num < bam_file->header->n_targets; t_num++) {
        std::future<void> future = thread_pool.push(filter_by_qname, bam_file->file->fn, t_num);
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

    close_samFile(bam_file);
    sam_close(bam_writer);

}