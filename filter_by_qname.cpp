#include <iostream>
#include <unordered_set>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "sam_utils.h"
#include "config.h"
#include "libs/cptl_stl.h"

config_t config;

struct fastq_entry_t {
    std::string qname;
    std::string seq, qual;

    fastq_entry_t(bam1_t* read) {
        qname = bam_get_qname(read);
        seq = get_sequence(read, true);
        qual = get_qualities(read, true);
    }
};

std::unordered_map<std::string, std::pair<fastq_entry_t*, fastq_entry_t*> > qnames;
std::ofstream fq1, fq2;
std::mutex mtx;

void write_fq(std::ofstream& fout, fastq_entry_t* fq_entry) {
    fout << "@" << fq_entry->qname << "\n";
    fout << fq_entry->seq << "\n+\n";
    fout << fq_entry->qual << std::endl;
}

void filter_by_qname(int id, char* bam_fname, char* contig) {

    open_samFile_t* bam_file = open_samFile(bam_fname, false);
    if (!config.cram_reference.empty()) {
        if (hts_set_fai_filename(bam_file->file, fai_path(config.cram_reference.c_str())) != 0) {
            throw "Failed to read reference " + config.cram_reference;
        }
    }

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig);

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (!is_primary(read)) continue;

        if (qnames.count(bam_get_qname(read))) {
            mtx.lock();
            auto& e = qnames[bam_get_qname(read)];
            if (read->core.flag & BAM_FREAD1) {
                e.first = new fastq_entry_t(read);
            } else {
                e.second = new fastq_entry_t(read);
            }

            if (e.first && e.second) {
                write_fq(fq1, e.first);
                write_fq(fq2, e.second);
                delete e.first;
                delete e.second;
                e.first = e.second = NULL;
            }
            mtx.unlock();
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
        qnames[qname] = {NULL, NULL};
    }

    fq1.open(workspace + "/retained-pairs_1.fq");
    fq2.open(workspace + "/retained-pairs_2.fq");

    config = parse_config(workdir + "/config.txt");

    ctpl::thread_pool thread_pool(config.threads);

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());

    std::vector<std::future<void> > futures;
    for (int t = 0; t < bam_file->header->n_targets; t++) {
        std::future<void> future = thread_pool.push(filter_by_qname, bam_file->file->fn, bam_file->header->target_name[t]);
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

    fq1.close();
    fq2.close();
}