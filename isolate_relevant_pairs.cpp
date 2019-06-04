#include <iostream>
#include <fstream>
#include <unordered_set>
#include <mutex>
#include <bitset>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include "libs/cptl_stl.h"

#define USE_BITSET


KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "config.h"

typedef unsigned long long ull;

const int KMER_LEN = 18;
const int NUMBER_OF_SEGS = 6;

const int KMER_BITS = KMER_LEN * 2;

const int SEG_LEN = KMER_LEN/NUMBER_OF_SEGS;
const int SEG_BITS = SEG_LEN * 2;

const ull KMER_MASK = (1ll << KMER_BITS)-1;

const int MASKED_KMER_LEN = KMER_LEN - SEG_LEN;
const int MASKED_KMER_BITS = MASKED_KMER_LEN * 2;
const ull MASKED_KMER_MASK = (1ll << MASKED_KMER_BITS)-1;

std::atomic<int> hits(0);

ull nucl_bm[256] = { 0 };
char bm_nucl[4] = { 'A', 'C', 'G', 'T' };
char nucl2chr[16];

#ifdef USE_BITSET
std::bitset<(1ll << MASKED_KMER_BITS)> segs[NUMBER_OF_SEGS];
#else
std::unordered_set<ull>  segs[NUMBER_OF_SEGS];
#endif

inline bool check(ull masked_kmer, int seg_n) {
#ifdef USE_BITSET
    return segs[seg_n].test(masked_kmer);
#else
    return segs[seg_n].count(masked_kmer);
#endif
}

inline void insert(ull masked_kmer, int seg_n) {
#ifdef USE_BITSET
    segs[seg_n].set(masked_kmer);
#else
    segs[seg_n].insert(masked_kmer);
#endif
}

std::mutex mtx;
config_t config;
std::ofstream fout;

std::unordered_set<std::string> host_names;


inline ull mask(ull kmer, int seg_n) {
    ull first_n_segs_mask = (1ll << SEG_BITS*seg_n)-1;
    ull first_n_segs = kmer & first_n_segs_mask;
    ull remaining_segs_shifted = (kmer >> SEG_BITS) & ~first_n_segs_mask;
    return (first_n_segs | remaining_segs_shifted) & MASKED_KMER_MASK;
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

void index_seq(char* seq, size_t len) {
    ull kmer = 0;
    for (int i = 0; i < len; i++) {
        ull nv = nucl_bm[seq[i]];
        kmer = ((kmer << 2) | nv) & KMER_MASK;

        if (i >= KMER_LEN-1) {
            for (int j = 0; j < NUMBER_OF_SEGS; j++) {
                ull seg = mask(kmer, j);

                if (valid_kmer(seg, MASKED_KMER_LEN)) {
                    insert(seg, j);
                }
            }
        }
    }
}

void isolate(int id, char* bam_fname, int t_num) {

    open_samFile_t* bam_file = open_samFile(bam_fname, false);

    std::string t_name = bam_file->header->target_name[t_num];
    bool is_host = host_names.count(t_name);

    char region[1000];
    sprintf(region, "%s:%d-%d", t_name.c_str(), 1, bam_file->header->target_len[t_num]);

    mtx.lock();
    std::cerr << "Detecting relevant pairs for " << bam_file->header->target_name[t_num] << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region);

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        uint32_t* c = bam_get_cigar(read);
        if (is_host && bam_cigar_op(c[0]) == BAM_CMATCH && read->core.n_cigar == 1) continue; // ignore if all matches

        if (read->core.flag & BAM_FDUP) continue;

        ull kmer = 0;
        const uint8_t* bam_seq = bam_get_seq(read);
        int hit = 0;
        int len = 0;
        for (int i = 0; i < read->core.l_qseq; i++) {
            mtx.lock();
            ull nv = nucl_bm[nucl2chr[bam_seqi(bam_seq, i)]];
            mtx.unlock();
            kmer = ((kmer << 2) | nv) & KMER_MASK;
            len++;
            if (len >= KMER_LEN) {
                for (int j = 0; j < NUMBER_OF_SEGS; j++) {
                    ull seg = mask(kmer, j);
                    if (check(seg, j)) {
                        hit++;
                        len -= 0;
                        break;
                    }
                }
            }

            if (hit >= 2) break;
        }

        if (hit >= 2) {
            hits++;
            mtx.lock();
            fout << bam_get_qname(read) << "\n";
            mtx.unlock();
        }
    }
    bam_destroy1(read);

    sam_itr_destroy(iter);

    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {

    nucl_bm['A'] = 0;
    nucl_bm['C'] = 1;
    nucl_bm['G'] = 2;
    nucl_bm['T'] = 3;
    nucl_bm['N'] = 0;

    std::string bam_fname = argv[1];
    std::string host_ref = argv[2];
    std::string virus_ref = argv[3];
    std::string workdir = argv[4];
    std::string workspace = argv[5];

    FILE* fastaf = fopen(host_ref.c_str(), "r");
    kseq_t* seq = kseq_init(fileno(fastaf));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        std::string qname = seq->name.s;
        host_names.insert(qname);
    }
    kseq_destroy(seq);
    fclose(fastaf);

    fastaf = fopen(virus_ref.c_str(), "r");
    seq = kseq_init(fileno(fastaf));
    while ((l = kseq_read(seq)) >= 0) {
        index_seq(seq->seq.s, seq->seq.l);
        get_rc(seq->seq.s, seq->seq.l);
        index_seq(seq->seq.s, seq->seq.l);
    }
    kseq_destroy(seq);
    fclose(fastaf);

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());

    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';

    config = parse_config(workdir + "/config.txt");

    ctpl::thread_pool thread_pool(config.threads);

    fout.open(workspace + "/qnames-to-keep");

    std::vector<std::future<void> > futures;
    for (int t = 0; t < bam_file->header->n_targets; t++) {
        std::future<void> future = thread_pool.push(isolate, bam_file->file->fn, t);
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

    fout.close();
}
