#include <iostream>
#include <unordered_set>
#include <mutex>
#include <bitset>
#include <unistd.h>
#include <zlib.h>
#include <htslib/kseq.h>
#include "libs/cptl_stl.h"

#define USE_BITSET


KSEQ_INIT(gzFile, gzread)

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

std::mutex mtx, mtx_out;
config_t config;
std::ofstream retained_fq1, retained_fq2;

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

char* kstring_to_cstr(kstring_t kstring) {
    char* cstr = (char*) malloc(kstring.l+1);
    strncpy(cstr, kstring.s, kstring.l);
    cstr[kstring.l] = '\0';
    return cstr;
}

struct read_t {
    char* name;
    char* seq;
    char* qual;

    read_t(kseq_t* kseq) {
        name = kstring_to_cstr(kseq->name);
        seq = kstring_to_cstr(kseq->seq);
        qual = kstring_to_cstr(kseq->qual);
    }

    void clear() {
        free(name);
        free(seq);
        free(qual);
    }
};
typedef std::pair<read_t, read_t> read_pair;
const int MAX_READS_TO_PROCESS = 10000;

bool is_virus_read(read_t read) {
    ull kmer = 0;
    int hit = 0, len = 0;

    const int seq1_len = strlen(read.seq);
    for (int i = 0; i < seq1_len; i++) {
        ull nv = nucl_bm[read.seq[i]];
        kmer = ((kmer << 2) | nv) & KMER_MASK;
        len++;

        if (len >= KMER_LEN) {
            print(kmer, KMER_LEN);
            for (int j = 0; j < NUMBER_OF_SEGS; j++) {
                ull seg = mask(kmer, j);
                if (check(seg, j)) {
                    hit++;
                    break;
                }
            }
        }

        if (hit >= 2) break;
    }

    return hit >= 2;
}

void isolate(int id, kseq_t* seq1, kseq_t* seq2) {

    std::vector<read_pair> read_pairs;
    do {
        read_pairs.clear();

        mtx.lock();
        for (int i = 0; i < MAX_READS_TO_PROCESS && kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0; i++) {
            read_t r1(seq1), r2(seq2);
            read_pairs.push_back({r1, r2});
        }
        mtx.unlock();

        std::vector<read_pair> to_write;
        for (read_pair& rp : read_pairs) {
            if (is_virus_read(rp.first) || is_virus_read(rp.second)) {
                to_write.push_back(rp);
            }
        }

        mtx_out.lock();
        for (read_pair& rp : to_write) {
            retained_fq1 << "@" << rp.first.name << std::endl;
            retained_fq1 << rp.first.seq << std::endl;
            retained_fq1 << "+" << rp.first.name << std::endl;
            retained_fq1 << rp.first.qual << std::endl;

            retained_fq2 << "@" << rp.second.name << std::endl;
            retained_fq2 << rp.second.seq << std::endl;
            retained_fq2 << "+" << rp.second.name << std::endl;
            retained_fq2 << rp.second.qual << std::endl;
        }
        mtx_out.unlock();

        for (read_pair& rp : read_pairs) {
            rp.first.clear();
            rp.second.clear();
        }
    } while (!read_pairs.empty());

}

int main(int argc, char* argv[]) {

    nucl_bm['A'] = 0;
    nucl_bm['C'] = 1;
    nucl_bm['G'] = 2;
    nucl_bm['T'] = 3;
    nucl_bm['N'] = 0;

    std::string fq1_fname = argv[1];
    std::string fq2_fname = argv[2];
    std::string host_ref = argv[3];
    std::string virus_ref = argv[4];
    std::string workdir = argv[5];
    std::string workspace = argv[6];

    gzFile fastaf = gzopen(virus_ref.c_str(), "r");
    kseq_t* seq = kseq_init(fastaf);
    while (kseq_read(seq) >= 0) {
        for (int i = 0; i < seq->seq.l; i++) {
            seq->seq.s[i] = toupper(seq->seq.s[i]);
        }
        index_seq(seq->seq.s, seq->seq.l);
        get_rc(seq->seq.s, seq->seq.l);
        index_seq(seq->seq.s, seq->seq.l);
    }
    kseq_destroy(seq);
    gzclose(fastaf);

    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';

    config = parse_config(workdir + "/config.txt");

    gzFile fq1f = gzopen(fq1_fname.c_str(), "r");
    gzFile fq2f = gzopen(fq2_fname.c_str(), "r");
    kseq_t* seq1 = kseq_init(fq1f);
    kseq_t* seq2 = kseq_init(fq2f);
    retained_fq1.open(workspace + "/retained-pairs_1.fq");
    retained_fq2.open(workspace + "/retained-pairs_2.fq");

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (int i = 0; i < config.threads; i++) {
        std::future<void> future = thread_pool.push(isolate, seq1, seq2);
        futures.push_back(std::move(future));
    }
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cout << s << std::endl;
        }
    }

    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fq1f);
    gzclose(fq2f);
    retained_fq1.close();
    retained_fq2.close();
}
