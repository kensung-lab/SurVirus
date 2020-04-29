#include <iostream>
#include <fstream>
#include <random>
#include <htslib/sam.h>

#include "sam_utils.h"
#include "utils.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"

#include <sparsehash/dense_hash_map>

typedef unsigned long long ull;

int KMER_LEN = 13;
int KMER_BITS = KMER_LEN * 2;
ull KMER_MASK = (1ll << KMER_BITS)-1;


chr_seqs_map_t chr_seqs;
contig_map_t contig_map;
config_t config;

ull nucl_bm[256] = { 0 };
char bm_nucl[4] = { 'A', 'C', 'G', 'T' };

std::mutex mtx;
std::mutex* mtx_kmers, * mtx_regions;
std::vector<region_t*>* kmer_index, * kmer_index_rc;
std::vector<read_realignment_t>* reads_per_region, * reads_per_region_rc;
google::dense_hash_map<std::string, std::pair<int, const uint32_t*> > cstring_to_cigar;

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

void index_regions(int id, std::vector<region_t*>* regions, int start, int end) {
    mtx.lock();
    std::cerr << "Thread " << id << ": " << "indexing regions " << start << " to " << end << " " << std::endl;
    mtx.unlock();

    char region_str[100000];
    for (int r = start; r < end; r++) {
        region_t* region = (*regions)[r];
        char* seq = chr_seqs.get_seq(contig_map.id2name(region->contig_id));
        for (int i = 0; i < region->len(); i++) {
            region_str[i] = toupper(seq[region->start+i]);
        }
        region_str[region->len()] = '\0';
        index_seq(region_str, region->len(), region, false);
        get_rc(region_str, strlen(region_str));
        index_seq(region_str, region->len(), region, true);
    }
}

void write_qnames_indices(std::string& workspace, std::vector<bam1_t *>& reads) {
    std::ofstream qnames_map_out(workspace + "/qnames-map");
    for (int i = 0; i < reads.size(); i++) {
        qnames_map_out << i << " " << bam_get_qname(reads[i]) << " " << get_sequence(reads[i]) << std::endl;
    }
    qnames_map_out.close();
}

std::atomic<int> remappings(0);
inline read_realignment_t remap_read(std::string seq, bam1_t* read, region_t* region, bool is_rc, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), chr_seqs.get_seq(contig_map.id2name(region->contig_id)) + region->start,
                  region->end - region->start, filter, &alignment, 0);

    if (!accept_alignment(alignment, seq, config.min_sc_size)) return read_realignment_t();

    std::string cigar = alignment_cigar_to_bam_cigar(alignment.cigar);
    mtx.lock();
    if (!cstring_to_cigar.count(cigar)) {
        cstring_to_cigar[cigar] = cigar_str_to_array(cigar);
    }
    auto cigar_v = cstring_to_cigar[cigar];
    mtx.unlock();

    remappings++;
    if (remappings % 1000000 == 0) {
        std::cerr << remappings << " remappings done." << std::endl;
    }

    return read_realignment_t(read, alignment.ref_begin, alignment.ref_end, cigar_v.first, cigar_v.second, alignment.sw_score, is_rc);
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

void associate_reads_to_regions(int id, std::vector<bam1_t*>* host_reads, int start, int end) {
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

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    for (int r = start+1; r <= end; r++) {
        bam1_t* read = (*host_reads)[r-1];
        std::vector<std::pair<read_realignment_t, int> > realignments;

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

                        read_realignment_t rr = remap_read(seq, read, region, false, aligner, filter);
                        if (rr.accepted()) realignments.push_back({rr, region->id});
                    } else if (abs(last_read_for_region[region->id]) != r) {
                        last_read_for_region[region->id] = r;
                        last_pos[region->id] = i;
                    }
                }

                for (region_t* region : kmer_index_rc[kmer]) {
                    if (last_read_for_region_rc[region->id] == r && i - last_pos_rc[region->id] >= KMER_LEN) {
                        last_read_for_region_rc[region->id] = -r;

                        read_realignment_t rr = remap_read(seq_rc, read, region, true, aligner, filter);
                        if (rr.accepted()) realignments.push_back({rr, region->id});
                    } else if (abs(last_read_for_region_rc[region->id]) != r) {
                        last_read_for_region_rc[region->id] = r;
                        last_pos_rc[region->id] = i;
                    }
                }
            }
        }

        uint16_t max_score = 0;
        for (std::pair<read_realignment_t, int>& rr : realignments) {
            max_score = std::max(max_score, rr.first.score);
        }
        for (std::pair<read_realignment_t, int>& rr_and_region : realignments) {
            if (rr_and_region.first.score >= max_score*0.75) {
                mtx_regions[rr_and_region.second].lock();
                add_realignment_to_region(rr_and_region.first, rr_and_region.second, config.min_sc_size, reads_per_region, reads_per_region_rc);
                mtx_regions[rr_and_region.second].unlock();
            }
        }
    }

    delete[] last_read_for_region;
    delete[] last_pos;
    delete[] last_read_for_region_rc;
    delete[] last_pos_rc;
}

void write_alignment_to_bytes(char bytes[16], uint32_t region_id, uint32_t read_id, uint32_t cigar_id, uint16_t ref_begin,
                              bool is_rc, uint16_t score) {

    cigar_id = (is_rc ? 0x80000000 : 0) | (cigar_id & 0x7FFFFFFF); // forcing is_rc inside the cigar_id to make 16 bytes

    memcpy(bytes, &region_id, 4);
    memcpy(bytes+4, &read_id, 4);
    memcpy(bytes+8, &cigar_id, 4);
    memcpy(bytes+12, &ref_begin, 2);
    memcpy(bytes+14, &score, 2);
}


bool is_virus_region(const region_t *region) {
    return chr_seqs.is_virus(contig_map.id2name(region->contig_id));
};

int main(int argc, char* argv[]) {

    nucl_bm['A'] = nucl_bm['a'] = 0;
    nucl_bm['C'] = nucl_bm['c'] = 1;
    nucl_bm['G'] = nucl_bm['g'] = 2;
    nucl_bm['T'] = nucl_bm['t'] = 3;
    nucl_bm['N'] = nucl_bm['n'] = 0;

    std::string host_reference_fname  = argv[1];
    std::string virus_reference_fname  = argv[2];
    std::string workdir = argv[3];
    std::vector<std::string> bam_fnames, workspaces;
    int max_is = 0;
    for (int i = 4; i < argc; i++) {
        std::string workspace = argv[i];
        stats_t stats = parse_stats(workspace + "/stats.txt");
        max_is = std::max(max_is, stats.max_is);
    }

    config = parse_config(workdir + "/config.txt");

    contig_map = contig_map_t(workdir);

    open_samFile_t* host_and_virus_file = open_samFile((workdir + "/host-side.cs.bam").c_str(), false, false); // just for the header
    std::unordered_map<int, int> contig_tid2id;
    for (int i = 0; i < host_and_virus_file->header->n_targets; i++) {
        int contig_id = contig_map.name2id(host_and_virus_file->header->target_name[i]);
        contig_tid2id[i] = contig_id;
    }

    chr_seqs.read_fasta_into_map(host_reference_fname);
    chr_seqs.read_fasta_into_map(virus_reference_fname, true, true);

    /* == LOAD HOST READS == */
    std::vector<bam1_t*> host_reads;
    load_reads(workdir + "/host-side.cs.bam", host_reads, false);
    /* == == */


    /* == EXTRACT REGIONS == */
    std::vector<region_t*> host_regions;
    extract_regions(host_reads, host_regions, max_is, contig_map, contig_tid2id, chr_seqs);

    // filter virus regions
    for (int i = 0; i < host_regions.size(); i++) {
        if (is_virus_region(host_regions[i])) {
            delete host_regions[i];
            host_regions[i] = nullptr;
        }
    }
    host_regions.erase(std::remove(host_regions.begin(), host_regions.end(), nullptr), host_regions.end());

    std::sort(host_regions.begin(), host_regions.end(), [](const region_t* reg1, const region_t* reg2) {
        if (reg1->contig_id == reg2->contig_id) return reg1->start < reg2->start;
        return reg1->contig_id < reg2->contig_id;
    });
    for (int i = 0; i < host_regions.size(); i++) {
        host_regions[i]->id = i;
    }

    std::ofstream host_regions_file(workdir + "/host-regions");
    for (region_t* region : host_regions) {
        host_regions_file << region->id << " " << contig_map.id2name(region->contig_id) << " " <<
                          region->start << " " << region->end << std::endl;
    }
    host_regions_file.close();
    /* == == */


    /* === COMPUTE READS-REGIONS ALIGNMENTS === */
    reads_per_region = new std::vector<read_realignment_t>[host_regions.size()];
    reads_per_region_rc = new std::vector<read_realignment_t>[host_regions.size()];
    cstring_to_cigar.set_empty_key("");

    // INDEX REGIONS
    int n_kmers = 1 << KMER_BITS;
    kmer_index = new std::vector<region_t *>[n_kmers];
    kmer_index_rc = new std::vector<region_t *>[n_kmers];
    mtx_kmers = new std::mutex[n_kmers];

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
    delete[] mtx_kmers;

    // ASSOCIATE READS TO REGIONS
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(host_reads), std::end(host_reads), rng);

    write_qnames_indices(workdir, host_reads);

    mtx_regions = new std::mutex[host_regions.size()];

    ctpl::thread_pool thread_pool2(config.threads);
    futures.clear();

    int reads_chunks = config.threads * 5;
    int reads_per_thread = host_reads.size() / reads_chunks;

    for (int i = 0; i < reads_chunks-1; i++) {
        std::future<void> future = thread_pool2.push(associate_reads_to_regions, &host_reads,
                                                     i * reads_per_thread, (i + 1) * reads_per_thread);
        futures.push_back(std::move(future));
    }
    future = thread_pool2.push(associate_reads_to_regions, &host_reads,
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

    delete[] mtx_regions;
    delete[] kmer_index;
    delete[] kmer_index_rc;

    // write file
    std::cerr << "WRITING FILE reads.score.bin" << std::endl;
    std::string scores_file_path = workdir + "/reads.scores.bin";
    FILE* scores_file_out_bin = fopen(scores_file_path.c_str(), "wb");
    google::dense_hash_map<bam1_t*, uint32_t> read_ids;
    google::dense_hash_map<std::string, uint32_t> cigar_ids;
    read_ids.set_empty_key(NULL);
    cigar_ids.set_empty_key("");
    for (int i = 0; i < host_reads.size(); i++) {
        read_ids[host_reads[i]] = i;
    }
    for (int i = 0; i < host_regions.size(); i++) {
        auto& v_fwd = reads_per_region[i];
        for (read_realignment_t& rr : v_fwd) {
            char line[16]; // each line is 16 bytes
            std::string cigar = cigar_array_to_str(rr.cigar_len, rr.cigar);
            if (cigar_ids.count(cigar) == 0) cigar_ids[cigar] = cigar_ids.size();
            uint32_t cigar_id = cigar_ids[cigar];
            write_alignment_to_bytes(line, i, read_ids[rr.read], cigar_id, rr.offset_start, rr.rc, rr.score);
            fwrite(line, 1, 16, scores_file_out_bin);
        }

        auto& v_rev = reads_per_region_rc[i];
        for (read_realignment_t& rr : v_rev) {
            char line[16]; // each line is 16 bytes
            std::string cigar = cigar_array_to_str(rr.cigar_len, rr.cigar);
            if (cigar_ids.count(cigar) == 0) cigar_ids[cigar] = cigar_ids.size();
            uint32_t cigar_id = cigar_ids[cigar];
            write_alignment_to_bytes(line, i, read_ids[rr.read], cigar_id, rr.offset_start, rr.rc, rr.score);
            fwrite(line, 1, 16, scores_file_out_bin);
        }
    }

    std::ofstream cigar_map_out(workdir + "/cigars-map");
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
}