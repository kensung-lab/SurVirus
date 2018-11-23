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
std::vector<std::ofstream*> mate_seqs_writers_by_tid; // indexed by tid, not contig_id
std::unordered_set<std::string> virus_names;

std::mutex mtx;

std::ofstream virus_clip_out, reference_clip_out;

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

void categorize(int id, std::string contig, std::string bam_fname, int target_len, std::vector<bam1_t*>* virus_side_reads,
                std::vector<bam1_t*>* reference_side_reads, samFile* anchor_writer, std::ofstream* clip_writer) {

    open_samFile_t* bam_file_open = open_samFile(bam_fname.data());

    samFile* bam_file = bam_file_open->file;
//    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
//    if (bam_file == NULL) {
//        throw "Unable to open BAM file.";
//    }
//
    hts_idx_t* idx = bam_file_open->idx;
//    hts_idx_t* idx = sam_index_load(bam_file, bam_fname.c_str());
//    if (idx == NULL) {
//        throw "Unable to open BAM index.";
//    }
//
    bam_hdr_t* header = bam_file_open->header;
//    bam_hdr_t* header = sam_hdr_read(bam_file);
//    if (header == NULL) {
//        throw "Unable to open BAM header.";
//    }

    int contig_id = contig_name2id[contig];

    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), 1, target_len);

    mtx.lock();
    std::cout << "Categorizing " << region << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(idx, header, region);
    bam1_t* read = bam_init1();

    std::unordered_map<disc_type_t, samFile*> writers;

    std::deque<bam1_t*> two_way_buffer, forward_buffer;

    int i = 0;
    while (sam_itr_next(bam_file, iter, read) >= 0 && i < MAX_BUFFER_SIZE-1) {
        if (!is_unmapped(read)) {
            bam1_t *read2 = bam_dup1(read);
            add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
            add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
            i++;
        }
    }

    while (!forward_buffer.empty()) {
        while (sam_itr_next(bam_file, iter, read) >= 0) {
            if (!is_unmapped(read)) {
                bam1_t *read2 = bam_dup1(read);
                bam1_t* to_destroy = add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
                bam_destroy1(to_destroy);
                add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
                break;
            }
        }

        bam1_t* read = forward_buffer.front();
        forward_buffer.pop_front();

//        if (is_clipped(read, true) && !is_primary(read) && !virus_names.count(header->target_name[read->core.tid])) {
//            std::string target_name = header->target_name[read->core.tid];
//            std::string sa_target_name = bam_aux2Z(bam_aux_get(read, "SA"));
//            sa_target_name = sa_target_name.substr(0, sa_target_name.find(','));
//            if (!virus_names.count(target_name) && virus_names.count(sa_target_name)) {
//                std::cout << bam_get_qname(read) << " " << target_name << " " << sa_target_name << " "
//                << get_sequence(read) << std::endl;
//            }
//            continue;
//        }

//        bam_aux_get(read, "MQ");
//        int64_t mq = get_mq(read);
//        if (read->core.qual < MIN_MAPQ && mq < MIN_MAPQ) continue;

        // clipped read
        if (is_clipped(read) && get_clip_len(read) >= config.min_sv_len) {
            // && check_SNP(read, two_way_buffer, config.avg_depth)) {
            mtx.lock();
            int ok = sam_write1(anchor_writer, header, read);
            std::string seq = get_sequence(read);
            *clip_writer << ">" << bam_get_qname(read) << "\n";
            if (is_left_clipped(read)) {
                *clip_writer << seq.substr(0, get_left_clip_len(read)) << "\n";
            } else if (is_right_clipped(read)) {
                *clip_writer << seq.substr(seq.length()- get_right_clip_len(read)-1, get_right_clip_len(read)) << "\n";
            }
            mtx.unlock();
            if (ok < 0) throw "Failed to write to " + std::string(anchor_writer->fn);
        }

//        // we accept one of the mates having mapq 0 only if they are on different chromosomes
//        if ((read->core.qual < MIN_MAPQ || mq < MIN_MAPQ)
//            && !is_dc_pair(read) && !is_mate_unmapped(read)) {
//            continue;
//        }

        if (is_dc_pair(read) && get_avg_qual(read) > 10 && !is_poly_ACGT(read, true)) {

            std::string target_name = header->target_name[read->core.tid];
            std::string mate_target_name = header->target_name[read->core.mtid];
            mtx.lock();
            if (virus_names.count(target_name) > 0 && virus_names.count(mate_target_name) == 0) { // virus side
                virus_side_reads->push_back(bam_dup1(read));
//                int ok = sam_write1(virus_side_writer, header, read);
//                if (ok < 0) throw "Failed to write to " + std::string(virus_side_writer->fn);
            } else if (virus_names.count(target_name) == 0 && virus_names.count(mate_target_name) > 0) { // reference side
                reference_side_reads->push_back(bam_dup1(read));
//                int ok = sam_write1(reference_side_writer, header, read);
//                if (ok < 0) throw "Failed to write to " + std::string(reference_side_writer->fn);
            }
            mtx.unlock();


//                if (read->core.qual >= mq && check_SNP(read, two_way_buffer, config.avg_depth)) { // stable end
//                    if ((bam_is_rev(read) && !is_right_clipped(read)) || (!bam_is_rev(read) && !is_left_clipped(read))) {
//                        int ok = sam_write1(bam_is_rev(read) ? ldc_writer : rdc_writer, header, read);
//                        if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
//                    }
//                }
//                if (read->core.qual <= mq) { // save read seq for remapping
//                    const uint8_t* read_seq = bam_get_seq(read);
//                    char read_seq_chr[MAX_READ_SUPPORTED];
//                    for (int i = 0; i < read->core.l_qseq; i++) {
//                        read_seq_chr[i] = get_base(read_seq, i);
//                    }
//                    read_seq_chr[read->core.l_qseq] = '\0';
//                    mtx.lock();
//                    if (bam_is_rev(read)) {
//                        get_rc(read_seq_chr);
//                    }
//                    std::string qname = bam_get_qname(read);
//                    if (is_samechr(read)) {
//                        if (read->core.isize > 0) qname += "_1";
//                        else qname += "_2";
//                    }
//                    *mate_seqs_writers_by_tid[read->core.mtid] << qname << " ";
//                    *mate_seqs_writers_by_tid[read->core.mtid] << read_seq_chr << "\n";
//                    mtx.unlock();
//                }
//            }
        }
    }

    close_samFile(bam_file_open);
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

    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file == NULL) {
        std::cerr << "Unable to open BAM file." << std::endl;
        return -1;
    }

    int code = sam_index_build(bam_fname.c_str(), 0);
    if (code != 0) {
        throw "Cannot index " + std::string(bam_fname.c_str());
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    tid_to_contig_id.resize(header->n_targets);
//    for (int i = 0; i < header->n_targets; i++) {
//        int contig_id = contig_name2id[std::string(header->target_name[i])];
//        std::string fname = std::to_string(contig_id) + "-MATESEQS";
//        mate_seqs_writers_by_tid.push_back(new std::ofstream(workdir + "/workspace/" + fname));
//    }

    std::vector<bam1_t*> virus_side_reads, reference_side_reads;
    samFile* virus_anchor_writer = open_bam_writer(workspace, "virus-anchors.bam", header);
    samFile* reference_anchor_writer = open_bam_writer(workspace, "reference-anchors.bam", header);
    virus_clip_out.open(workspace + "/virus-clips.fa");
    reference_clip_out.open(workspace + "/reference-clips.fa");

    std::vector<std::future<void> > futures;
    for (int i = 0; i < header->n_targets; i++) {
        bool is_virus = virus_names.count(header->target_name[i]);
        std::future<void> future = thread_pool.push(categorize, header->target_name[i], bam_fname, header->target_len[i],
                                                    &virus_side_reads, &reference_side_reads,
                                                    is_virus ? virus_anchor_writer : reference_anchor_writer,
                                                    is_virus ? &virus_clip_out : &reference_clip_out);
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

    // get shared qnames
    std::unordered_set<std::string> virus_qnames;
    for (bam1_t* r : virus_side_reads) {
        virus_qnames.insert(bam_get_qname(r));
    }
    std::unordered_set<std::string> shared_qnames;
    for (bam1_t* r : reference_side_reads) {
        if (virus_qnames.count(bam_get_qname(r))) {
            shared_qnames.insert(bam_get_qname(r));
        }
    }

    // extract and sort shared reads
    std::vector<bam1_t*> final_virus_reads, final_reference_reads;
    for (bam1_t* r : virus_side_reads) {
        if (shared_qnames.count(bam_get_qname(r))) {
            final_virus_reads.push_back(r);
        }
    }
    std::sort(final_virus_reads.begin(), final_virus_reads.end(), [](const bam1_t* r1, const bam1_t* r2) {
        return strcmp(bam_get_qname(r1), bam_get_qname(r2)) < 0;});
    for (bam1_t* r : reference_side_reads) {
        if (shared_qnames.count(bam_get_qname(r))) {
            final_reference_reads.push_back(r);
        }
    }
    std::sort(final_reference_reads.begin(), final_reference_reads.end(), [](const bam1_t* r1, const bam1_t* r2) {
        return strcmp(bam_get_qname(r1), bam_get_qname(r2)) < 0;});

    std::ofstream virus_fq(workspace + "/virus-side.fq"), reference_fq(workspace + "/reference-side.fq");
    for (bam1_t* r : final_reference_reads) {
        reference_fq << print_fq(r);
    }
    for (bam1_t* r : final_virus_reads) {
        virus_fq << print_fq(r);
    }

    sam_close(virus_anchor_writer);
    sam_close(reference_anchor_writer);
    reference_clip_out.close();
    virus_clip_out.close();

//    for (int i = 0; i < header->n_targets; i++) {
//        mate_seqs_writers_by_tid[i]->close();
//        delete mate_seqs_writers_by_tid[i];
//    }
}
