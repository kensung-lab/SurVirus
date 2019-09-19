#ifndef SURVEYOR_SAM_UTILS_H
#define SURVEYOR_SAM_UTILS_H

#include <iostream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstring>
#include <htslib/sam.h>

extern int MIN_CLIP_LEN;
extern int MAX_READ_SUPPORTED;

int get_mate_endpos(bam1_t* r);

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}


bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}

bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid;
}

bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000;
}

bool is_valid(bam1_t* r) {
    return !is_unmapped(r) && !is_mate_unmapped(r) && is_primary(r);
}

int get_left_clip_len(bam1_t* r, bool include_HC = false) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' || (include_HC && bam_cigar_opchr(cigar[0]) == 'H') ?
            bam_cigar_oplen(cigar[0]) : 0;
}

bool is_left_clipped(bam1_t* r, bool include_HC = false) {
    return get_left_clip_len(r, include_HC) >= MIN_CLIP_LEN;
}

int get_right_clip_len(bam1_t *r, bool include_HC = false) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' || (include_HC && bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'H') ?
            bam_cigar_oplen(cigar[r->core.n_cigar-1]) : 0;
}

bool is_right_clipped(bam1_t* r, bool include_HC = false) {
    return get_right_clip_len(r, include_HC) >= MIN_CLIP_LEN;
}

std::string get_cigar_code(const uint32_t* cigar, int n_cigar) {
    std::stringstream ss;
    for (int i = 0; i < n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}
std::string get_cigar_code(bam1_t* r) {
    return get_cigar_code(bam_get_cigar(r), r->core.n_cigar);
}

int get_mate_endpos(bam1_t* r) {
    uint8_t *mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0, pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}


char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}

int get_avg_qual(bam1_t* r, int start, int end) {
    if (start >= end) return 0;
    int q = 0;
    for (int i = start; i < end; i++) {
        q += bam_get_qual(r)[i];
    }
    return q/(end-start);
}
int get_avg_qual(bam1_t* r, bool aligned_only = true) {
    int start = aligned_only ? get_left_clip_len(r) : 0;
    int end = r->core.l_qseq - (aligned_only ? get_right_clip_len(r) : 0);
    return get_avg_qual(r, start, end);
}

bool is_poly_ACGT(bam1_t* r, bool aligned_only = false) {
    int a = 0, c = 0, g = 0, t = 0;
    int start = aligned_only ? get_left_clip_len(r) : 0;
    int end = r->core.l_qseq - (aligned_only ? get_right_clip_len(r) : 0);
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = start; i < end; i++) {
        char base = get_base(bam_seq, i);
        if (base == 'A') a++;
        else if (base == 'C') c++;
        else if (base == 'G') g++;
        else if (base == 'T') t++;
    }

    int maxc = std::max(std::max(a,c), std::max(g,t));
    return double(maxc)/(end-start) >= 0.8;
}


void get_rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

void get_rc(char *read, int len) {
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

std::string get_sequence(bam1_t* r, bool original_seq = false) { // if original_seq == true, return the sequence found in fastx file
    char seq[MAX_READ_SUPPORTED];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (original_seq && bam_is_rev(r)) get_rc(seq, r->core.l_qseq);
    return std::string(seq);
}

bool is_low_complexity(bam1_t* read, bool rc, int start, int end) {

    int count[256];
    count[0] = 0;
    for (uint8_t i = 1; i < 16; i*=2) {
        for (uint8_t j = 1; j < 16; j *= 2) {
            uint8_t c = (i << 4) | j;
            count[c] = 0;
        }
    }

    std::string seq = get_sequence(read, true);
    if (rc) get_rc(seq);

    char chr2nucl[256];
    chr2nucl['A'] = 1; chr2nucl['C'] = 2; chr2nucl['G'] = 4; chr2nucl['T'] = 8; chr2nucl['N'] = 15;

    for (uint8_t i = start+1; i < end; i++) {
        uint8_t twobases = (chr2nucl[seq[i-1]] << 4) | chr2nucl[seq[i]];
        count[int(twobases)]++;
    }

    uint8_t top1 = 0, top2 = 0;
    for (uint8_t i = 1; i < 16; i*=2) {
        for (uint8_t j = 1; j < 16; j*=2) {
            uint8_t c = (i << 4) | j;
            if (count[c] > count[top1]) {
                top2 = top1;
                top1 = c;
            } else if (count[c] > count[top2]) {
                top2 = c;
            }
        }
    }

    bool is_lc = count[top1] + count[top2] >= (end-start+1)*0.75;
    return is_lc;
}

std::string get_left_clip(bam1_t* read) {
    std::string seq = get_sequence(read);
    return seq.substr(0, get_left_clip_len(read));
}
std::string get_right_clip(bam1_t* read) {
    std::string seq = get_sequence(read);
    return seq.substr(seq.length()- get_right_clip_len(read)-1, get_right_clip_len(read));
}

std::pair<int, const uint32_t*> cigar_str_to_array(std::string& cigar_str) {
    std::vector<uint32_t> opv;

    int pos = 0, prev = 0;
    std::string bam_ops = BAM_CIGAR_STR;
    while ((pos = cigar_str.find_first_of(bam_ops, prev)) != std::string::npos) {
        opv.push_back(bam_cigar_gen(std::stoi(cigar_str.substr(prev, pos-prev)), bam_ops.find(cigar_str[pos])));
        prev = pos+1;
    }

    uint32_t* opa = new uint32_t[opv.size()];
    std::copy(opv.begin(), opv.end(), opa);
    return std::make_pair(opv.size(), opa);
}
std::string cigar_array_to_str(int cigar_len, const uint32_t* cigar) {
    std::stringstream ss;
    for (int i = 0; i < cigar_len; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}

struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() {}

    open_samFile_t(samFile* file, bam_hdr_t* header, hts_idx_t* idx) : file(file), header(header), idx(idx) {}
};

open_samFile_t* open_samFile(const char* fname, bool index_file = false) {
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
    if (sam_file->idx == NULL) {
        throw "Unable to open index for " + std::string(fname);
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
    hts_idx_destroy(f->idx);
    bam_hdr_destroy(f->header);
    sam_close(f->file);
    delete f;
}

samFile* open_bam_writer(std::string workspace, std::string name, bam_hdr_t* header, bool sam = false) {
    std::string filename = workspace + "/" + name;
    samFile* writer = sam_open(filename.c_str(), sam ? "w" : "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + filename;
    }
    return writer;
}

template<typename T>
inline T max(T a, T b, T c) { return std::max(std::max(a,b), c); }

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

template<typename T>
inline T max(T a, T b, T c, T d, T e) { return max(std::max(a,b), std::max(c,d), e); }

template<typename T>
inline T min(T a, T b, T c) { return std::min(std::min(a,b), c); }

template<typename T>
inline T min(T a, T b, T c, T d) { return std::min(std::min(a,b), std::min(c,d)); }




#endif //SURVEYOR_SAM_UTILS_H
