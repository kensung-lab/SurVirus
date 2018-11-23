#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include <iostream>
#include <atomic>
#include <htslib/sam.h>
#include "sam_utils.h"
#include "config.h"

typedef uint8_t sv_type_t;
static struct sv_types_t {
    sv_type_t DEL, INV, DUP, TRA, INS, NOV;
    size_t n_types = 6;

    sv_types_t() : DEL(1), INV(2), DUP(3), TRA(4), INS(5), NOV(6) {}
} SV_TYPES;

std::string svt_to_str(sv_type_t svt) {
    if (svt == SV_TYPES.DEL) return "DEL";
    else if (svt == SV_TYPES.INV) return "INV";
    else if (svt == SV_TYPES.DUP) return "DUP";
    else if (svt == SV_TYPES.TRA) return "TRA";
    else if (svt == SV_TYPES.INS) return "INS";
    else if (svt == SV_TYPES.NOV) return "NOV";
    else return "";
}

sv_type_t str_to_svt(std::string str) {
    if (str == "DEL") return SV_TYPES.DEL;
    else if (str == "INV") return SV_TYPES.INV;
    else if (str == "DUP") return SV_TYPES.DUP;
    else if (str == "TRA") return SV_TYPES.TRA;
    else if (str == "INS") return SV_TYPES.INS;
    else if (str == "NOV") return SV_TYPES.NOV;
    else return 0;
};

sv_type_t disct_to_svt(disc_type_t dt) {
    if (dt == DISC_TYPES.DC) return SV_TYPES.TRA;
    else if (dt == DISC_TYPES.LI) return SV_TYPES.DEL;
    else if (dt == DISC_TYPES.OW) return SV_TYPES.DUP;
    else if (dt == DISC_TYPES.SS) return SV_TYPES.INV;
    else if (dt == DISC_TYPES.SI) return SV_TYPES.INS;
    else return SV_TYPES.NOV;
}


struct anchor_t {
    static constexpr const char* pattern = "%c:%d:%d:%d:%d";

    char dir;
    int contig_id, start, end;
    int sc_reads;

    anchor_t() {}

    anchor_t(char dir, int contig_id, int start, int end, int sc_reads) : dir(dir), contig_id(contig_id), start(start),
                                                                          end(end), sc_reads(sc_reads) {}

    anchor_t(char* s) {
        sscanf(s, pattern, &dir, &contig_id, &start, &end, &sc_reads);
    }

    int pos() { return dir == 'L' ? start : end; }

    static bool check_clipped(anchor_t& clipped, anchor_t& other) {
        if (clipped.dir == 'L') {
            return clipped.pos() <= other.pos()+10;
        } else {
            return clipped.pos() >= other.pos()-10;
        }
    }

    static bool can_merge(anchor_t& a1, anchor_t& a2, int max_is) {
        if (a1.contig_id != a2.contig_id || a1.dir != a2.dir) return false;
        if (a1.sc_reads >= 2 && !check_clipped(a1, a2)) return false;
        if (a2.sc_reads >= 2 && !check_clipped(a2, a1)) return false;
        return std::max(a1.end,  a2.end)-std::min(a1.start, a2.start) <= max_is;
    }

    static anchor_t merge(anchor_t& a1, anchor_t& a2) {
        return anchor_t(a1.dir, a1.contig_id, std::min(a1.start, a2.start), std::max(a1.end,  a2.end), a1.sc_reads+a2.sc_reads);
    }

    static int distance(anchor_t& a1, anchor_t& a2) {
        int overlap = std::min(a1.end, a2.end) - std::max(a1.start, a2.start);
        if (overlap > 0) {
            int l2 = std::min(a1.size(), a2.size());
            return -100*overlap/l2;
        } else {
            return std::min(std::abs(a1.start-a2.end), std::abs(a2.start-a1.end));
        }
    }

    std::string to_str() {
        char buffer[100];
        sprintf(buffer, pattern, dir, contig_id, start, end, sc_reads);
        return std::string(buffer);
    }

    int size() {
        return end-start+1;
    }
};

bool operator < (const anchor_t& a1, const anchor_t& a2) {
    if (a1.start != a2.start) return a1.start < a2.start;
    return a1.end < a2.end;
}
struct cluster_t {
    static constexpr const char* scan_pattern = "%[^-]-%[^-]-%d";
    static constexpr const char* print_pattern = "%s-%s-%d";

    int id; // not necessarily set, use if needed
    anchor_t a1, a2;
    disc_type_t dt;
    int disc_pairs;
    bool dead = false;

    cluster_t(const anchor_t &a1, const anchor_t &a2, disc_type_t dt, int disc_pairs) : a1(a1), a2(a2), dt(dt),
                                                                                        disc_pairs(disc_pairs) {}

    cluster_t(std::string& s, disc_type_t dt) : dt(dt) {
        char anchor1[100], anchor2[100];
        sscanf(s.c_str(), scan_pattern, anchor1, anchor2, &disc_pairs);
        a1 = anchor_t(anchor1);
        a2 = anchor_t(anchor2);
    }

    static bool can_merge(cluster_t* c1, cluster_t* c2, int max_is) {
        return anchor_t::can_merge(c1->a1, c2->a1, max_is) && anchor_t::can_merge(c1->a2, c2->a2, max_is);
    }
    static int distance(cluster_t* c1, cluster_t* c2) {
        if (c1->a1.contig_id != c2->a1.contig_id || c1->a1.dir != c2->a1.dir ||
            c1->a2.contig_id != c2->a2.contig_id || c1->a2.dir != c2->a2.dir) {
            return INT32_MAX;
        }
        return anchor_t::distance(c1->a1, c2->a1) + anchor_t::distance(c1->a2, c2->a2);
    }

    static cluster_t* merge(cluster_t* c1, cluster_t* c2) {
        return new cluster_t(anchor_t::merge(c1->a1, c2->a1), anchor_t::merge(c1->a2, c2->a2), c1->dt,
                           c1->disc_pairs+c2->disc_pairs);
    }

    std::string to_str() {
        char buffer[1000];
        sprintf(buffer, print_pattern, a1.to_str().c_str(), a2.to_str().c_str(), disc_pairs);
        return std::string(buffer);
    }
};


struct cc_distance_t {
    int distance;
    cluster_t* c1,* c2;

    cc_distance_t(int distance, cluster_t *c1, cluster_t *c2) : distance(distance), c1(c1), c2(c2) {}
};
bool operator < (const cc_distance_t& ccd1, const cc_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance > ccd2.distance;
}


struct breakpoint_t {
    static constexpr const char* pattern = "%c:%d:%d:%d:%d:%d";

    char dir;
    int contig_id, start, end;
    int sc_reads;
    int spanning_reads = 0;

    breakpoint_t() {}

    breakpoint_t(anchor_t anchor) : dir(anchor.dir), contig_id(anchor.contig_id),
                                    start(anchor.start), end(anchor.end), sc_reads(anchor.sc_reads) {}

    breakpoint_t(char* s) {
        sscanf(s, pattern, &dir, &contig_id, &start, &end, &sc_reads, &spanning_reads);
    }

    int pos() { return dir == 'L' ? start : end; }

    std::string to_str() {
        char buffer[100];
        sprintf(buffer, pattern, dir, contig_id, start, end, sc_reads, spanning_reads);
        return std::string(buffer);
    }
};

std::atomic<int> pred_id(0);

struct prediction_t {
    static constexpr const char* scan_pattern = "%d,%[^,],%[^,],%[^,],%d,%d,%d,%lf,%lf";
    static constexpr const char* print_pattern = "%d,%s,%s,%s,%d,%d,%d,%lf,%lf";

    int id;
    breakpoint_t bp1, bp2;
    sv_type_t sv_type;
    int disc_pairs;
    int size = INT32_MAX, conf_ival = 0;
    double pval = -1.0, shift_pval = 1.0;

    prediction_t(cluster_t* c, disc_type_t dt) : id(pred_id++), bp1(c->a1), bp2(c->a2), sv_type(disct_to_svt(dt)),
                                                 disc_pairs(c->disc_pairs) {}

    prediction_t(std::string& line) {
        char breakpoint1[1000], breakpoint2[1000];
        char svt[10];
        sscanf(line.c_str(), scan_pattern, &id, breakpoint1, breakpoint2, svt, &disc_pairs, &size, &conf_ival, &pval, &shift_pval);
        if (id >= pred_id) pred_id = id+1;
        bp1 = breakpoint_t(breakpoint1);
        bp2 = breakpoint_t(breakpoint2);
        sv_type = str_to_svt(svt);
    }

    // len is simply the difference between the breakpoints
    int len() { return bp1.contig_id == bp2.contig_id ? bp2.pos()-bp1.pos() : INT32_MAX; }

    // size is: len if soft-clipped, the estimated size if not
    int get_size() { return bp1.sc_reads > 0 && sv_type != SV_TYPES.INS ? len() :
                            (size == INT32_MAX ? 0 : size); }

    std::string to_str() {
        char buffer[10000];
        sprintf(buffer, print_pattern, id, bp1.to_str().c_str(), bp2.to_str().c_str(), svt_to_str(sv_type).c_str(),
                disc_pairs, size, conf_ival, pval, shift_pval);
        return std::string(buffer);
    }
};

#endif //SURVEYOR_CLUSTER_H
