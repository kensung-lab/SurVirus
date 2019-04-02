#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include "config.h"


enum strand_t {
    F, R
};
struct anchor_t {
    strand_t strand;
    int contig_id, start, end;
    bool dead = false;

    anchor_t() {}

    anchor_t(strand_t strand, int contig_id, int start, int end) : strand(strand), contig_id(contig_id), start(start),
                                                                          end(end) {}

    static bool can_merge(anchor_t& a1, anchor_t& a2, int max_is) {
        if (a1.contig_id != a2.contig_id || a1.strand != a2.strand) return false;
        return std::max(a1.end,  a2.end)-std::min(a1.start, a2.start) <= max_is;
    }

    static anchor_t* merge(anchor_t* a1, anchor_t* a2) {
        return new anchor_t(a1->strand, a1->contig_id, std::min(a1->start, a2->start), std::max(a1->end,  a2->end));
    }

    static int distance(anchor_t& a1, anchor_t& a2) {
        int overlap = std::min(a1.end, a2.end) - std::max(a1.start, a2.start);
        if (overlap > 0) {
            int l2 = std::min(a1.len(), a2.len());
            return -100*overlap/l2;
        } else {
            return std::min(std::abs(a1.start-a2.end), std::abs(a2.start-a1.end));
        }
    }

    int len() {
        return end-start+1;
    }
};
bool operator < (const anchor_t& a1, const anchor_t& a2) {
    if (a1.start != a2.start) return a1.start < a2.start;
    return a1.end < a2.end;
}


struct aa_distance_t {
    int distance;
    anchor_t* a1, * a2;

    aa_distance_t(int distance, anchor_t *a1, anchor_t *a2) : distance(distance), a1(a1), a2(a2) {}
};
bool operator < (const aa_distance_t& aad1, const aa_distance_t& aad2) { // reverse op for priority queue
    return aad1.distance > aad2.distance;
}


#endif //SURVEYOR_CLUSTER_H
