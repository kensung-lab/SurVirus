#ifndef SURVEYOR_SW_UTILS_H
#define SURVEYOR_SW_UTILS_H

#include "libs/ssw.h"
#include "libs/ssw_cpp.h"

int get_matches(StripedSmithWaterman::Alignment& alignment) {
    int m = 0;
    for (uint32_t c : alignment.cigar) {
        if (cigar_int_to_op(c) == 'M') {
            m += cigar_int_to_len(c);
        }
    }
    return m - alignment.mismatches;
}

int get_gaps(StripedSmithWaterman::Alignment& alignment) {
    int g = 0;
    for (uint32_t c : alignment.cigar) {
        if (cigar_int_to_op(c) == 'D' || cigar_int_to_op(c) == 'I') {
            g += cigar_int_to_len(c);
        }
    }
    return g - alignment.mismatches;
}

#endif //SURVEYOR_SW_UTILS_H
