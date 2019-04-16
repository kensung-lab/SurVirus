#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <htslib/sam.h>

char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
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
    char seq[10000];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (original_seq && bam_is_rev(r)) get_rc(seq, r->core.l_qseq);
    return std::string(seq);
}


int main(int argc, char* argv[]) {

    std::string workspace = argv[1];
    std::string target_qname = argv[2];

    std::string scores_file_path = workspace + "/reads.scores.bin";
    FILE* scores_file_in_bin = fopen(scores_file_path.c_str(), "rb");
    if (scores_file_in_bin) {
        std::set<int> target_read_ids;
        std::ifstream qnames_map_in(workspace + "/qnames-map");
        int i; std::string qname, seq;
        while (qnames_map_in >> i >> qname >> seq) {
            if (qname == target_qname) {
                target_read_ids.insert(i);
            }
        }
        qnames_map_in.close();

        std::unordered_map<uint32_t, std::string> id_to_cigar;
        std::ifstream cigars_map_in(workspace + "/cigars-map");
        std::string cigar_str;
        while (cigars_map_in >> i >> cigar_str) {
            id_to_cigar[i] = cigar_str;
        }

        std::unordered_map<uint32_t, std::string> id_to_region;
        std::ifstream regions_map_in(workspace + "/host-regions");
        std::string chr, start, end;
        while (regions_map_in >> i >> chr >> start >> end) {
            id_to_region[i] = chr + ":" + start + "-" + end;
        }

        char line[16];
        while (fread(line, 16, 1, scores_file_in_bin)) {
            uint32_t region_id = *((uint32_t *) (line + 0));
            uint32_t read_id = *((uint32_t *) (line + 4));
            if (!target_read_ids.count(read_id)) continue;

            uint32_t cigar_id = *((uint32_t *) (line + 8));
            bool is_rc = cigar_id & 0x80000000;
            std::string cigar_str = id_to_cigar[cigar_id & 0x7FFFFFFF];
            uint16_t offset_start = *((uint16_t *) (line + 12));
            uint16_t score = *((uint16_t *) (line + 14));
            std::cout << id_to_region[region_id] << " " << offset_start << " " << " " << cigar_str << " " <<
            score << " " << is_rc << std::endl;
        }
        fclose(scores_file_in_bin);
    }

}