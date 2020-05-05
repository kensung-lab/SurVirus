import pyfaidx, os, gzip, pysam
import numpy as np
from random_pos_generator import RandomPositionGenerator

MAX_READS = 1000
MAX_ACCEPTABLE_IS = 20000
GEN_DIST_SIZE = 10000

def get_sampling_regions_from_bam(reference, bam_file, wgs, covered_regions_bed_file):
    sampling_regions = []
    if covered_regions_bed_file:
        with open(covered_regions_bed_file) as sr_inf:
            for line in sr_inf:
                sl = line.split()
                sampling_regions += [(sl[0], int(sl[1]), int(sl[2]))]
        return sampling_regions

    if wgs:
        reference_fa = pyfaidx.Fasta(reference)
        for k in reference_fa.keys():
            sampling_regions += [(k, 1, len(reference_fa[k]))]
    else:
        curr_cluster = [None, 0, 0, 0]
        for r in bam_file.fetch(until_eof=True):
            if r.is_reverse or not r.is_proper_pair or r.is_unmapped or r.is_secondary \
                    or r.is_supplementary: continue

            if r.reference_name != curr_cluster[0] or r.reference_start > curr_cluster[2]:
                if curr_cluster[0] and curr_cluster[3] >= 10:
                    sampling_regions += [(curr_cluster[0], curr_cluster[1], curr_cluster[2])]
                curr_cluster = [r.reference_name, r.reference_start, r.reference_start+r.template_length, 1]
            else:
                curr_cluster[2] = max(curr_cluster[1], r.reference_start+r.template_length)
                curr_cluster[3] += 1

    return sampling_regions


def get_max_is_from_bam(reference_path, bam_files, wgs, covered_regions_bed_file):
    sampling_regions = get_sampling_regions_from_bam(reference_path, bam_files[0], wgs, covered_regions_bed_file)

    rand_pos_gen = RandomPositionGenerator(sampling_regions)
    random_positions = []
    for i in range(1, GEN_DIST_SIZE * 2):
        random_positions.append(rand_pos_gen.next())

    max_read_len = 0
    max_is_list = list()
    for file_index, bam_file in enumerate(bam_files):
        read_len = 0
        for i, read in enumerate(bam_file.fetch(until_eof=True)):
            if i > MAX_READS: break
            read_len = max(read_len, read.query_length)
        max_read_len = max(max_read_len, read_len)

        general_dist = []
        for e in random_positions:
            if len(general_dist) > GEN_DIST_SIZE: break

            chr, pos = e
            used = set()
            i = 0
            for read in bam_file.fetch(contig=chr, start=pos):
                if (read.reference_start, read.next_reference_start) in used: continue  # try not to use duplicates

                if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and 0 < read.template_length:
                    if i > 100: break
                    i += 1
                    general_dist.append(read.template_length)
                    used.add((read.reference_start, read.next_reference_start))

        mean_is = int(np.mean(general_dist))
        higher_stddev_is = int(np.sqrt(np.mean([(x - mean_is) ** 2 for x in general_dist if x > mean_is])))

        max_is = mean_is + 5 * higher_stddev_is
        max_is_list.append(max_is)

    return max_read_len, max_is_list


READS_TO_MAP = GEN_DIST_SIZE * 10

def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rb')
    else:
        return open(filename, 'r')

def get_max_is_from_fq(workdir, fq1, fq2, reference, bwa_exec, threads):
    with open_by_suffix(fq1) as fq1_f, open_by_suffix(fq2) as fq2_f, \
        open("%s/head_1.fq" % workdir, "w") as head_fq1, open("%s/head_2.fq" % workdir, "w") as head_fq2:
        for i in xrange(READS_TO_MAP * 4):
            line1, line2 = next(fq1_f, None), next(fq2_f, None)
            if not line1 or not line2: break
            head_fq1.write(line1)
            head_fq2.write(line2)

    bwa_cmd = "%s mem -t %d %s %s/head_1.fq %s/head_2.fq > %s/head.sam" \
              % (bwa_exec, threads, reference, workdir, workdir, workdir)
    os.system(bwa_cmd)

    max_read_len = 0
    with pysam.AlignmentFile("%s/head.sam" % workdir) as head_f:
        general_dist = []
        for read in head_f.fetch(until_eof=True):
            if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and 0 < read.template_length:
                general_dist.append(read.template_length)
                max_read_len = max(max_read_len, read.query_length)

            if len(general_dist) > GEN_DIST_SIZE: break

        mean_is = int(np.mean(general_dist))
        higher_stddev_is = int(np.sqrt(np.mean([(x - mean_is) ** 2 for x in general_dist if x > mean_is])))

        max_is = mean_is + 5 * higher_stddev_is
    return max_read_len, max_is