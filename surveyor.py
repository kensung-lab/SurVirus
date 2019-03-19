import argparse, os
import pysam, pyfaidx
from random_pos_generator import RandomPositionGenerator
import numpy as np
from shutil import copyfile

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

cmd_parser = argparse.ArgumentParser(description='SurVirus, a virus integration caller.')
cmd_parser.add_argument('bam_files', help='Input bam files, separated by a semi-colon.')
cmd_parser.add_argument('workdir', help='Working directory for SurVirus to use.')
cmd_parser.add_argument('host_reference', help='Reference of the host organism in FASTA format.')
cmd_parser.add_argument('virus_reference', help='References of a list of viruses in FASTA format.')
cmd_parser.add_argument('host_and_virus_reference', help='Joint references of host and viruses.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--samtools', help='Samtools path.', default='samtools')
cmd_parser.add_argument('--bwa', default='bwa', help='BWA path.')
cmd_parser.add_argument('--minSVLen', type=int, default=18, help='Min SV length.')
cmd_parser.add_argument('--maxSCDist', type=int, default=10, help='Max SC distance.')
cmd_args = cmd_parser.parse_args()


# Create config file in workdir

config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("min_sv_len %d\n" % cmd_args.minSVLen)
config_file.write("max_sc_dist %d\n" % cmd_args.maxSCDist)
config_file.close();


# Find read length

bam_names = cmd_args.bam_files.split(';')
bam_files = [pysam.AlignmentFile(bam_name) for bam_name in bam_names]


# Generate general distribution of insert sizes
contig_map = open("%s/contig_map" % cmd_args.workdir, "w")
num_contigs = 0

reference_fa = pyfaidx.Fasta(cmd_args.host_and_virus_reference)
i = 1
for k in reference_fa.keys():
    contig_map.write("%s %d\n" % (k, i));
    i += 1
    num_contigs += 1
contig_map.close();


# count viruses
reference_fa = pyfaidx.Fasta(cmd_args.virus_reference)
n_viruses = len(reference_fa.keys())


reference_fa = pyfaidx.Fasta(cmd_args.host_reference)

rand_pos_gen = RandomPositionGenerator(reference_fa)
random_positions = []
for i in range(1,GEN_DIST_SIZE*2):
    if i % 100000 == 0: print i, "random positions generated."
    random_positions.append(rand_pos_gen.next())

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)


for file_index, bam_file in enumerate(bam_files):
    read_len = 0
    for i, read in enumerate(bam_file.fetch(until_eof=True)):
        if i > MAX_READS: break
        read_len = max(read_len, read.query_length)

    general_dist = []
    avg_depth = 0
    samplings = 0
    rnd_i = 0
    while len(general_dist) < GEN_DIST_SIZE:
        chr, pos = random_positions[rnd_i]
        rnd_i += 1

        if pos > len(reference_fa[chr])-10000:
            continue

        samplings += 1
        i = 0
        for read in bam_file.fetch(contig=chr, start=pos, end=pos+10000):
            if read.reference_start < pos:
                avg_depth += 1
            if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
            0 < read.template_length < MAX_ACCEPTABLE_IS and 'S' not in read.cigarstring and 'S' not in read.get_tag('MC'):
                if i > 100: break
                i += 1
                general_dist.append(read.template_length)

    mean_is = np.mean(general_dist)
    stddev_is = np.std(general_dist)
    avg_depth = float(avg_depth)/samplings;

    print "Average depth:", avg_depth

    general_dist = [x for x in general_dist if abs(x-mean_is) < 8*stddev_is]

    mean_is = int(np.mean(general_dist))
    lower_stddev_is = int(np.sqrt(np.mean([(mean_is-x)**2 for x in general_dist if x < mean_is])))
    higher_stddev_is = int(np.sqrt(np.mean([(x-mean_is)**2 for x in general_dist if x > mean_is])))

    min_is, max_is = mean_is-5*lower_stddev_is, mean_is+5*higher_stddev_is
    print mean_is, min_is, max_is

    folder = "%s/bam_%d/" % (cmd_args.workdir, file_index)
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open("%s/stats.txt" % folder, "w") as stat_file:
        stat_file.write("filename %s\n" % bam_file.filename)
        stat_file.write("min_is %d\n" % min_is)
        stat_file.write("avg_is %d\n" % mean_is)
        stat_file.write("max_is %d\n" % max_is)
        stat_file.write("read_len %d\n" % read_len)

for file_index, bam_file in enumerate(bam_files):
    bam_workspace = "%s/bam_%d/" % (cmd_args.workdir, file_index)

    isolate_cmd = "./isolate_relevant_pairs %s %s %s %s" % (bam_file.filename, cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    print "Executing:", isolate_cmd
    os.system(isolate_cmd)

    filter_by_qname_cmd = "./filter_by_qname %s %s %s" % (bam_file.filename, cmd_args.workdir, bam_workspace)
    print "Executing:", filter_by_qname_cmd
    os.system(filter_by_qname_cmd)

    samtools_sortbyqname_cmd = "%s sort -n -@ %d -o %s/retained-pairs.qname-sorted.bam %s/retained-pairs.bam" % \
                               (cmd_args.samtools, cmd_args.threads, bam_workspace, bam_workspace)
    print "Executing:", samtools_sortbyqname_cmd
    os.system(samtools_sortbyqname_cmd)

    samtools_dump_cmd = "%s fastq -1 %s/retained-pairs_1.fq -2 %s/retained-pairs_2.fq %s/retained-pairs.qname-sorted.bam" % \
                        (cmd_args.samtools, bam_workspace, bam_workspace, bam_workspace)
    print "Executing:", samtools_dump_cmd
    os.system(samtools_dump_cmd)

    bwa_cmd = "%s mem -t %d %s %s/retained-pairs_1.fq %s/retained-pairs_2.fq | %s view -b > %s/retained-pairs-remapped.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, \
                 bam_workspace, bam_workspace, cmd_args.samtools, bam_workspace)
    print "Executing:", bwa_cmd
    os.system(bwa_cmd)

    samtools_sort_cmd = "%s sort -@ %d %s/retained-pairs-remapped.bam -o %s/retained-pairs-remapped.sorted.bam" \
                        % (cmd_args.samtools, cmd_args.threads, bam_workspace, bam_workspace)
    print "Executing:", samtools_sort_cmd
    os.system(samtools_sort_cmd)

    read_categorizer_cmd = "./reads_categorizer %s %s %s" % (cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    print "Executing:", read_categorizer_cmd
    os.system(read_categorizer_cmd)

    bwa_cmd = "%s mem -t %d -h %d %s %s/virus-side.fq | %s view -b -F 2308 > %s/virus-side.bam" \
              % (cmd_args.bwa, cmd_args.threads, n_viruses, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    print "Executing:", bwa_cmd
    os.system(bwa_cmd)

    bwa_cmd = "%s mem -t %d -h 100 %s %s/host-side.fq | %s view -b -F 2308 > %s/host-side.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    print "Executing:", bwa_cmd
    os.system(bwa_cmd)

    bwa_cmd = "%s mem -t %d %s %s/virus-clips.fa | %s view -b -F 2308 > %s/virus-clips.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    print "Executing:", bwa_cmd
    os.system(bwa_cmd)

    bwa_cmd = "%s mem -t %d -h %d %s %s/host-clips.fa | %s view -b -F 2308 > %s/host-clips.bam" \
              % (cmd_args.bwa, cmd_args.threads, n_viruses, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    print "Executing:", bwa_cmd
    os.system(bwa_cmd)

    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-side.sorted.bam" % bam_workspace,
               "%s/host-side.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-side.sorted.bam" % bam_workspace,
               "%s/virus-side.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-anchors.sorted.bam" % bam_workspace,
               "%s/host-anchors.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-anchors.sorted.bam" % bam_workspace,
               "%s/virus-anchors.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-clips.sorted.bam" % bam_workspace,
               "%s/host-clips.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-clips.sorted.bam" % bam_workspace,
               "%s/virus-clips.bam" % bam_workspace)

    readsx = bam_workspace + "/readsx"
    if not os.path.exists(readsx):
        os.makedirs(readsx)

    remapper_cmd = "./remapper %s %s %s %s > %s/results.txt" \
                   % (cmd_args.host_reference, cmd_args.virus_reference, cmd_args.workdir, bam_workspace, bam_workspace)
    print "Executing:", remapper_cmd
    os.system(remapper_cmd)
