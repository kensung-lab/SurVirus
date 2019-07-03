import argparse, os
import pysam, pyfaidx
from random_pos_generator import RandomPositionGenerator
import numpy as np

MAX_READS = 1000
GEN_DIST_SIZE = 10000
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
cmd_parser.add_argument('--bedtools', default='bedtools', help='Bedtools path. Necessary if --wgs is not set and'
                                                               '--coveredRegionsBed is not provided.')
cmd_parser.add_argument('--wgs', action='store_true', help='The reference genome is uniformly covered by reads.'
                                                           'SurVirus needs to sample read pairs, and this option lets'
                                                           'it sample them from all over the genome.')
cmd_parser.add_argument('--coveredRegionsBed', default='',
                        help='In case only a few regions of the genome are covered by reads, this directs SurVirus'
                             'on where to sample reads. In case --wgs is not used and this is not provided,'
                             'SurVirus will compute such file by itself.')
cmd_parser.add_argument('--minClipSize', type=int, default=30, help='Min size for a clip to be considered.')
cmd_parser.add_argument('--maxSCDist', type=int, default=10, help='Max SC distance.')
cmd_args = cmd_parser.parse_args()

bam_names = cmd_args.bam_files.split(';')
bam_files = [pysam.AlignmentFile(bam_name) for bam_name in bam_names]


# Create config file in workdir

config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("min_sc_size %d\n" % cmd_args.minClipSize)
config_file.write("max_sc_dist %d\n" % cmd_args.maxSCDist)
config_file.close();


def execute(cmd):
    print "Executing:", cmd
    os.system(cmd)


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


sampling_regions = []
if cmd_args.wgs:
    reference_fa = pyfaidx.Fasta(cmd_args.host_reference)
    for k in reference_fa.keys():
        sampling_regions += [(k, 1, len(reference_fa[k]))]
else:
    sampling_regions_file = cmd_args.coveredRegionsBed
    if not sampling_regions_file:
        bamtobed_cmd = "%s bamtobed -i %s | cut -f-3 | uniq > %s/reads.bed" % (cmd_args.bedtools, bam_names[0], cmd_args.workdir)
        print "Executing:", bamtobed_cmd
        os.system(bamtobed_cmd)

        bedmerge_cmd = "%s merge -i %s/reads.bed -c 1 -o count > %s/reads.merged.bed" \
                       % (cmd_args.bedtools, cmd_args.workdir, cmd_args.workdir)
        print "Executing:", bedmerge_cmd
        os.system(bedmerge_cmd)

        sampling_regions_file = "%s/sampling_regions.bed" % cmd_args.workdir
        with open(sampling_regions_file, "w") as sr_outf, open("%s/reads.merged.bed" % cmd_args.workdir) as mr_inf:
            for line in mr_inf:
                sl = line.split()
                if int(sl[3]) >= 20:
                    sr_outf.write("%s\t%s\t%s\n" % tuple(sl[:3]))

    with open(sampling_regions_file) as sr_inf:
        for line in sr_inf:
            sl = line.split()
            sampling_regions += [(sl[0], int(sl[1]), int(sl[2]))]

rand_pos_gen = RandomPositionGenerator(sampling_regions)

print "Generating random genomic positions..."
random_positions = []
for i in range(1,GEN_DIST_SIZE*2):
    if i % 100000 == 0: print i, "random positions generated."
    random_positions.append(rand_pos_gen.next())

for file_index, bam_file in enumerate(bam_files):
    read_len = 0
    for i, read in enumerate(bam_file.fetch(until_eof=True)):
        if i > MAX_READS: break
        read_len = max(read_len, read.query_length)

    general_dist = []
    avg_depth = 0
    samplings = 0
    rnd_i = 0
    for e in random_positions:
        if len(general_dist) > GEN_DIST_SIZE: break

        print len(general_dist)

        chr, pos = e
        rnd_i += 1

        samplings += 1
        i = 0
        used = set()
        for read in bam_file.fetch(contig=chr, start=pos, end=pos+10000):
            if read.reference_start < pos:
                avg_depth += 1
            if (read.reference_start, read.next_reference_start) in used: continue # try not to use duplicates

            if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
            0 < read.template_length < MAX_ACCEPTABLE_IS and 'S' not in read.cigarstring and 'S' not in read.get_tag('MC'):
                if i > 100: break
                i += 1
                general_dist.append(read.template_length)
                used.add((read.reference_start, read.next_reference_start))

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


print "Regenerating random genomic positions..."
random_positions = []
for i in range(1,GEN_DIST_SIZE*2):
    if i % 100000 == 0: print i, "random positions generated."
    random_positions.append(rand_pos_gen.next())

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)


bam_workspaces = []
for file_index, bam_file in enumerate(bam_files):
    bam_workspace = "%s/bam_%d/" % (cmd_args.workdir, file_index)
    bam_workspaces.append(bam_workspace)

    isolate_cmd = "./isolate_relevant_pairs %s %s %s %s %s" % (bam_file.filename, cmd_args.host_reference,
                                                            cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    execute(isolate_cmd)

    filter_by_qname_cmd = "./filter_by_qname %s %s %s" % (bam_file.filename, cmd_args.workdir, bam_workspace)
    execute(filter_by_qname_cmd)

    samtools_sortbyqname_cmd = "%s sort -n -@ %d -o %s/retained-pairs.qname-sorted.bam %s/retained-pairs.bam" % \
                               (cmd_args.samtools, cmd_args.threads, bam_workspace, bam_workspace)
    execute(samtools_sortbyqname_cmd)

    samtools_dump_cmd = "%s fastq -1 %s/retained-pairs_1.fq -2 %s/retained-pairs_2.fq %s/retained-pairs.qname-sorted.bam" % \
                        (cmd_args.samtools, bam_workspace, bam_workspace, bam_workspace)
    execute(samtools_dump_cmd)

    bwa_cmd = "%s mem -t %d %s %s/retained-pairs_1.fq %s/retained-pairs_2.fq | %s view -b -F 2304 > %s/retained-pairs-remapped.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, \
                 bam_workspace, bam_workspace, cmd_args.samtools, bam_workspace)
    execute(bwa_cmd)

    samtools_sort_cmd = "%s sort -@ %d %s/retained-pairs-remapped.bam -o %s/retained-pairs-remapped.sorted.bam" \
                        % (cmd_args.samtools, cmd_args.threads, bam_workspace, bam_workspace)
    execute(samtools_sort_cmd)

    extract_clips_cmd = "./extract_clips %s %s %s" % (cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    execute(extract_clips_cmd)

    # map virus clips
    bwa_aln_cmd = "%s aln -t %d %s %s/virus-clips.fa -f %s/virus-clips.sai" \
                  % (cmd_args.bwa, cmd_args.threads, cmd_args.host_reference, bam_workspace, bam_workspace)
    bwa_samse_cmd = "%s samse %s %s/virus-clips.sai %s/virus-clips.fa | %s view -b -F 2304 > %s/virus-clips.full.bam" \
                    % (cmd_args.bwa, cmd_args.host_reference, bam_workspace, bam_workspace, cmd_args.samtools, bam_workspace)
    execute(bwa_aln_cmd)
    execute(bwa_samse_cmd)

    filter_unmapped_cmd = "%s view -b -F 4 %s/virus-clips.full.bam > %s/virus-clips.aln.bam" \
                          % (cmd_args.samtools, bam_workspace, bam_workspace)
    execute(filter_unmapped_cmd)

    dump_unmapped_fa = "%s fasta -f 4 %s/virus-clips.full.bam > %s/virus-clips.unmapped.fa" \
                       % (cmd_args.samtools, bam_workspace, bam_workspace)
    execute(dump_unmapped_fa)

    bwa_mem_cmd = "%s mem -t %d %s %s/virus-clips.unmapped.fa | %s view -b -F 2308 > %s/virus-clips.mem.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    execute(bwa_mem_cmd)

    cat_cmd = "%s cat %s/virus-clips.aln.bam %s/virus-clips.mem.bam -o %s/virus-clips.bam" \
              % (cmd_args.samtools, bam_workspace, bam_workspace, bam_workspace)
    execute(cat_cmd)

    # map host clips
    bwa_aln_cmd = "%s aln -t %d %s %s/host-clips.fa -f %s/host-clips.sai" \
                  % (cmd_args.bwa, cmd_args.threads, cmd_args.virus_reference, bam_workspace, bam_workspace)
    bwa_samse_cmd = "%s samse %s %s/host-clips.sai %s/host-clips.fa | %s view -b -F 2304 > %s/host-clips.full.bam" \
                    % (cmd_args.bwa, cmd_args.virus_reference, bam_workspace, bam_workspace, cmd_args.samtools,
                       bam_workspace)
    execute(bwa_aln_cmd)
    execute(bwa_samse_cmd)

    filter_unmapped_cmd = "%s view -b -F 4 %s/host-clips.full.bam > %s/host-clips.aln.bam" \
                          % (cmd_args.samtools, bam_workspace, bam_workspace)
    execute(filter_unmapped_cmd)

    dump_unmapped_fa = "%s fasta -f 4 %s/host-clips.full.bam > %s/host-clips.unmapped.fa" \
                       % (cmd_args.samtools, bam_workspace, bam_workspace)
    execute(dump_unmapped_fa)

    bwa_mem_cmd = "%s mem -t %d %s %s/host-clips.unmapped.fa | %s view -b -F 2308 > %s/host-clips.mem.bam" \
                  % (cmd_args.bwa, cmd_args.threads, cmd_args.virus_reference, bam_workspace,
                     cmd_args.samtools, bam_workspace)
    execute(bwa_mem_cmd)

    cat_cmd = "%s cat %s/host-clips.aln.bam %s/host-clips.mem.bam -o %s/host-clips.bam" \
              % (cmd_args.samtools, bam_workspace, bam_workspace, bam_workspace)
    execute(cat_cmd)

    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-clips.sorted.bam" % bam_workspace,
               "%s/host-clips.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-clips.sorted.bam" % bam_workspace,
               "%s/virus-clips.bam" % bam_workspace)

    read_categorizer_cmd = "./reads_categorizer %s %s %s" % (cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    execute(read_categorizer_cmd)

    bwa_cmd = "%s mem -t %d -h %d %s %s/virus-side.fq | %s view -b -F 2308 > %s/virus-side.bam" \
              % (cmd_args.bwa, cmd_args.threads, n_viruses, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    execute(bwa_cmd)

    bwa_cmd = "%s mem -t %d -h 100 %s %s/host-side.fq | %s view -b -F 2308 > %s/host-side.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    execute(bwa_cmd)

    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-side.sorted.bam" % bam_workspace,
               "%s/host-side.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-side.sorted.bam" % bam_workspace,
               "%s/virus-side.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-anchors.sorted.bam" % bam_workspace,
             "%s/host-anchors.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-anchors.sorted.bam" % bam_workspace,
               "%s/virus-anchors.bam" % bam_workspace)

readsx = cmd_args.workdir + "/readsx"
if not os.path.exists(readsx):
    os.makedirs(readsx)


bam_files_and_workspaces = " ".join([bam_files[i].filename + " " + bam_workspaces[i] for i in xrange(len(bam_workspaces))])
remapper_cmd = "./remapper %s %s %s %s > %s/results.txt 2> %s/log.txt" \
               % (cmd_args.host_reference, cmd_args.virus_reference, cmd_args.workdir,
                  bam_files_and_workspaces, cmd_args.workdir, cmd_args.workdir)
execute(remapper_cmd)

filter_cmd = "./filter %s > %s/results.t1.txt" % (cmd_args.workdir, cmd_args.workdir)
execute(filter_cmd)

filter_cmd = "./filter %s --print-rejected > %s/results.discarded.txt" % (cmd_args.workdir, cmd_args.workdir)
execute(filter_cmd)


print "Finding alternative locations..."

bwa_cmd = "%s mem -t %d -h 1000 %s %s/host_bp_seqs.fa | %s view -b -F 2308 > %s/host_bp_seqs.bam" \
          % (cmd_args.bwa, cmd_args.threads, cmd_args.host_reference, cmd_args.workdir, cmd_args.samtools, cmd_args.workdir)
execute(bwa_cmd)

with pysam.AlignmentFile("%s/host_bp_seqs.bam" % cmd_args.workdir) as bp_seqs_bam, \
        open("%s/results.alternative.txt" % cmd_args.workdir, 'w') as altf:
    for r in bp_seqs_bam.fetch(until_eof=True):
        if not r.has_tag('XA'): continue

        xa_regions = r.get_tag('XA').split(';')[:-1]
        pos = r.reference_start if r.is_reverse else r.reference_end
        altf.write("ID=%s %s:%c%d\n" % (r.qname, r.reference_name, "-" if r.is_reverse else "+", pos))
        for xa_region in xa_regions:
            xa_region_split = xa_region.split(',')
            pos = int(xa_region_split[1])
            if pos > 0: pos += r.query_length
            altf.write("ID=%s %s:%c%d\n" % (r.qname, xa_region_split[0], "+" if pos > 0 else "-", abs(pos)))

