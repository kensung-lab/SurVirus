import argparse, os
import pysam, pyfaidx
import max_is_calc

cmd_parser = argparse.ArgumentParser(description='SurVirus, a virus integration caller.')
cmd_parser.add_argument('input_files', help='Input files, separated by a comma.')
cmd_parser.add_argument('workdir', help='Working directory for SurVirus to use.')
cmd_parser.add_argument('host_reference', help='Reference of the host organism in FASTA format.')
cmd_parser.add_argument('virus_reference', help='References of a list of viruses in FASTA format.')
cmd_parser.add_argument('host_and_virus_reference', help='Joint references of host and viruses.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--bwa', default='bwa', help='BWA path.')
cmd_parser.add_argument('--samtools', help='Samtools path.', default='samtools')
cmd_parser.add_argument('--dust', help='Dust path.', default='dust')
cmd_parser.add_argument('--wgs', action='store_true', help='The reference genome is uniformly covered by reads.'
                                                           'SurVirus needs to sample read pairs, and this option lets'
                                                           'it sample them from all over the genome.')
cmd_parser.add_argument('--covered_regions_bed', default='',
                        help='In case only a few regions of the genome are covered by reads, this directs SurVirus'
                             'on where to sample reads. In case --wgs is not used and this is not provided,'
                             'SurVirus will compute such file by itself.')
cmd_parser.add_argument('--minClipSize', type=int, default=20, help='Min size for a clip to be considered.')
cmd_parser.add_argument('--maxSCDist', type=int, default=10, help='Max SC distance.')
cmd_parser.add_argument('--fq', action='store_true', help='Input is in fastq format.')
cmd_parser.add_argument('--cram-reference', help='Can optionally provide a reference for decoding the input file(s) if '
                                                 'in CRAM.')
cmd_args = cmd_parser.parse_args()

def execute(cmd):
    print "Executing:", cmd
    os.system(cmd)

input_names = cmd_args.input_files.split(',')
if cmd_args.fq:
    if len(input_names) != 2: # TODO: accept multiple pairs of fq files
        print "Two colon-separated fastq files are required."
        exit(1)

# Create config file in workdir
config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("min_sc_size %d\n" % cmd_args.minClipSize)
config_file.write("max_sc_dist %d\n" % cmd_args.maxSCDist)
if cmd_args.cram_reference:
    config_file.write("cram_reference %s\n" % cmd_args.cram_reference)

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

bam_workspaces = []
if cmd_args.fq:
    bam_workspace = "%s/bam_0/" % (cmd_args.workdir)
    if not os.path.exists(bam_workspace):
        os.makedirs(bam_workspace)
    bam_workspaces.append(bam_workspace)

    max_read_len, max_is = \
        max_is_calc.get_max_is_from_fq(cmd_args.workdir, input_names[0], input_names[1], cmd_args.host_and_virus_reference, \
                                            cmd_args.bwa, cmd_args.threads)
    with open("%s/stats.txt" % bam_workspace, "w") as stat_file:
        stat_file.write("max_is %d\n" % max_is)
    config_file.write("read_len %d\n" % max_read_len)
    config_file.close();

    isolate_cmd = "./isolate_relevant_pairs_fq %s %s %s %s %s %s " % (input_names[0], input_names[1], cmd_args.host_reference,
                                                        cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    execute(isolate_cmd)
else:
    if cmd_args.cram_reference:
        bam_files = [pysam.AlignmentFile(bam_name, reference_filename=cmd_args.cram_reference) for bam_name in input_names]
    else:
        bam_files = [pysam.AlignmentFile(bam_name) for bam_name in input_names]
    max_read_len, max_is_list = \
        max_is_calc.get_max_is_from_bam(cmd_args.host_reference, bam_files, cmd_args.wgs, cmd_args.covered_regions_bed)
    config_file.write("read_len %d\n" % max_read_len)
    config_file.close();

    for file_index, bam_file in enumerate(bam_files):
        bam_workspace = "%s/bam_%d/" % (cmd_args.workdir, file_index)
        if not os.path.exists(bam_workspace):
            os.makedirs(bam_workspace)
        bam_workspaces.append(bam_workspace)

        max_is = max_is_list[file_index]
        with open("%s/stats.txt" % bam_workspace, "w") as stat_file:
            stat_file.write("max_is %d\n" % max_is)

        isolate_cmd = "./isolate_relevant_pairs %s %s %s %s %s" % (bam_file.filename, cmd_args.host_reference,
                                                                   cmd_args.virus_reference, cmd_args.workdir,
                                                                   bam_workspace)
        execute(isolate_cmd)

        filter_by_qname_cmd = "./filter_by_qname %s %s %s" % (bam_file.filename, cmd_args.workdir, bam_workspace)
        execute(filter_by_qname_cmd)


def map_clips(prefix, reference):
    bwa_aln_cmd = "%s aln -t %d %s %s.fa -f %s.sai" \
                  % (cmd_args.bwa, cmd_args.threads, reference, prefix, prefix)
    bwa_samse_cmd = "%s samse %s %s.sai %s.fa | %s view -b -F 2304 > %s.full.bam" \
                    % (cmd_args.bwa, reference, prefix, prefix, cmd_args.samtools, prefix)
    execute(bwa_aln_cmd)
    execute(bwa_samse_cmd)

    filter_unmapped_cmd = "%s view -b -F 4 %s.full.bam > %s.aln.bam" \
                          % (cmd_args.samtools, prefix, prefix)
    execute(filter_unmapped_cmd)

    dump_unmapped_fa = "%s fasta -f 4 %s.full.bam > %s.unmapped.fa" \
                       % (cmd_args.samtools, prefix, prefix)
    execute(dump_unmapped_fa)

    bwa_mem_cmd = "%s mem -t %d %s %s.unmapped.fa | %s view -b -F 2308 > %s.mem.bam" \
                  % (cmd_args.bwa, cmd_args.threads, reference, prefix,
                     cmd_args.samtools, prefix)
    execute(bwa_mem_cmd)

    cat_cmd = "%s cat %s.aln.bam %s.mem.bam -o %s.bam" \
              % (cmd_args.samtools, prefix, prefix, prefix)
    execute(cat_cmd)

    pysam.sort("-@", str(cmd_args.threads), "-o", "%s.cs.bam" % prefix, "%s.bam" % prefix)

for bam_workspace in bam_workspaces:
    bwa_cmd = "%s mem -t %d %s %s/retained-pairs_1.fq %s/retained-pairs_2.fq | %s view -b -F 2304 > %s/retained-pairs.remapped.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, \
                 bam_workspace, bam_workspace, cmd_args.samtools, bam_workspace)
    execute(bwa_cmd)

    samtools_sort_cmd = "%s sort -@ %d %s/retained-pairs.remapped.bam -o %s/retained-pairs.remapped.cs.bam" \
                        % (cmd_args.samtools, cmd_args.threads, bam_workspace, bam_workspace)
    execute(samtools_sort_cmd)

    extract_clips_cmd = "./extract_clips %s %s %s" % (cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    execute(extract_clips_cmd)

    # map virus clips
    map_clips("%s/virus-clips" % bam_workspace, cmd_args.host_reference)

    # map host clips
    map_clips("%s/host-clips" % bam_workspace, cmd_args.virus_reference)

    read_categorizer_cmd = "./reads_categorizer %s %s %s" % (cmd_args.virus_reference, cmd_args.workdir, bam_workspace)
    execute(read_categorizer_cmd)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-anchors.cs.bam" % bam_workspace,
             "%s/host-anchors.bam" % bam_workspace)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-anchors.cs.bam" % bam_workspace,
               "%s/virus-anchors.bam" % bam_workspace)

    bwa_cmd = "%s mem -t %d -h %d %s %s/virus-side.fq | %s view -b -F 2308 > %s/virus-side.bam" \
              % (cmd_args.bwa, cmd_args.threads, n_viruses, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    execute(bwa_cmd)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/virus-side.cs.bam" % bam_workspace, "%s/virus-side.bam" % bam_workspace)

    bwa_cmd = "%s mem -t %d -h 100 %s %s/host-side.fq | %s view -b -F 2308 > %s/host-side.bam" \
              % (cmd_args.bwa, cmd_args.threads, cmd_args.host_and_virus_reference, bam_workspace,
                 cmd_args.samtools, bam_workspace)
    execute(bwa_cmd)
    pysam.sort("-@", str(cmd_args.threads), "-o", "%s/host-side.cs.bam" % bam_workspace, "%s/host-side.bam" % bam_workspace)

readsx = cmd_args.workdir + "/readsx"
if not os.path.exists(readsx):
    os.makedirs(readsx)

bam_workspaces_str = " ".join(bam_workspaces)
merge_retained_reads_cmd = "./merge_retained_reads %s %s" % (cmd_args.workdir, bam_workspaces_str)
execute(merge_retained_reads_cmd)

build_rr_associations_cmd = "./build_region-reads_associations %s %s %s %s" \
                            % (cmd_args.host_reference, cmd_args.virus_reference, cmd_args.workdir, bam_workspaces_str)
execute(build_rr_associations_cmd)

remapper_cmd = "./remapper %s %s %s %s > %s/results.txt 2> %s/log.txt" \
               % (cmd_args.host_reference, cmd_args.virus_reference, cmd_args.workdir,
                  bam_workspaces_str, cmd_args.workdir, cmd_args.workdir)
execute(remapper_cmd)

with open(cmd_args.workdir + "/results.txt") as results_file:
    for line in results_file:
        id = line.split()[0]
        bam_prefix = "%s/%s" % (readsx, id)
        pysam.sort("-@", str(cmd_args.threads), "-o", "%s.cs.bam" % bam_prefix, "%s.bam" % bam_prefix)
        os.rename("%s.cs.bam" % bam_prefix, "%s.bam" % bam_prefix)
        execute("%s index %s" % (cmd_args.samtools, "%s.bam" % bam_prefix))

bp_consensus_cmd = "./bp_region_consensus_builder %s %s %s %s" \
                   % (cmd_args.host_reference, cmd_args.virus_reference, cmd_args.workdir, bam_workspaces_str)
execute(bp_consensus_cmd)

dust_cmd = "%s %s/host_bp_seqs.fa > %s/host_bp_seqs.masked.bed" % (cmd_args.dust, cmd_args.workdir, cmd_args.workdir)
execute(dust_cmd)

dust_cmd = "%s %s/virus_bp_seqs.fa > %s/virus_bp_seqs.masked.bed" % (cmd_args.dust, cmd_args.workdir, cmd_args.workdir)
execute(dust_cmd)

filter_cmd = "./filter %s > %s/results.t1.txt" % (cmd_args.workdir, cmd_args.workdir)
execute(filter_cmd)

filter_cmd = "./filter %s --remapped > %s/results.remapped.t1.txt" % (cmd_args.workdir, cmd_args.workdir)
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
        id, dir = r.qname.split("_")
        is_rev = r.is_reverse and dir == "R" or not r.is_reverse and dir == "L"
        pos = r.reference_start if is_rev else r.reference_end
        altf.write("ID=%s %s:%c%d\n" % (id, r.reference_name, "-" if is_rev else "+", pos))
        for xa_region in xa_regions:
            xa_region_split = xa_region.split(',')
            pos = int(xa_region_split[1])
            is_rev = pos < 0 and dir == "R" or pos > 0 and dir == "L"
            pos = abs(pos)
            if not is_rev: pos += r.query_length
            altf.write("ID=%s %s:%c%d\n" % (id, xa_region_split[0], "-" if is_rev else "+", pos))

