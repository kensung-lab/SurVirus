# SurVirus

## Compiling

SurVirus has been compiled and tested with gcc 4.9.3, so we recommend this version or higher.

First of all, the required external libraries (downloaded with the source code) must be compiled with
```
./build_htslib.sh
```

Then, run
```
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

## Preparing the references

SurVirus needs three references:
1) the host genome
2) the virus(es) reference, one fasta sequence per virus
3) a concatenation of the host genome and the viruses genomes

Each of the references should be indexed with bwa and samtools. For example, suppose the host genome is contained in a file host.fa, and the virus genomes are in virus.fa. You should run
```
bwa index host.fa
samtools faidx host.fa

bwa index virus.fa
samtools faidx virus.fa

cat host.fa virus.fa > host+virus.fa
bwa index host+virus.fa
samtools faidx host+virus.fa
```

## Required software

Python 2 and libraries NumPy (http://www.numpy.org/), PyFaidx (https://github.com/mdshw5/pyfaidx) and PySam (https://github.com/pysam-developers/pysam) are required. 

Recent versions of samtools, bedtools, bwa and sdust are required. The latest versions at the moment of writing are 1.10 for samtools, 2.29 for bedtools, 0.7.17 for bwa.
For sdust, we recommend the implementation at https://github.com/lh3/sdust

## Running

The bare minimum command for running SurVirus is 
```
python surveyor input_files /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference 
```

input_files can be a list of comma-separated bam_files, for example
```
python surveyor	input1.bam,input2.bam /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference
```
Note that if multiple BAM files are present, they must all be relative to the same sample.

input_files can also be a pair of comma-separated fastq files, containing read pairs
```
python surveyor reads_1.fq,reads_2.fq /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference --fq
```
Note that in this second case, the flag --fq must be provided.
In the next release we will target support for multiple pairs of fastq files.

If data is whole-genome sequencing, the flag --wgs should be provided. This is not mandatory (i.e. SurVirus will run anyway), but recommended.

If samtools, bedtools, bwa or sdust are not in your PATH, or if you wish to provide an alternative location for either of them, you can do so with the --samtools, --bedtools, --bwa and --dust flags
```
python surveyor input_files /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference --samtools /path/to/samtools --bedtools /path/to/bedtools --bwa /path/to/bwa --dust /path/to/sdust
```

Finally, an useful flag is --threads, which you can use to specify the number of threads you wish to assign to SurVirus. The default is 1.

## Output

The final output will be placed in the workdir, under the name results.remapped.t1.txt.

The following line is an example of a predicted integration:
```
ID=0 chr13:-73788865 type16:-3383 SUPPORTING_PAIRS=700 SPLIT_READS=725 HOST_PBS=0.927897 COVERAGE=0.684993
```

The ID is simply a unique number. The second and the third fields are the host and the virus coordinates for the integration. In this example, the sample contains a junction between chr13:73788865, negative strand, and type16:3383, negative strand.
Supporting pairs and split reads are self-explanatory, the first is the number of read pairs that supports the integration while the second is the number of split reads that overlap the junction (i.e. they map partly to the host breakpoint and partly to the virus breakpoint).

HOST_PBS and COVERAGE are quality metrics. Since they are used in the filtering of the results, the user can safely ignore them most of the time. 
HOST_PBS is the fraction of base matches between the supporting reads that are mapped to the host breakpoint and the reference sequence. In the example, a value of 0.927897 means that when performing a local alignment between the supporting reads and the host reference sequence near the breakpoint, we produce ~92.8% base matches (as opposed to mismatches and gaps).
SurVirus internally analyses the distribution of insert sizes in the input sample, and determines the maximum insert size (maxIS) that is considered not to be an artifact (i.e. what is usually referred to as discordant by most SV callers). When remapping reads, SurVirus considers a maxIS bp-long window next to the suspected breakpoint, for both virus and host. The rationale behind this is that if a read is more than maxIS bp away from a breakpoint, its pair would not be able to be chimeric, as the mate would not be able to cross the breakpoint.
COVERAGE is the fraction of such maxIS window that is covered by reads, averaged between host and virus.
