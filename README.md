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

Each of the references should be indexed with bwa and samtools, i.e.
```
bwa index ref.fa
samtools faidx ref.fa
```

## Required software

Python 2 and libraries NumPy (http://www.numpy.org/), PyFaidx (https://github.com/mdshw5/pyfaidx) and PySam (https://github.com/pysam-developers/pysam) are required. 

Recent versions of samtools, bedtools, bwa and sdust are required. The latest versions at the moment of writing are 1.10 for samtools, 2.29 for bedtools, 0.7.17 for bwa.
For sdust, we recommend the implementation at https://github.com/lh3/sdust

## Running

The bare minimum command for running SurVirus is 
```
python surveyor input_files /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference 
```

input_files can be a list of comma-separated bam_files, for example
```
python surveyor	input1.bam,input2.bam /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference
```
Note that if multiple BAM files are present, they must all be relative to the same sample.

input_files can also be a pair of comma-separated fastq files, containing read pairs
```
python surveyor reads_1.fq,reads_2.fq /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference --fq
```
Note that in this second case, the flag --fq must be provided.
In the next release we will target support for multiple pairs of fastq files.

If data is whole-genome sequencing, the flag --wgs should be provided. This is not mandatory (i.e. SurVirus will run anyway), but recommended.
