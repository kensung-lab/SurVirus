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


