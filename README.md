# Setup

## docker
```
git clone https://github.com/aquaskyline/16GT.git
cd 16GT
docker build --no-cache .
docker images
```

## use the respective "IMAGE ID" displayed above as  <docker-id> below
```
docker run -it --privileged <docker-id> /bin/bash
```

## once inside the docker image, index the reference
```
cd /16GT/SOAP3-dp
./soap3-dp-builder <path-to-ref-gen-fasta>
./BGS-Build <path-to-ref-gen-fasta>.index
```


## variant call using aligned/indexed bam file
```
cd /16GT
./bam2snapshot -i <path-to-ref-gen-fasta>.index -b <aligned-bam-file> -o <output-prefix>
./snapshotSnpcaller  -i <path-to-ref-gen-fasta>.index  -o <output-prefix>
perl txt2vcf.pl <output-prefix>.txt <pro-id> <path-to-ref-gen-fasta> > <output>.vcf
perl filterVCF.pl <output>.vcf > <output>.filtered.vcf
```

# 16GT
16GT is a variant caller utilizing a 16-genotype probabilistic model to unify SNP and indel calling in a single algorithm.
16GT is easy to use. The default parameters will fit most of the use cases with human genome.
For the detailed parameters for each module, please run the module to get an info.

## Quick start
> Inputs: genome.fa alignments.bam, Output: <output>.vcf

#### 0. Install
```
git clone https://github.com/aquaskyline/16GT
cd 16GT
make
# Tested in Ubuntu 14.04 and CentOS 6.7 with GCC 4.7.2
```
#### 1. Build reference index
```
git clone https://github.com/aquaskyline/SOAP3-dp.git
cd SOAP3-dp
make SOAP3-Builder
make BGS-Build
soap3-dp-builder genome.fa
BGS-Build genome.fa.index
```
#### 2. Convert BAM to SNAPSHOT
```
bam2snapshot -i genome.fa.index -b alignments.bam -o output/prefix
```
#### 3. Call <output>
```
snapshotSnpcaller -i genome.fa.index -o output/prefix
perl txt2vcf.pl output/prefix.txt sampleName genome.fa > <output>.vcf
perl filterVCF.pl <output>.vcf dbSNP.vcf.gz > <output>.filtered.vcf
```

## Exome variant calling
> Inputs: genome.fa alignement.bam region.bed, Outputs: region.bin <output>.vcf

```
RegionIndexBuilder genome.fa.index region.bed region.bin -bed/-gff
bam2snapshot -i genome.fa.index -b alignments.bam -o output/prefix -e region.bin
snapshotSnpcaller -i genome.fa.index -o output/prefix -e region.bin
```

## License
GPLv3

