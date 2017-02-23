# 16GT
16GT is a variant caller utilizing a 16-genotype probabilistic model to unify SNP and indel calling in a single algorithm.

## Quick start
> Inputs: genome.fa alignments.bam, Output: variants.vcf

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
#### 3. Call variants
```
snapshotSnpcaller -i genome.fa.index -o output/prefix
perl txt2vcf.pl output/prefix.txt sampleName genome.fa > variants.vcf
perl filterVCF.pl variants.vcf > variants.filtered.vcf
```

## Exome variant calling
> Inputs: genome.fa alignement.bam region.bed, Outputs: region.bin variants.vcf

```
RegionIndexBuilder genome.fa.index region.bed region.bin -bed/-gff
bam2snapshot -i genome.fa.index -b alignments.bam -o output/prefix -e region.bin
snapshotSnpcaller -i genome.fa.index -o output/prefix -e region.bin
```

## License
GPLv3

