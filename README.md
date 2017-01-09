# 16GT
16GT is a variant caller utilizing a 16-genotype probabilistic model to unify SNP and indel calling in a single algorithm.

## Quick start
> Inputs: genome.fa alignments.bam
> Output: variants.vcf

#### 0. Install
```
git clone https://github.com/aquaskyline/16GT
make
```
#### 1. Build reference index (optional)
```
git clone https://github.com/aquaskyline/SOAP3-dp.git
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
```

## License
GPLv3

