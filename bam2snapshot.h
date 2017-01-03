#ifndef _BAM2SNAPSHOT_H
#define _BAM2SNAPSHOT_H

#include <assert.h>
#include <time.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "SNP.h"
#include "definitions.h"
#include "dependencies.h"
#include "SnapshotHandler.h"
#include "CounterReader.h"
#include "sam.h"
#include "indexFunction.h"

#ifndef BAM_CIGAR_STR
#define BAM_CIGAR_STR  "MIDNSHP=XB"
#endif

#ifndef bam_cigar_op
#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#endif

#ifndef bam_cigar_oplen
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#endif

#ifndef bam_cigar_opchr
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#endif

#define Max_Num_Of_Chrs  (100)
#define Max_Read_Length  (65535)
#define Max_Num_Of_Input_Files (65536)

typedef struct AlignmentQC
{
  unsigned long long totalReads;
  unsigned long long unalignedReads;
  unsigned long long softclippedReads;
  unsigned long long totalReadLength;
  unsigned long long totalAlignedLength;
  unsigned long long targetAlignedLength;
} AlignmentQC;

void resetAlignmentQC(AlignmentQC &alignmentQC);
void updateAlignmentQCForReads(AlignmentQC &alignmentQC, int readLength, char *cigar, unsigned int startPos,
                               char isAlign, unsigned int *exomePosBitVec);
void writeAlignmentQCToFile(AlignmentQC alignmentQC, char *outputFilePrefix);

static const int bamcode2snpcode[16] =
{
  0, 0, 1, 0,
  2, 0, 0, 0,
  3, 0, 0, 0,
  0, 0, 0, 2
};

struct Alignment
{
  unsigned startPos;
  short strand;
  short readLength;
  unsigned char *read;
  char *qualities, *cigar;
};

void initAlignment(Alignment *align, int maxReadLength, int maxCigarLength);
void freeAlignment(Alignment *align);
void clearSnapshot(unsigned dnaLength, unsigned numOfCPUThreads, SnpBundle *snpBundle);
void freeSnapshot(SnpBundle *snpBundle);
void bamFilesToSnapshot(SnpBundle *snpBundle, unsigned dnaLength,
                        char *inputFilenames[], int inputFileTypes[], int n_inputFiles, char *outputFilenamePrefix,
                        unsigned char recalScore[], unsigned char weightMap[], unsigned snpTrimSize,
                        unsigned softClipThreshold,
                        unsigned mapqThreshold, unsigned indelWeightThreshold,
                        Annotation *annotations, int numOfSeq, unsigned int *ambMap, Translate *translate,
                        AlignmentQC &alignmentQC, unsigned int *exomePosBitVec,
                        const time_t *startTime);

#endif
