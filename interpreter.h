#ifndef _INTERPRETER_H_
#define _INTERPRETER_H_

#include <stdio.h>
#include "SNP_Meta.h"

void printHeader(FILE *fout);

void
printPossibleSNP(MetaReference *reference, MetaSnpCounter *snpCounter, SnpCallingInfo *scInfo, unsigned char isIndel,
                 char indelType1, unsigned char indelLength1, unsigned char *indelPattern1, unsigned int hqCount1,
                 unsigned int lqCount1, char indelType2, unsigned char indelLength2, unsigned char *indelPattern2,
                 unsigned int hqCount2, unsigned int lqCount2, GenotypeLikelihood *genotypeLikelihood,
                 StrandBias *strandBias, BaseQualityBias *baseQualityBias, DeltaStrandCount *deltaStrandCount,
                 GCCount *gcCount, GCount *gCount, AverageStrandCount *avgStrandCount, ReadQuality *rQuality,
                 unsigned short polyrun, unsigned int indelHqCount, unsigned int indelLqCount, FILE *fout);

#endif
