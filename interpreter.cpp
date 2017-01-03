#include <stdio.h>

#include "interpreter.h"

static char genotypeString[] = "AA\0\0AC\0\0AG\0\0AT\0\0CC\0\0CG\0\0CT\0\0GG\0\0GT\0\0TT\0\0AX\0\0CX\0\0GX\0\0TX\0\0XX\0\0XY\0\0NN\0\0";
static char sym[] = {'A', 'C', 'G', 'T', 'X', 'Y', 'N'};

void printHeader(FILE *fout) {
  fprintf(fout, "Chr\tPos\tRef\tPreBase\t");
  fprintf(fout, "pAA\tpAC\tpAG\tpAT\tpCC\tpCG\tpCT\tpGG\tpGT\tpTT\tpAX\tpCX\tpGX\tpTX\tpXX\tpXY\t");
  fprintf(fout, "GT\tOptb\tSuboptb\tpD\t");
  fprintf(fout, "rfProb\trfMode\t");
  fprintf(fout, "optIn/Del\tPatt\tHQ\tLQ\t");
  fprintf(fout, "suboptIn/Del\tPatt\tHQ\tLQ\t");
  fprintf(fout, "Depth\tA\tC\tG\tT\tA+\tA-\tC+\tC-\tG+\tG-\tT+\tT-\t");
  fprintf(fout, "AA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\tAX\tCX\tGX\tTX\tXX\tXY\t");
  fprintf(fout, "SB-A\tSB-C\tSB-G\tSB-T\t");
  fprintf(fout, "SB1-A\tSB1-C\tSB1-G\tSB1-T\t");
  fprintf(fout, "SB2-A\tSB2-C\tSB2-G\tSB2-T\t");
  fprintf(fout, "BQB-A\tBQB-C\tBQB-G\tBQB-T\t");
  fprintf(fout, "BQB1-A\tBQB1-C\tBQB1-G\tBQB1-T\t");
  fprintf(fout, "BQB2-A\tBQB2-C\tBQB2-G\tBQB2-T\t");
  fprintf(fout, "RPB\t");
  fprintf(fout, "LGC\tRGC\t");
  fprintf(fout, "LG\tRG\t");
  fprintf(fout, "ASC\t");
  fprintf(fout, "LQ\tRQ\t");
  fprintf(fout, "polyrun\t");
  fprintf(fout, "totalHQ\ttotalLQ\n"); 
  
}

__inline__ void printBase(FILE *fout, unsigned char base) {
  fprintf(fout, "%c", sym[base]);
}

__inline__ void printGenotype(FILE *fout, unsigned int genotype) {
  fprintf(fout, "%s\t", &(genotypeString[genotype << 2]));
}

__inline__ void printGenotypeLikelihood(FILE *fout, double *lh) {
  for (unsigned int i = 0; i < GENOTYPE_NUM; ++i) {
    fprintf(fout, "%.4e\t", lh[i]);
  }
}

__inline__ void printStrandBias(FILE *fout, StrandBias *sb) {
  for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
    fprintf(fout, "%.4e\t", sb->bias[i]);
  }
}

__inline__ void printBaseQualityBias(FILE *fout, BaseQualityBias *bqb) {
  for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
    fprintf(fout, "%.4e\t", bqb->bias[i]);
  }
}

__inline__ void printReadPositionBias(FILE *fout, DeltaStrandCount *rpb) {
  fprintf(fout, "%.4f\t", rpb->lrRatio);
}

__inline__ void printGCCount(FILE *fout, GCCount *gc) {
  fprintf(fout, "%.4f\t%.4f\t", gc->left, gc->right);
}

__inline__ void printGCount(FILE *fout, GCount *g) {
  fprintf(fout, "%u\t%u\t", g->left, g->right);
}

__inline__ void printAverageStrandCount(FILE *fout, AverageStrandCount *asc) {
  fprintf(fout, "%.4f\t", asc->lrRatio);
}

__inline__ void printReadQuality(FILE *fout, ReadQuality *rq) {
  fprintf(fout, "%.4f\t%.4f\t", rq->left, rq->right);
}

__inline__ void printIndel(FILE *fout, char indelType, unsigned char indelLength, unsigned char *pattern,
                           unsigned int hqCount, unsigned int lqCount) {
  if (indelType == '*') {
    fprintf(fout, "*\t*\t*\t*\t");
  } else {
    fprintf(fout, "%u%c\t", indelLength, indelType);
    for (unsigned int i = 0; i < indelLength; ++i) {
      printBase(fout, (pattern[i >> 2] >> ((i & 3) << 1)) & 3);
    }
    fprintf(fout, "\t");
    fprintf(fout, "%u\t%u\t", hqCount, lqCount);
  }
}

__inline__ void printSnpCounter(FILE *fout, MetaSnpCounter *snpCounter) {
  unsigned int depth = 0;
  for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
    depth += snpCounter->W[i];
  }
  fprintf(fout, "%u\t", depth);
  for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
    fprintf(fout, "%hu\t", snpCounter->W[i]);
  }
  for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
    fprintf(fout, "%hu\t%hu\t", snpCounter->F[i], snpCounter->R[i]);
  }
}

void
printPossibleSNP(MetaReference *reference, MetaSnpCounter *snpCounter, SnpCallingInfo *scInfo, unsigned char isIndel,
                 char indelType1, unsigned char indelLength1, unsigned char *indelPattern1, unsigned int hqCount1,
                 unsigned int lqCount1, char indelType2, unsigned char indelLength2, unsigned char *indelPattern2,
                 unsigned int hqCount2, unsigned int lqCount2, GenotypeLikelihood *genotypeLikelihood,
                 StrandBias *strandBias, BaseQualityBias *baseQualityBias, DeltaStrandCount *deltaStrandCount,
                 GCCount *gcCount, GCount *gCount, AverageStrandCount *avgStrandCount, ReadQuality *rQuality,
                 unsigned short polyrun, unsigned int indelHqCount, unsigned int indelLqCount, FILE *fout) {
  fprintf(fout, "%s\t", reference->chrName);
  fprintf(fout, "%u\t", reference->chrPos);
  fprintf(fout, "%c\t", sym[reference->refChar]);
  if (isIndel) {
    
    fprintf(fout, "%c\t", sym[reference->refChar]);
  } else {
    fprintf(fout, "*\t");
  }
  printGenotypeLikelihood(fout, scInfo->flh);
  printGenotype(fout, scInfo->genotype);
  fprintf(fout, "%.4e\t%.4e\t%.4e\t", scInfo->bestD, scInfo->secondD, scInfo->pD);
  fprintf(fout, "%.4e\t%d\t", scInfo->rfProb, scInfo->rfMode);
  printIndel(fout, indelType1, indelLength1, indelPattern1, hqCount1, lqCount1);
  printIndel(fout, indelType2, indelLength2, indelPattern2, hqCount2, lqCount2);
  printSnpCounter(fout, snpCounter);
  printGenotypeLikelihood(fout, genotypeLikelihood->likelihood);
  if (isIndel) {
    printStrandBias(fout, &(strandBias[0]));
    printStrandBias(fout, &(strandBias[1]));
    printStrandBias(fout, &(strandBias[2]));

    printBaseQualityBias(fout, &(baseQualityBias[0]));
    printBaseQualityBias(fout, &(baseQualityBias[1]));
    printBaseQualityBias(fout, &(baseQualityBias[2]));
  } else {
    printStrandBias(fout, &(strandBias[0]));
    fprintf(fout, "*\t*\t*\t*\t");
    fprintf(fout, "*\t*\t*\t*\t");

    printBaseQualityBias(fout, &(baseQualityBias[0]));
    fprintf(fout, "*\t*\t*\t*\t");
    fprintf(fout, "*\t*\t*\t*\t");
  }
  printReadPositionBias(fout, deltaStrandCount);
  printGCCount(fout, gcCount);
  printGCount(fout, gCount);
  printAverageStrandCount(fout, avgStrandCount);
  printReadQuality(fout, rQuality);
  fprintf(fout, "%hu\t", polyrun);
  fprintf(fout, "%u\t%u\n", indelHqCount, indelLqCount);
}
