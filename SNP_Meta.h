#ifndef _SNP_META_H_
#define _SNP_META_H_

#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"

#define ALPHABET_SIZE 4
#define MAX_SEQ_NAME_LENGTH 256
#define SNP_META_WINDOW_SIZE (100)
#define MAX_PATTERN_LENGTH 1024


#define GENOTYPE_16

#ifdef GENOTYPE_10
#define GENOTYPE_NUM 10
#define COUNTER_NUM (ALPHABET_SIZE)
#endif
#ifdef GENOTYPE_16
#define GENOTYPE_NUM 16
#define COUNTER_NUM (ALPHABET_SIZE + 2)
#endif

#define UNIDENTIFIED_GENOTYPE GENOTYPE_NUM

typedef struct InputStat
{
  unsigned int W[COUNTER_NUM];
  unsigned int F[COUNTER_NUM];
  unsigned int R[COUNTER_NUM];


  unsigned int outW[COUNTER_NUM];
  unsigned int outF[COUNTER_NUM];
  unsigned int outR[COUNTER_NUM];
} InputStat;

typedef struct GenotypeInputStat
{
  unsigned short W[COUNTER_NUM];
} GenotypeInputStat;

typedef struct OutputStat
{
  float Genotype_Likelihood[GENOTYPE_NUM];
  float Strand_Bias[ALPHABET_SIZE];
  float Base_Quality_Bias[ALPHABET_SIZE];
} OutputStat;

#pragma pack(1)
typedef struct MetaReference
{
  unsigned int amb;
  char refChar;
  char chrName[MAX_SEQ_NAME_LENGTH + 1];
  unsigned long long chrPos;
} MetaReference;
#pragma pack()

typedef InputStat MetaSnpCounter;

typedef struct GenotypeLikelihood
{

  double likelihood[GENOTYPE_NUM];
} GenotypeLikelihood;

typedef struct StrandBias
{
  float bias[ALPHABET_SIZE];
} StrandBias;

typedef struct BaseQualityBias
{
  float bias[ALPHABET_SIZE];
} BaseQualityBias;

typedef struct ReadPositionBias
{
  float lrRatio;
} DeltaStrandCount;

typedef struct AverageStrandCount
{
  float lrRatio;
} AverageStrandCount;

struct LeftRightFloatCount
{
  float left;
  float right;
};

struct LeftRightShortCount
{
  unsigned short left;
  unsigned short right;
};

struct LeftRightCharCount
{
  unsigned char left;
  unsigned char right;
};


typedef struct LeftRightFloatCount GCCount;

typedef struct LeftRightFloatCount AverageDepth;

typedef struct LeftRightFloatCount ReadQuality;

typedef struct LeftRightCharCount GCount;

typedef struct SnpCallingInfo
{
  double pD;
  int genotype;
  int secondGenotype;
  double flh[GENOTYPE_NUM];
  double bestD;
  double secondD;
  float rfProb;
  int rfMode;
} SnpCallingInfo;

typedef struct IncExcList
{
  unsigned int amb;
  char refChar;
  char altChar[5];
} IncExcList;

#pragma pack(1)
typedef struct MetaWindowInfo
{
  unsigned int weightedCount;
  unsigned int posStrandCount;
  unsigned int negStrandCount;
  unsigned char gcStrandCount;
  unsigned char gStrandCount;
  unsigned short polyrun;
  unsigned int indelHqCount;
  unsigned int indelLqCount;
  char isValid;
  char preBase;
} MetaWindowInfo;
#pragma pack()

typedef struct IndelInfo
{
  unsigned char *array;
  unsigned int size;
  unsigned int limit;
} IndelInfo;

__inline__ void checkIndelInfoSpace ( IndelInfo *indelInfo, unsigned int nextSize )
{
  if ( indelInfo->size + nextSize > indelInfo->limit )
    {
      indelInfo->limit *= 2;
      indelInfo->array = ( unsigned char * ) realloc ( indelInfo->array, indelInfo->limit * sizeof(unsigned char) );
    }
}

__inline__ void jumpIndelInfoIndex ( IndelInfo *indelInfo, unsigned int &indelInfoIndex )
{
  indelInfoIndex++;
  unsigned char indelLength = indelInfo->array[indelInfoIndex++];
  indelInfoIndex += ( indelLength + 3 ) / 4;
  indelInfoIndex += 4;
}

#endif
