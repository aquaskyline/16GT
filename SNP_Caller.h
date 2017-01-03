#ifndef _SNP_CALLER_H_
#define _SNP_CALLER_H_

#include <stdio.h>
#include "SNP_Meta.h"

enum SNPMetaBufferStatus
{
  NOT_READY, READY, FINISHED
};

typedef struct SNPMetaBuffer
{
  MetaReference *reference;
  MetaSnpCounter *snpCounter;
  GenotypeLikelihood *genotypeLikelihood;
  StrandBias *strandBias;
  BaseQualityBias *baseQualityBias;
  DeltaStrandCount *deltaStrandCount;
  GCCount *gcCount;
  AverageStrandCount *avgStrandCount;
  AverageDepth *avgDepth;
  SnpCallingInfo *scInfo;
  MetaWindowInfo *mwInfo;
  IndelInfo *indelInfo;
  unsigned int batchSize;
  unsigned int startIdx;
  SNPMetaBufferStatus status;
} SNPMetaBuffer;

typedef struct SnpCallerWrapperObj
{
  SNPMetaBuffer *buffer0;
  SNPMetaBuffer *buffer1;
  FILE *snpOutput;
  unsigned int *attriSize;
} SnpCallerWrapperObj;

void filterSNP(MetaReference *reference, MetaSnpCounter *snpCounter,
               GenotypeLikelihood *genotypeLikelihood, SnpCallingInfo *scInfo,
               MetaWindowInfo *mwInfo, IndelInfo *indelInfo,
               unsigned int batchSize, unsigned int startIdx,
               FILE *snpOutput, unsigned int *attriSize);

SNPMetaBuffer *createSNPMetaBuffer(MetaReference *reference,
                                   MetaSnpCounter *snpCounter,
                                   GenotypeLikelihood *genotypeLikelihood,
                                   SnpCallingInfo *scInfo,
                                   MetaWindowInfo *mwInfo,
                                   IndelInfo *indelInfo);

void freeSNPMetaBuffer(SNPMetaBuffer *buffer);

void setReadySNPMetaBufferStatus(SNPMetaBuffer *buffer,
                                 unsigned int batchSize,
                                 unsigned int startIdx);

void setFinishSNPMetaBufferStatus(SNPMetaBuffer *buffer);

void waitFinishSNPMetaBuffer(SNPMetaBuffer *buffer);

SNPMetaBuffer *getIdleSNPMetaBuffer(SNPMetaBuffer *buffer0, SNPMetaBuffer *buffer1);

void startFilterSNPThread(SNPMetaBuffer *buffer0, SNPMetaBuffer *buffer1,
                          FILE *snpOutput,
                          unsigned int *attriSize,
                          pthread_t &filterSnpThread);

void calRandomForestpPedictProb(const char *snp_noRF_filename,
                                const char *snp_filename,
                                unsigned int attriSize, char verbose);

void closeFilterSNPThread(pthread_t &thread);

#endif
