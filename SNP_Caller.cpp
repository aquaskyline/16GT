#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "SNP_Caller.h"
#include "SNPFunctions.h"
#include "interpreter.h"

struct FourBasesDepth_s_ {
  short depth;
  short id;
};

static __inline__ void sort4_(FourBasesDepth_s_ *d) {
#define mind(x, y) ((x.depth<y.depth)?x:y)
#define maxd(x, y) ((x.depth<y.depth)?y:x)
#define SWAP(x, y) { FourBasesDepth_s_ tmp = maxd(d[x], d[y]); d[y] = mind(d[x], d[y]); d[x] = tmp; } 
  SWAP(0, 1);
  SWAP(2, 3);
  SWAP(0, 2);
  SWAP(1, 3);
  SWAP(1, 2);
#undef SWAP
#undef maxd
#undef mind
}

static const unsigned char offsetASiteLHSBBQB_[12] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44};

static __inline__ float offsetFloatSingle_(void *ptr, unsigned int a1) {
  return *(float *) (((char *) ptr) + (int) offsetASiteLHSBBQB_[a1]);
}

void resetBufferStatus(SNPMetaBuffer *buffer) {
  buffer->status = NOT_READY;
}

SNPMetaBuffer *createSNPMetaBuffer(MetaReference *reference,
                                   MetaSnpCounter *snpCounter,
                                   GenotypeLikelihood *genotypeLikelihood,
                                   SnpCallingInfo *scInfo,
                                   MetaWindowInfo *mwInfo,
                                   IndelInfo *indelInfo) {
  SNPMetaBuffer *buffer = (SNPMetaBuffer *) malloc(sizeof(SNPMetaBuffer));
  buffer->reference = reference;
  buffer->snpCounter = snpCounter;
  buffer->genotypeLikelihood = genotypeLikelihood;
  buffer->scInfo = scInfo;
  buffer->mwInfo = mwInfo;
  buffer->indelInfo = indelInfo;
  buffer->status = NOT_READY;
  return buffer;
}

void freeSNPMetaBuffer(SNPMetaBuffer *buffer) {
  free(buffer);
}

void setReadySNPMetaBufferStatus(SNPMetaBuffer *buffer,
                                 unsigned int batchSize,
                                 unsigned int startIdx) {
  buffer->batchSize = batchSize;
  buffer->startIdx = startIdx;
  buffer->status = READY;
}

void setFinishSNPMetaBufferStatus(SNPMetaBuffer *buffer) {
  buffer->status = FINISHED;
}

void waitFinishSNPMetaBuffer(SNPMetaBuffer *buffer) {
  while (1) {
    if (buffer->status == NOT_READY || buffer->status == FINISHED) {
      break;
    }
    sleep(1);
  }
}

SNPMetaBuffer *getIdleSNPMetaBuffer(SNPMetaBuffer *buffer0, SNPMetaBuffer *buffer1) {
  while (1) {
    if (buffer0->status == NOT_READY) {
      return buffer0;
    }

    if (buffer1->status == NOT_READY) {
      return buffer1;
    }

    if (buffer0->status == FINISHED || buffer1->status == FINISHED) {
      return NULL;
    }

    sleep(1);
  }
}

void outputPossibleSNPs(MetaReference *reference, unsigned char preBase,
                        MetaSnpCounter *snpCounter, SnpCallingInfo *scInfo,
                        IndelInfo *indelInfo, unsigned int &indelInfoIndex,
                        GenotypeLikelihood *genotypeLikelihood,
                        StrandBias *strandBias,
                        BaseQualityBias *baseQualityBias,
                        DeltaStrandCount *deltaStrandCount,
                        GCCount *gcCount,
                        GCount *gCount,
                        AverageStrandCount *avgStrandCount,
                        ReadQuality *rQuality,
                        unsigned short polyrun,
                        unsigned int indelHqCount, unsigned int indelLqCount,
                        unsigned char isIndel,
                        FILE *snpOutput) {
  unsigned char buffer[1024];
  unsigned short bufferSize = 0;
  memcpy(buffer + bufferSize, reference, sizeof(MetaReference));
  bufferSize += sizeof(MetaReference);
  memcpy(buffer + bufferSize, &preBase, sizeof(unsigned char));
  bufferSize += sizeof(unsigned char);
  memcpy(buffer + bufferSize, snpCounter, sizeof(MetaSnpCounter));
  bufferSize += sizeof(MetaSnpCounter);

#ifdef GENOTYPE_16
  if (isIndel) {
    char indelType = indelInfo->array[indelInfoIndex++];
    if (indelType != 'I' && indelType != 'D') {
      fprintf(stderr, "Error in Indel Type %d\n", indelType);
      exit(1);
    }
    unsigned char indelLength = indelInfo->array[indelInfoIndex++];
    memcpy(buffer + bufferSize, &indelType, sizeof(char));
    bufferSize += sizeof(char);
    memcpy(buffer + bufferSize, &indelLength, sizeof(unsigned char));
    bufferSize += sizeof(unsigned char);
    memcpy(buffer + bufferSize, &(indelInfo->array[indelInfoIndex]), sizeof(unsigned char) * (indelLength + 3) / 4);
    bufferSize += sizeof(unsigned char) * (indelLength + 3) / 4;
    indelInfoIndex += (indelLength + 3) / 4;
    memcpy(buffer + bufferSize, &(indelInfo->array[indelInfoIndex]), sizeof(unsigned short) * 2);
    bufferSize += sizeof(unsigned short) * 2;
    indelInfoIndex += 4;

    if (isIndel == 3) {
      char indelType = indelInfo->array[indelInfoIndex++];
      if (indelType != 'I' && indelType != 'D') {
        fprintf(stderr, "Error in Indel Type %d\n", indelType);
        exit(1);
      }
      unsigned char indelLength = indelInfo->array[indelInfoIndex++];
      memcpy(buffer + bufferSize, &indelType, sizeof(char));
      bufferSize += sizeof(char);
      memcpy(buffer + bufferSize, &indelLength, sizeof(unsigned char));
      bufferSize += sizeof(unsigned char);
      memcpy(buffer + bufferSize, &(indelInfo->array[indelInfoIndex]), sizeof(unsigned char) * (indelLength + 3) / 4);
      bufferSize += sizeof(unsigned char) * (indelLength + 3) / 4;
      indelInfoIndex += (indelLength + 3) / 4;
      memcpy(buffer + bufferSize, &(indelInfo->array[indelInfoIndex]), sizeof(unsigned short) * 2);
      bufferSize += sizeof(unsigned short) * 2;
      indelInfoIndex += 4;
    } else {
      char indelType = '*';
      memcpy(buffer + bufferSize, &indelType, sizeof(char));
      bufferSize += sizeof(char);
    }
  } else {
    char indelType = '*';
    memcpy(buffer + bufferSize, &indelType, sizeof(char));
    bufferSize += sizeof(char);
  }
#endif
  memcpy(buffer + bufferSize, scInfo, sizeof(SnpCallingInfo));
  bufferSize += sizeof(SnpCallingInfo);
  memcpy(buffer + bufferSize, genotypeLikelihood, sizeof(GenotypeLikelihood));
  bufferSize += sizeof(GenotypeLikelihood);
  if (isIndel) {
    memcpy(buffer + bufferSize, strandBias, sizeof(StrandBias) * 3);
    bufferSize += sizeof(StrandBias) * 3;
    memcpy(buffer + bufferSize, baseQualityBias, sizeof(BaseQualityBias) * 3);
    bufferSize += sizeof(BaseQualityBias) * 3;
  } else {
    memcpy(buffer + bufferSize, strandBias, sizeof(StrandBias));
    bufferSize += sizeof(StrandBias);
    memcpy(buffer + bufferSize, baseQualityBias, sizeof(BaseQualityBias));
    bufferSize += sizeof(BaseQualityBias);
  }
  memcpy(buffer + bufferSize, deltaStrandCount, sizeof(DeltaStrandCount));
  bufferSize += sizeof(DeltaStrandCount);
  memcpy(buffer + bufferSize, gcCount, sizeof(GCCount));
  bufferSize += sizeof(GCCount);
  memcpy(buffer + bufferSize, gCount, sizeof(GCount));
  bufferSize += sizeof(GCount);
  memcpy(buffer + bufferSize, avgStrandCount, sizeof(AverageStrandCount));
  bufferSize += sizeof(AverageStrandCount);
  memcpy(buffer + bufferSize, rQuality, sizeof(ReadQuality));
  bufferSize += sizeof(ReadQuality);
  memcpy(buffer + bufferSize, &polyrun, sizeof(unsigned short));
  bufferSize += sizeof(unsigned short);
  memcpy(buffer + bufferSize, &indelHqCount, sizeof(unsigned int));
  bufferSize += sizeof(unsigned int);
  memcpy(buffer + bufferSize, &indelLqCount, sizeof(unsigned int));
  bufferSize += sizeof(unsigned int);

  if (bufferSize > 1024) {
    fprintf(stderr, "Buffer overflow\n");
  }

  fwrite(buffer, sizeof(unsigned char), bufferSize, snpOutput);
}

void computeMeta(unsigned int index, unsigned int batchSize, MetaReference *reference,
                 MetaSnpCounter *snpCounter, SnpCallingInfo *scInfo,
                 MetaWindowInfo *mwInfo,
                 IndelInfo *indelInfo, unsigned int &indelInfoIndex,
                 GenotypeLikelihood *genotypeLikelihood,
                 StrandBias *strandBias,
                 BaseQualityBias *baseQualityBias,
                 DeltaStrandCount *deltaStrandCount,
                 GCCount *gcCount, GCount *gCount,
                 AverageStrandCount *avgStrandCount,
                 ReadQuality *rQuality,
                 FILE *snpOutput) {
  char isIndel = 0;
  unsigned char nextBp = 0;
#ifdef GENOTYPE_16
  
  
  if (scInfo[index].genotype > 9 && mwInfo[index + SNP_META_WINDOW_SIZE].isValid >= 2) {
    char indelType = indelInfo->array[indelInfoIndex];
    unsigned char indelLength = indelInfo->array[indelInfoIndex + 1];
    if (indelType == 'D') {
      nextBp = indelLength;
    }
    isIndel = mwInfo[index + SNP_META_WINDOW_SIZE].isValid;
  } else if (scInfo[index].genotype > 9) {
    if (scInfo[index].genotype != UNIDENTIFIED_GENOTYPE) {
      fprintf(stderr, "Error in selecting Genotype : %u\n", reference[index].amb);
      exit(1);
    }
  } else if (mwInfo[index + SNP_META_WINDOW_SIZE].isValid == 2) {
    jumpIndelInfoIndex(indelInfo, indelInfoIndex);
  } else if (mwInfo[index + SNP_META_WINDOW_SIZE].isValid == 3) {
    jumpIndelInfoIndex(indelInfo, indelInfoIndex);
    jumpIndelInfoIndex(indelInfo, indelInfoIndex);
  }
#endif

  if (isIndel) {
    if (index != 0) {
      computeStrandBias(&(snpCounter[index - 1]), &(strandBias[0]));
    } else {
      memset(&(strandBias[0]), 0, sizeof(StrandBias));
    }
    if (index + nextBp < batchSize) {
      computeStrandBias(&(snpCounter[index + nextBp]), &(strandBias[1]));
    } else {
      memset(&(strandBias[1]), 0, sizeof(StrandBias));
    }
    computeStrandBias(&(snpCounter[index]), &(strandBias[2]));
    if (index != 0) {
      computeBaseQualityBias(&(snpCounter[index - 1]), &(baseQualityBias[0]));
    } else {
      memset(&(baseQualityBias[0]), 0, sizeof(BaseQualityBias));
    }
    if (index + nextBp < batchSize) {
      computeBaseQualityBias(&(snpCounter[index + nextBp]), &(baseQualityBias[1]));
    } else {
      memset(&(baseQualityBias[1]), 0, sizeof(BaseQualityBias));
    }
    computeBaseQualityBias(&(snpCounter[index]), &(baseQualityBias[2]));
  } else {
    computeStrandBias(&(snpCounter[index]), &(strandBias[0]));
    computeBaseQualityBias(&(snpCounter[index]), &(baseQualityBias[0]));
  }

  int strandSum[2] = {0};
  int deltaSum[2] = {0};
  int depthSum[2] = {0};
  unsigned char gcSum[2] = {0};
  unsigned char gSum[2] = {0};
  unsigned char validSum[2] = {0};
  float avgSC[2] = {0.0f};
  float avgD[2] = {0.0f};

  unsigned int i;
  for (i = index - SNP_META_WINDOW_SIZE; i < index; ++i) {
    if (mwInfo[i].isValid) {
      unsigned int strandCount0 = mwInfo[i].posStrandCount + mwInfo[i].negStrandCount;
      strandSum[0] += strandCount0;
      depthSum[0] += mwInfo[i].weightedCount;
      gcSum[0] += mwInfo[i].gcStrandCount;
      validSum[0] += 1;

      if (mwInfo[i + 1].isValid) {
        unsigned int strandCount1 = mwInfo[i + 1].posStrandCount + mwInfo[i + 1].negStrandCount;
        deltaSum[0] += (strandCount0 > strandCount1 ? strandCount0 - strandCount1 : strandCount1 - strandCount0);
      }
    }
  }

  for (i = index + nextBp; i < index + nextBp + SNP_META_WINDOW_SIZE; ++i) {
    if (mwInfo[i].isValid) {
      unsigned int strandCount0 = mwInfo[i].posStrandCount + mwInfo[i].negStrandCount;
      strandSum[1] += strandCount0;
      depthSum[1] += mwInfo[i].weightedCount;
      gcSum[1] += mwInfo[i].gcStrandCount;
      validSum[1] += 1;

      if (mwInfo[i + 1].isValid) {
        unsigned int strandCount1 = mwInfo[i + 1].posStrandCount + mwInfo[i + 1].negStrandCount;
        deltaSum[1] += (strandCount0 > strandCount1 ? strandCount0 - strandCount1 : strandCount1 - strandCount0);
      }
    }
  }

  if (deltaSum[0] | deltaSum[1]) {
    deltaStrandCount->lrRatio = (deltaSum[0] - deltaSum[1] + 0.0f) / (deltaSum[0] + deltaSum[1]);
  } else {
    deltaStrandCount->lrRatio = 0.0f;
  }

  if (validSum[0]) {
    gcCount->left = (gcSum[0] + 0.0f) / validSum[0];
    avgSC[0] = (strandSum[0] + 0.0f) / validSum[0];
    avgD[0] = (depthSum[0] + 0.0f) / validSum[0];
  } else {
    gcCount->left = 0.0f;
  }

  if (validSum[1]) {
    gcCount->right = (gcSum[1] + 0.0f) / validSum[1];
    avgSC[1] = (strandSum[1] + 0.0f) / validSum[1];
    avgD[1] = (depthSum[1] + 0.0f) / validSum[1];
  } else {
    gcCount->right = 0.0f;
  }

  if (fabs(avgSC[0] + avgSC[1]) > FLT_EPSILON) {
    avgStrandCount->lrRatio = (avgSC[0] - avgSC[1]) / (avgSC[0] + avgSC[1]);
    if (fabs(avgSC[0]) > FLT_EPSILON) {
      rQuality->left = avgD[0] / avgSC[0];
    } else {
      rQuality->left = 0.0f;
    }

    if (fabs(avgSC[1]) > FLT_EPSILON) {
      rQuality->right = avgD[1] / avgSC[1];
    } else {
      rQuality->right = 0.0f;
    }
  } else {
    avgStrandCount->lrRatio = 0.0f;
  }

  for (i = index - 3; i < index; ++i) {
    if (mwInfo[i].isValid) {
      gSum[0] += mwInfo[i].gStrandCount;
    }
  }

  for (i = index + nextBp; i < index + nextBp + 3; ++i) {
    if (mwInfo[i].isValid) {
      gSum[1] += mwInfo[i].gStrandCount;
    }
  }
  gCount->left = gSum[0];
  gCount->right = gSum[1];

  FourBasesDepth_s_ v[4];

  v[0].depth = snpCounter[index].W[0];
  v[0].id = 0;
  v[1].depth = snpCounter[index].W[1];
  v[1].id = 1;
  v[2].depth = snpCounter[index].W[2];
  v[2].id = 2;
  v[3].depth = snpCounter[index].W[3];
  v[3].id = 3;
  sort4_(v);

  outputPossibleSNPs(&(reference[index]), mwInfo[index + SNP_META_WINDOW_SIZE].preBase,
                     &(snpCounter[index]),
                     &(scInfo[index]),
                     indelInfo, indelInfoIndex,
                     &(genotypeLikelihood[index]),
                     strandBias, baseQualityBias, deltaStrandCount,
                     gcCount, gCount, avgStrandCount,
                     rQuality, mwInfo[index + SNP_META_WINDOW_SIZE].polyrun,
                     mwInfo[index + SNP_META_WINDOW_SIZE].indelHqCount,
                     mwInfo[index + SNP_META_WINDOW_SIZE].indelLqCount,
                     isIndel,
                     snpOutput);
}

void filterSNP(MetaReference *reference, MetaSnpCounter *snpCounter,
               GenotypeLikelihood *genotypeLikelihood, SnpCallingInfo *scInfo,
               MetaWindowInfo *mwInfo, IndelInfo *indelInfo,
               unsigned int batchSize, unsigned int startIdx,
               FILE *snpOutput, unsigned int *attriSize) {
  unsigned int filterTable[ALPHABET_SIZE] = {0, 4, 7, 9};

  StrandBias strandBias[3];
  BaseQualityBias baseQualityBias[3];
  DeltaStrandCount deltaStrandCount;
  GCCount gcCount;
  GCount gCount;
  AverageStrandCount avgStrandCount;
  ReadQuality rQuality;

  unsigned int indelInfoIndex = 0;
  unsigned int i;
  unsigned int bound = startIdx + batchSize;
  unsigned int mwIdx;
  for (i = startIdx; i < bound; ++i) {
    mwIdx = i + SNP_META_WINDOW_SIZE;
    if (mwInfo[mwIdx].isValid &&
        !(filterTable[reference[i].refChar] == scInfo[i].genotype || scInfo[i].genotype >= UNIDENTIFIED_GENOTYPE)) {
      (*attriSize)++;
      computeMeta(i, batchSize, reference, snpCounter,
                  scInfo, mwInfo,
                  indelInfo, indelInfoIndex,
                  genotypeLikelihood,
                  strandBias, baseQualityBias, &deltaStrandCount,
                  &gcCount, &gCount, &avgStrandCount,
                  &rQuality,
                  snpOutput);
    } else if (mwInfo[mwIdx].isValid == 2) {
      jumpIndelInfoIndex(indelInfo, indelInfoIndex);
    } else if (mwInfo[mwIdx].isValid == 3) {
      jumpIndelInfoIndex(indelInfo, indelInfoIndex);
      jumpIndelInfoIndex(indelInfo, indelInfoIndex);
    }
    if (indelInfoIndex > indelInfo->size) {
      fprintf(stderr, "Error in reading IndelInfo : size = %u, [%u, %u]\n", indelInfo->size, i, bound);
      exit(1);
    }
  }
}

void filterSNPPipeline(SNPMetaBuffer *buffer0, SNPMetaBuffer *buffer1,
                       FILE *snpOutput,
                       unsigned int *attriSize) {
  if (!buffer0 || !buffer1) {
    fprintf(stderr, "Null Pointer Exception : SNPMetaBuffer\n");
    exit(1);
  }

  SNPMetaBuffer *buffer;

  while (1) {
    while (1) {
      if (buffer0->status == NOT_READY && buffer1->status == NOT_READY) {
        sleep(1);
      } else if (buffer0->status == READY && buffer1->status == READY) {
        
        
        if (buffer0->reference->amb < buffer1->reference->amb) {
          buffer = buffer0;
        } else {
          buffer = buffer1;
        }
        break;
      } else if (buffer0->status == READY) {
        buffer = buffer0;
        break;
      } else if (buffer1->status == READY) {
        buffer = buffer1;
        break;
      } else 
      {
        return;
      }
    }

    filterSNP(buffer->reference, buffer->snpCounter,
              buffer->genotypeLikelihood, buffer->scInfo,
              buffer->mwInfo, buffer->indelInfo,
              buffer->batchSize, buffer->startIdx,
              snpOutput, attriSize);

    buffer->status = NOT_READY;
  }
}

void *filterSNPWrapper(void *ptr) {
  SnpCallerWrapperObj *obj = (SnpCallerWrapperObj *) ptr;

  filterSNPPipeline(obj->buffer0, obj->buffer1, obj->snpOutput, obj->attriSize);

  free(obj);
  return NULL;
}

void startFilterSNPThread(SNPMetaBuffer *buffer0, SNPMetaBuffer *buffer1,
                          FILE *snpOutput,
                          unsigned int *attriSize,
                          pthread_t &thread) {
  SnpCallerWrapperObj *obj = (SnpCallerWrapperObj *) malloc(sizeof(SnpCallerWrapperObj));
  obj->buffer0 = buffer0;
  obj->buffer1 = buffer1;
  obj->snpOutput = snpOutput;
  obj->attriSize = attriSize;

  int t;
  t = pthread_create(&thread, NULL, filterSNPWrapper, obj);

  if (t) {
    fprintf(stderr, "Error in creating SNP Calling Thread.\n");
  }
}

void calRandomForestpPedictProb(const char *snp_noRF_filename,
    const char *snp_filename,
    unsigned int attriSize, char verbose) {

  FILE *snpOutput = fopen(snp_noRF_filename, "rb");
  FILE *newSnpOutput = fopen(snp_filename, "w");

  printf("Called Entries : %u\n", attriSize);
  if (snpOutput == NULL) {
    fprintf(stderr, "[SNP] %s does not exist\n", snp_noRF_filename);
    exit(1);
  }
  if (newSnpOutput == NULL) {
    fprintf(stderr, "[SNP] %s does not exist\n", snp_filename);
    exit(1);
  }

  float rfProb = 0;
  int rfMode = 0;

  MetaReference reference;
  unsigned char preBase;
  MetaSnpCounter snpCounter;
  char indelType[2];
  unsigned char indelLength[2];
  unsigned char indelPattern[2][MAX_PATTERN_LENGTH];
  unsigned short hqCount[2], lqCount[2];

  unsigned char isIndel = 0;

  if (verbose) { printf("Outputing Results...\n"); }
  printHeader(newSnpOutput);

  for (unsigned int i = 0; i < attriSize; i++) {

    fread(&reference, sizeof(MetaReference), 1, snpOutput);
    fread(&preBase, sizeof(unsigned char), 1, snpOutput);
    fread(&snpCounter, sizeof(MetaSnpCounter), 1, snpOutput);

    isIndel = 0;
    indelType[0] = '*';
    indelType[1] = '*';
    for (unsigned int j = 0; j < 2; j++) {
      fread(&(indelType[j]), sizeof(char), 1, snpOutput);

      if (indelType[j] == '*') {
        break;
      } else {
        isIndel = 1;
        fread(&(indelLength[j]), sizeof(unsigned char), 1, snpOutput);
        fread(indelPattern[j], sizeof(unsigned char), (indelLength[j] + 3) / 4, snpOutput);
        fread(&(hqCount[j]), sizeof(unsigned short), 1, snpOutput);
        fread(&(lqCount[j]), sizeof(unsigned short), 1, snpOutput);
      }
    }

    SnpCallingInfo scInfo;
    GenotypeLikelihood genotypeLikelihood;
    StrandBias strandBias[3];
    BaseQualityBias baseQualityBias[3];
    DeltaStrandCount deltaStrandCount;
    GCCount gcCount;
    GCount gCount;
    AverageStrandCount avgStrandCount;
    ReadQuality rQuality;
    unsigned short polyrun;
    unsigned int indelHqCount;
    unsigned int indelLqCount;

    fread(&scInfo, sizeof(SnpCallingInfo), 1, snpOutput);

    scInfo.rfProb = rfProb;
    scInfo.rfMode = rfMode;

    fread(&genotypeLikelihood, sizeof(GenotypeLikelihood), 1, snpOutput);
    if (isIndel) {
      fread(strandBias, sizeof(StrandBias), 3, snpOutput);
      fread(baseQualityBias, sizeof(BaseQualityBias), 3, snpOutput);
    } else {
      fread(strandBias, sizeof(StrandBias), 1, snpOutput);
      fread(baseQualityBias, sizeof(BaseQualityBias), 1, snpOutput);
    }
    fread(&deltaStrandCount, sizeof(DeltaStrandCount), 1, snpOutput);
    fread(&gcCount, sizeof(GCCount), 1, snpOutput);
    fread(&gCount, sizeof(GCount), 1, snpOutput);
    fread(&avgStrandCount, sizeof(AverageStrandCount), 1, snpOutput);
    fread(&rQuality, sizeof(ReadQuality), 1, snpOutput);
    fread(&polyrun, sizeof(unsigned short), 1, snpOutput);
    fread(&indelHqCount, sizeof(unsigned int), 1, snpOutput);
    fread(&indelLqCount, sizeof(unsigned int), 1, snpOutput);

    printPossibleSNP(&reference, &snpCounter, &scInfo, isIndel, indelType[0], indelLength[0], indelPattern[0],
        hqCount[0], lqCount[0], indelType[1], indelLength[1], indelPattern[1], hqCount[1], lqCount[1],
        &genotypeLikelihood, strandBias, baseQualityBias,
        &deltaStrandCount, &gcCount, &gCount, &avgStrandCount, &rQuality, polyrun, indelHqCount,
        indelLqCount, newSnpOutput);
  }

  fclose(snpOutput);
  fclose(newSnpOutput);
}

void closeFilterSNPThread(pthread_t &filterSnpThread) {
  if (pthread_join(filterSnpThread, NULL)) {
    fprintf(stderr, "[SNP] Error in joining SNP Calling Thread.\n");
  }
}


