#include <cstdio>
#include <cstring>

#include "CounterReader.h"

void addInsSeqs(InsertSeqs *insSeqs, int insertLength, char *insertSeqs, unsigned int hqCount, unsigned int lqCount) {
  int newSize = insSeqs->size + 1;
  insSeqs->insertSeq = (char **) realloc(insSeqs->insertSeq, newSize * sizeof(char *));
  insSeqs->hqCount = (unsigned int *) realloc(insSeqs->hqCount, newSize * sizeof(unsigned int));
  insSeqs->lqCount = (unsigned int *) realloc(insSeqs->lqCount, newSize * sizeof(unsigned int));

  insSeqs->insertSeq[insSeqs->size] = (char *) malloc((insertLength + 1) * sizeof(char));
  strcpy(insSeqs->insertSeq[insSeqs->size], insertSeqs);
  insSeqs->hqCount[insSeqs->size] = hqCount;
  insSeqs->lqCount[insSeqs->size] = lqCount;
  insSeqs->size++;
  return;
}

void addDelSeqs(DeleteSeqs *delSeqs, char deleteLength, unsigned int hqCount, unsigned int lqCount) {
  int newSize = delSeqs->size + 1;
  delSeqs->deleteLength = (char *) realloc(delSeqs->deleteLength, newSize * sizeof(char));
  delSeqs->hqCount = (unsigned int *) realloc(delSeqs->hqCount, newSize * sizeof(unsigned int));
  delSeqs->lqCount = (unsigned int *) realloc(delSeqs->lqCount, newSize * sizeof(unsigned int));

  delSeqs->deleteLength[delSeqs->size] = deleteLength;
  delSeqs->hqCount[delSeqs->size] = hqCount;
  delSeqs->lqCount[delSeqs->size] = lqCount;
  delSeqs->size++;
  return;
}

void freeSnpCounterInfo(SnpCounterInfo *info) {
  if (info == NULL) { return; }
  if (info->insSeqs != NULL) {
    for (int i = 0; i < info->insSeqs->size; ++i) { free(info->insSeqs->insertSeq[i]); }
    free(info->insSeqs->insertSeq);
    free(info->insSeqs->hqCount);
    free(info->insSeqs->lqCount);
    free(info->insSeqs);
  }
  if (info->delSeqs != NULL) {
    free(info->delSeqs->deleteLength);
    free(info->delSeqs->hqCount);
    free(info->delSeqs->lqCount);
    free(info->delSeqs);
  }
  free(info);
}

SnpCounterInfo *getCounterInfo(SnpBundle *snpBundle, unsigned int textLength, unsigned int pos, bool isIndel) {
  const char dnaMap[] = {'A', 'C', 'G', 'T'};
  SnpCounter *snpCounter = snpBundle->snpCounter;
  SnpOverflowCounterArray *snpOverflowCounterArray = snpBundle->snpOverflowCounterArray;
  unsigned int numCpuThreads = snpBundle->numOfCPUThreads;
  MemoryPool *pool = snpBundle->snpMemoryPool;

  SnpCounterInfo *counterInfo = (SnpCounterInfo *) malloc(sizeof(SnpCounterInfo));
  memset(counterInfo, 0, sizeof(SnpCounterInfo));

  counterInfo->pos = pos;

  unsigned short w1 = snpCounter[pos].weightedCount1;
  unsigned short w2 = snpCounter[pos].weightedCount2;
  if (((w1 & w2) & 0xC000) == 0xC000) {
    int region = findRegionByPosition(pos, textLength, numCpuThreads);
    unsigned int arrayIndex = snpCounter[pos].data.arrayIndex;
    memcpy(&(counterInfo->counter), &(snpOverflowCounterArray[region].counters[arrayIndex]),
           sizeof(SnpOverflowCounter));

    if (isIndel && (w1 & 0x2000)) {
      char insSeqBuf[1024];
      counterInfo->insSeqs = (InsertSeqs *) calloc(1, sizeof(InsertSeqs));
      if (w1 & 0x1000) {
        unsigned int arraySize = (w1 & 0xFFF);
        unsigned char *insertInfo = (unsigned char *) getAddress(pool,
                                                                 snpOverflowCounterArray[region].counters[arrayIndex].insertion.ptr);
        unsigned int infoIndex = 0;
        while (infoIndex < arraySize) {
          unsigned char insertLength = insertInfo[infoIndex++];
          if (insertLength == 0) {
            fprintf(stderr, "WARN: Zero insert length [%u]\n", pos);
            fprintf(stderr, "Region: %u, ArrayIndex: %u\n", region, arrayIndex);
          }

          for (int i = 0; i < insertLength; ++i) {
            char baseBit = (insertInfo[infoIndex + (i >> 2)] >> ((i & 3) << 1)) & 3;
            insSeqBuf[i] = dnaMap[baseBit];
          }
          insSeqBuf[insertLength] = 0;
          infoIndex += (insertLength + 3) / 4;
#ifndef AMPLICON
          unsigned short hqCount;
          hqCount = insertInfo[infoIndex++];
          hqCount |= insertInfo[infoIndex++] << 8;
          unsigned short lqCount;
          lqCount = insertInfo[infoIndex++];
          lqCount |= insertInfo[infoIndex++] << 8;
#else
          unsigned int hqCount;
          hqCount = insertInfo[infoIndex++];
          hqCount |= insertInfo[infoIndex++] << 8;
          hqCount |= insertInfo[infoIndex++] << 16;
          hqCount |= insertInfo[infoIndex++] << 24;
          unsigned int lqCount;
          lqCount = insertInfo[infoIndex++];
          lqCount |= insertInfo[infoIndex++] << 8;
          lqCount |= insertInfo[infoIndex++] << 16;
          lqCount |= insertInfo[infoIndex++] << 24;
#endif
          addInsSeqs(counterInfo->insSeqs, insertLength, insSeqBuf, hqCount, lqCount);
        }
      } else {
        unsigned char insertLength = ((w1 >> 8) & 0xF);
        unsigned char hqCount = (w1 & 0xFF);
        unsigned char lqCount = snpOverflowCounterArray[region].counters[arrayIndex].insertion.counters.lqCount;
        for (int i = 0; i < insertLength; ++i) {
          char baseBit = (snpOverflowCounterArray[region].counters[arrayIndex].insertion.counters.insertSeq[i >> 2]
              >> ((i & 3) << 1)) & 3;
          insSeqBuf[i] = dnaMap[baseBit];
        }
        insSeqBuf[insertLength] = 0;
        addInsSeqs(counterInfo->insSeqs, insertLength, insSeqBuf, hqCount, lqCount);
      }
    }
    if (isIndel && (w2 & 0x2000)) {
      counterInfo->delSeqs = (DeleteSeqs *) calloc(1, sizeof(DeleteSeqs));
      if (w2 & 0x1000) {
        unsigned int arraySize = (w2 & 0xFFF);
        unsigned char *deleteInfo = (unsigned char *) getAddress(pool,
                                                                 snpOverflowCounterArray[region].counters[arrayIndex].deletion.ptr);
        unsigned int infoIndex = 0;
        while (infoIndex < arraySize) {
          unsigned char deletionLength = deleteInfo[infoIndex++];
          if (deletionLength == 0) {
            fprintf(stderr, "WARN: Zero delete length [%u]\n", pos);
            fprintf(stderr, "Region: %u, ArrayIndex: %u\n", region, arrayIndex);
          }
#ifndef AMPLICON
          unsigned short hqCount;
          hqCount = deleteInfo[infoIndex++];
          hqCount |= deleteInfo[infoIndex++] << 8;
          unsigned short lqCount;
          lqCount = deleteInfo[infoIndex++];
          lqCount |= deleteInfo[infoIndex++] << 8;
#else
          unsigned int hqCount;
          hqCount = deleteInfo[infoIndex++];
          hqCount |= deleteInfo[infoIndex++] << 8;
          hqCount |= deleteInfo[infoIndex++] << 16;
          hqCount |= deleteInfo[infoIndex++] << 24;
          unsigned int lqCount;
          lqCount = deleteInfo[infoIndex++];
          lqCount |= deleteInfo[infoIndex++] << 8;
          lqCount |= deleteInfo[infoIndex++] << 16;
          lqCount |= deleteInfo[infoIndex++] << 24;
#endif
          addDelSeqs(counterInfo->delSeqs, deletionLength, hqCount, lqCount);
        }
      } else {
        unsigned char deletionLength = snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.deletionLength;
        unsigned char hqCount = snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.hqCount;
        unsigned char lqCount = snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.lqCount;
        addDelSeqs(counterInfo->delSeqs, deletionLength, hqCount, lqCount);
      }
    }
  } else {
    char baseBit1 = (w1 >> 14) & 3;
    counterInfo->counter.weightedCount[baseBit1] = w1 & 0x3FFF;
    counterInfo->counter.posStrandCount[baseBit1] = snpCounter[pos].data.strandCounts.posStrandCount1;
    counterInfo->counter.negStrandCount[baseBit1] = snpCounter[pos].data.strandCounts.negStrandCount1;

    char baseBit2 = (w2 >> 14) & 3;
    if (baseBit2 != baseBit1) {
      counterInfo->counter.weightedCount[baseBit2] = w2 & 0x3FFF;
      counterInfo->counter.posStrandCount[baseBit2] = snpCounter[pos].data.strandCounts.posStrandCount2;
      counterInfo->counter.negStrandCount[baseBit2] = snpCounter[pos].data.strandCounts.negStrandCount2;
    }
  }

  return counterInfo;
}

void printSnpCounterInfo(SnpCounterInfo *counter) {
  const char dnaMap[] = {'A', 'C', 'G', 'T'};
  printf("%u\n", counter->pos);
  for (int i = 0; i < ALPHABET_SIZE; ++i) {
    printf("%c %u %u %u\n", dnaMap[i], counter->counter.weightedCount[i], counter->counter.posStrandCount[i],
           counter->counter.negStrandCount[i]);
  }
  if (counter->insSeqs != NULL) {
    for (int j = 0; j < counter->insSeqs->size; ++j) {
      printf("I %s %u %u\n", counter->insSeqs->insertSeq[j], counter->insSeqs->hqCount[j],
             counter->insSeqs->lqCount[j]);
    }
  }
  if (counter->delSeqs != NULL) {
    for (int j = 0; j < counter->delSeqs->size; ++j) {
      printf("D %u %u %u\n", counter->delSeqs->deleteLength[j], counter->delSeqs->hqCount[j],
             counter->delSeqs->lqCount[j]);
    }
  }
}

int getCounterSupport(SnpBundle *snpBundle, unsigned int textLength, unsigned int pos) {
  SnpCounter *snpCounter = snpBundle->snpCounter;
  SnpOverflowCounterArray *snpOverflowCounterArray = snpBundle->snpOverflowCounterArray;
  unsigned int numCpuThreads = snpBundle->numOfCPUThreads;
  MemoryPool *pool = snpBundle->snpMemoryPool;

  unsigned int support = 0;

  unsigned short w1 = snpCounter[pos].weightedCount1;
  unsigned short w2 = snpCounter[pos].weightedCount2;
  if (((w1 & w2) & 0xC000) == 0xC000) {
    int region = findRegionByPosition(pos, textLength, numCpuThreads);
    unsigned int arrayIndex = snpCounter[pos].data.arrayIndex;
    for (int i = 0; i < ALPHABET_SIZE; ++i) {
      support += snpOverflowCounterArray[region].counters[arrayIndex].posStrandCount[i];
      support += snpOverflowCounterArray[region].counters[arrayIndex].negStrandCount[i];
    }
    if (w1 & 0x2000) {
      if (w1 & 0x1000) {
        unsigned int arraySize = (w1 & 0xFFF);
        unsigned char *insertInfo = (unsigned char *) getAddress(pool,
                                                                 snpOverflowCounterArray[region].counters[arrayIndex].insertion.ptr);
        unsigned int infoIndex = 0;
        while (infoIndex < arraySize) {
          unsigned char insertLength = insertInfo[infoIndex++];
          if (insertLength == 0) {
            fprintf(stderr, "WARN: Zero insert length [%u]\n", pos);
            fprintf(stderr, "Region: %u, ArrayIndex: %u\n", region, arrayIndex);
          }

          infoIndex += (insertLength + 3) / 4;
#ifndef AMPLICON
          unsigned short hqCount;
          hqCount = insertInfo[infoIndex++];
          hqCount |= insertInfo[infoIndex++] << 8;
          unsigned short lqCount;
          lqCount = insertInfo[infoIndex++];
          lqCount |= insertInfo[infoIndex++] << 8;
#else
          unsigned int hqCount;
          hqCount = insertInfo[infoIndex++];
          hqCount |= insertInfo[infoIndex++] << 8;
          hqCount |= insertInfo[infoIndex++] << 16;
          hqCount |= insertInfo[infoIndex++] << 24;
          unsigned int lqCount;
          lqCount = insertInfo[infoIndex++];
          lqCount |= insertInfo[infoIndex++] << 8;
          lqCount |= insertInfo[infoIndex++] << 16;
          lqCount |= insertInfo[infoIndex++] << 24;
#endif
          
          
        }
      } else {
        
        
      }
    }
    if (w2 & 0x2000) {
      if (w2 & 0x1000) {
        unsigned int arraySize = (w2 & 0xFFF);
        unsigned char *deleteInfo = (unsigned char *) getAddress(pool,
                                                                 snpOverflowCounterArray[region].counters[arrayIndex].deletion.ptr);
        unsigned int infoIndex = 0;
        while (infoIndex < arraySize) {
          unsigned char deletionLength = deleteInfo[infoIndex++];
          if (deletionLength == 0) {
            fprintf(stderr, "WARN: Zero delete length [%u]\n", pos);
            fprintf(stderr, "Region: %u, ArrayIndex: %u\n", region, arrayIndex);
          }
#ifndef AMPLICON
          unsigned short hqCount;
          hqCount = deleteInfo[infoIndex++];
          hqCount |= deleteInfo[infoIndex++] << 8;
          unsigned short lqCount;
          lqCount = deleteInfo[infoIndex++];
          lqCount |= deleteInfo[infoIndex++] << 8;
#else
          unsigned int hqCount;
          hqCount = deleteInfo[infoIndex++];
          hqCount |= deleteInfo[infoIndex++] << 8;
          hqCount |= deleteInfo[infoIndex++] << 16;
          hqCount |= deleteInfo[infoIndex++] << 24;
          unsigned int lqCount;
          lqCount = deleteInfo[infoIndex++];
          lqCount |= deleteInfo[infoIndex++] << 8;
          lqCount |= deleteInfo[infoIndex++] << 16;
          lqCount |= deleteInfo[infoIndex++] << 24;
#endif
          support += hqCount + lqCount;
        }
      } else {
        unsigned char hqCount = snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.hqCount;
        unsigned char lqCount = snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.lqCount;
        support += hqCount + lqCount;
      }
    }
  } else {
    char baseBit1 = (w1 >> 14) & 3;
    char baseBit2 = (w2 >> 14) & 3;
    support += snpCounter[pos].data.strandCounts.posStrandCount1;
    support += snpCounter[pos].data.strandCounts.negStrandCount1;

    if (baseBit2 != baseBit1) {
      support += snpCounter[pos].data.strandCounts.posStrandCount2;
      support += snpCounter[pos].data.strandCounts.negStrandCount2;
    }
  }

  return support;
}

