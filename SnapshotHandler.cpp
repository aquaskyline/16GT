#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "SnapshotHandler.h"
#include "SNP_Meta.h"

size_t writeSnpCounter(SnpCounter *snpCounter, unsigned int textLength,
                       ExomeRegion *region, unsigned int regionSize,
                       FILE *outputSnpCounter,
                       char mode) {
  fwrite(&mode, sizeof(char), 1, outputSnpCounter);
  fwrite(&textLength, sizeof(unsigned int), 1, outputSnpCounter);
  if (mode & 1) {
    fwrite(snpCounter, sizeof(SnpCounter), textLength, outputSnpCounter);
    return sizeof(char) + sizeof(unsigned int) + textLength * sizeof(SnpCounter);
  } else if (mode & 2) {
    unsigned int *bitvector = (unsigned int *) malloc(((textLength + 31) >> 5) * sizeof(unsigned int));
    memset(bitvector, 0, ((textLength + 31) >> 5) * sizeof(unsigned int));
    for (unsigned int i = 0; i < regionSize; ++i) {
      unsigned int start = region[i].startPos > SNP_META_WINDOW_SIZE ? region[i].startPos - SNP_META_WINDOW_SIZE : 0;
      unsigned int end =
          region[i].endPos < textLength - SNP_META_WINDOW_SIZE ? region[i].endPos + SNP_META_WINDOW_SIZE : textLength -
                                                                                                           1;
      for (unsigned int j = start; j <= end; ++j) {
        bitvector[j >> 5] |= (1 << (31 - (j & 31)));
      }
      
    }
    unsigned int dumpCnt = 0;
    for (unsigned int i = 0; i < textLength; ++i) {
      if (bitvector[i >> 5] & (1 << (31 - (i & 31)))) {
        fwrite(&i, sizeof(unsigned int), 1, outputSnpCounter);
        fwrite(&(snpCounter[i]), sizeof(SnpCounter), 1, outputSnpCounter);
        ++dumpCnt;
      }
    }
    free(bitvector);
    return sizeof(char) + sizeof(unsigned int) + dumpCnt * (sizeof(unsigned int) + sizeof(SnpCounter));
  }
}

size_t writeDscCounter(SnpDirectionalSoftclipCounter *dscCounter, unsigned int dscSize,
                       ExomeRegion *region, unsigned int regionSize,
                       char isExome,
                       FILE *dscResult) {
  unsigned int dscCnt = dscSize;
  if (isExome) {
    unsigned int idx = 0;
    dscCnt = 0;
    for (unsigned int e = 0; e < regionSize; ++e) {
      while (idx < dscSize && dscCounter[idx].ambPosition < region[e].startPos) { ++idx; }
      while (idx < dscSize && region[e].startPos <= dscCounter[idx].ambPosition &&
             dscCounter[idx].ambPosition <= region[e].endPos) {
        fwrite(&(dscCounter[idx]), sizeof(SnpDirectionalSoftclipCounter), 1, dscResult);
        ++dscCnt;
        ++idx;
      }
    }
    return dscCnt * sizeof(SnpDirectionalSoftclipCounter);
  } else {
    fwrite(dscCounter, sizeof(SnpDirectionalSoftclipCounter), dscSize, dscResult);
    return dscCnt * sizeof(SnpDirectionalSoftclipCounter);
  }
}

void writeHeader(FILE *file, SnpBundle *snpStatBundle,
                 ExomeRegion *region, unsigned int regionSize,
                 char isExome) {
  unsigned int textLength = snpStatBundle->textLength;
  SnpOverflowCounterArray *snpOverflowCounterArray = snpStatBundle->snpOverflowCounterArray;
  int numOfCPUThreads = snpStatBundle->numOfCPUThreads;
  unsigned int dscSize = snpStatBundle->snpDirectionalSoftclipCounterSize;
  MemoryPool *mp = snpStatBundle->snpMemoryPool;
  SnpDirectionalSoftclipCounter *dscCounter = snpStatBundle->snpDirectionalSoftclipCounter;

  int pageSize = sysconf(_SC_PAGE_SIZE);
  if (pageSize < 1024) {
    fprintf(stderr, "System Page Size is too small\n");
    exit(1);
  }
  fwrite(&pageSize, sizeof(int), 1, file);
  unsigned int headerSize = pageSize;
  size_t subfile_start;
  size_t subfile_length;
  size_t curfile_end = headerSize;

#define writeSubFileHeader(start, length, File, pageSize) \
  fwrite ( &start, sizeof (size_t), 1, File); \
  fwrite ( &length, sizeof (size_t), 1, File); \
  (length) = ((length) + (pageSize) - 1 ) / (pageSize) * (pageSize); \

  
  subfile_start = curfile_end;
  subfile_length = sizeof(unsigned int) + sizeof(size_t) + sizeof(char) * mp->capacity;
  writeSubFileHeader (subfile_start, subfile_length, file, pageSize);
  curfile_end += subfile_length;

  
  subfile_start = curfile_end;
  subfile_length = sizeof(char) + sizeof(unsigned int);
  if (isExome) {
    unsigned int *bitvector = (unsigned int *) malloc(((textLength + 31) >> 5) * sizeof(unsigned int));
    memset(bitvector, 0, ((textLength + 31) >> 5) * sizeof(unsigned int));
    for (unsigned int i = 0; i < regionSize; ++i) {
      unsigned int start = region[i].startPos > SNP_META_WINDOW_SIZE ? region[i].startPos - SNP_META_WINDOW_SIZE : 0;
      unsigned int end =
          region[i].endPos < textLength - SNP_META_WINDOW_SIZE ? region[i].endPos + SNP_META_WINDOW_SIZE : textLength -
                                                                                                           1;
      for (unsigned int j = start; j <= end; ++j) {
        bitvector[j >> 5] |= (1 << (31 - (j & 31)));
      }
    }
    unsigned int dumpCnt = 0;
    for (unsigned int i = 0; i < textLength; ++i) {
      if (bitvector[i >> 5] & (1 << (31 - (i & 31)))) {
        ++dumpCnt;
      }
    }
    free(bitvector);
    subfile_length += dumpCnt * (sizeof(unsigned int) + sizeof(SnpCounter));
  } else {
    subfile_length += sizeof(SnpCounter) * snpStatBundle->textLength;
  }
  writeSubFileHeader (subfile_start, subfile_length, file, pageSize);
  curfile_end += subfile_length;

  
  subfile_start = curfile_end;
  subfile_length = sizeof(unsigned int);
  for (unsigned int s = 0; s < numOfCPUThreads; ++s) {
    unsigned int overflowSize = snpOverflowCounterArray[s].size;
    subfile_length += sizeof(unsigned int) + sizeof(unsigned int) + sizeof(SnpOverflowCounter) * overflowSize;
  }
  writeSubFileHeader (subfile_start, subfile_length, file, pageSize);
  curfile_end += subfile_length;

  
  subfile_start = curfile_end;
  subfile_length = isExome ? 0 : ((textLength + 31) >> 5) * sizeof(unsigned int);
  writeSubFileHeader (subfile_start, subfile_length, file, pageSize);
  curfile_end += subfile_length;

  
  subfile_start = curfile_end;
  subfile_length = 0;
  unsigned int dscCnt = dscSize;
  if (isExome) {
    dscCnt = 0;
    unsigned int idx = 0;
    for (unsigned int e = 0; e < regionSize; ++e) {
      while (idx < dscSize && dscCounter[idx].ambPosition < region[e].startPos) { ++idx; }
      while (idx < dscSize && region[e].startPos <= dscCounter[idx].ambPosition &&
             dscCounter[idx].ambPosition <= region[e].endPos) {
        ++dscCnt;
        ++idx;
      }
    }
  }
  subfile_length = sizeof(SnpDirectionalSoftclipCounter) * dscCnt;
  writeSubFileHeader (subfile_start, subfile_length, file, pageSize);
  curfile_end += subfile_length;

  char amplicon = 0;
  amplicon = 0;

  fwrite(&amplicon, sizeof(char), 1, file);

  
  int paddingLength = (pageSize - sizeof(int) - sizeof(size_t) * 10 - sizeof(char));
  char pad = 0;
  for (int i = 0; i < paddingLength; ++i) { fwrite(&pad, sizeof(char), 1, file); }

  return;
}

void writeSnpInfoSnapshot(SnpBundle *snpStatBundle,
                          ExomeRegion *region, unsigned int regionSize,
                          unsigned char isExome,
                          const char *snapshotFilename) {
  unsigned int textLength = snpStatBundle->textLength;
  SnpCounter *snpCounter = snpStatBundle->snpCounter;
  SnpOverflowCounterArray *snpOverflowCounterArray = snpStatBundle->snpOverflowCounterArray;
  int numOfCPUThreads = snpStatBundle->numOfCPUThreads;
  unsigned int *invalidSnpCounterPos = snpStatBundle->invalidSnpCounterPos;
  MemoryPool *pool = snpStatBundle->snpMemoryPool;
  int pageSize = sysconf(_SC_PAGE_SIZE);
  char padding = 0;
  int paddingLength = 0;

#define padding(size) \
  paddingLength = pageSize - ( ( size ) % pageSize ) ; \
  for ( int i=0; i< paddingLength; ++i ) \
  { fwrite ( &padding, sizeof ( char ), 1, snapshot ); }

  
  FILE *snapshot = fopen(snapshotFilename, "wb");
  writeHeader(snapshot, snpStatBundle, region, regionSize, isExome);

  
  fwrite(&(pool->curPtr), sizeof(unsigned int), 1, snapshot);
  fwrite(&(pool->capacity), sizeof(size_t), 1, snapshot);
  fwrite(pool->address, sizeof(char), pool->capacity, snapshot);
  padding (sizeof(unsigned int) + sizeof(size_t) + pool->capacity * sizeof(char));

  
  size_t counterSize = writeSnpCounter(snpCounter, textLength, region, regionSize, snapshot, isExome ? 2 : 1);
  padding (counterSize);

  
  counterSize = sizeof(unsigned int);
  fwrite(&(numOfCPUThreads), sizeof(unsigned int), 1, snapshot);
  for (unsigned int s = 0; s < numOfCPUThreads; ++s) {
    unsigned int overflowSize = snpOverflowCounterArray[s].size;
    fwrite(&(overflowSize), sizeof(unsigned int), 1, snapshot);
    fwrite(&(snpOverflowCounterArray[s].limit), sizeof(unsigned int), 1, snapshot);
    fwrite(snpOverflowCounterArray[s].counters, sizeof(SnpOverflowCounter), overflowSize, snapshot);
    counterSize += sizeof(unsigned int) + sizeof(unsigned int) + sizeof(SnpOverflowCounter) * overflowSize;
  }
  padding (counterSize);

  
  if (!isExome) {
    fwrite(invalidSnpCounterPos, sizeof(unsigned int), (textLength + 31) >> 5, snapshot);
    padding (sizeof(unsigned int) * ((textLength + 31) >> 5));
  }

  
  SnpDirectionalSoftclipCounter *dscCounter = snpStatBundle->snpDirectionalSoftclipCounter;
  unsigned int dscSize = snpStatBundle->snpDirectionalSoftclipCounterSize;
  counterSize = writeDscCounter(dscCounter, dscSize, region, regionSize, isExome, snapshot);
  padding (counterSize);

  fclose(snapshot);
}

void readHeader(FILE *file, int &pageSize,
                size_t *subfile_start, size_t *subfile_length) {
  int sysPageSize = sysconf(_SC_PAGE_SIZE);
  if (sysPageSize < 1024) {
    fprintf(stderr, "System Page Size is too small\n");
    exit(1);
  }
  fread(&pageSize, sizeof(int), 1, file);
  if (sysPageSize != pageSize) {
    fprintf(stderr, "System paga size is not matched to the snapshot\n");
    exit(1);
  }
  for (int i = 0; i < 5; ++i) {
    fread(&subfile_start[i], sizeof(size_t), 1, file);
    fread(&subfile_length[i], sizeof(size_t), 1, file);
    
  }

  return;
}

void readSnpInfoSnapshot(SnpBundle *snpStatBundle, unsigned int &textLength,
                         const char *snapshotFilename) {
  size_t subfile_start[5], subfile_length[5];
  int pageSize;

  FILE *snapshot = fopen(snapshotFilename, "rb");
  if (!snapshot) {
    fprintf(stderr, "Error opening %s\n", snapshotFilename);
    exit(1);
  }
  readHeader(snapshot, pageSize, subfile_start, subfile_length);
  fseek(snapshot, subfile_start[1], SEEK_SET);
  char mode = 0;
  fread(&mode, sizeof(char), 1, snapshot);
  fread(&textLength, sizeof(unsigned int), 1, snapshot);
  snpStatBundle->textLength = textLength;
  snpStatBundle->snpCounter = (SnpCounter *) malloc(textLength * sizeof(SnpCounter));
  if (snpStatBundle->snpCounter == NULL) {
    fprintf(stderr, "Memory Allocation Failed!\n\n");
    exit(1);
  }
  memset(snpStatBundle->snpCounter, 0, textLength * sizeof(SnpCounter));
  if (mode & 1) {
    fread(snpStatBundle->snpCounter, sizeof(SnpCounter), textLength, snapshot);
  } else if (mode & 2) {
    int bRead = subfile_length[1] / (sizeof(unsigned int) + sizeof(SnpCounter));
    unsigned int pos;
    while (bRead) {
      fread(&pos, sizeof(unsigned int), 1, snapshot);
      if (bRead <= 0) {
        break;
      }
      fread(&(snpStatBundle->snpCounter[pos]), sizeof(SnpCounter), 1, snapshot);
      --bRead;
    }
  }

  fseek(snapshot, subfile_start[2], SEEK_SET);
  fread(&(snpStatBundle->numOfCPUThreads), sizeof(unsigned int), 1, snapshot);

  snpStatBundle->snpOverflowCounterArray = (SnpOverflowCounterArray *) malloc(
      snpStatBundle->numOfCPUThreads * sizeof(SnpOverflowCounterArray));
  for (unsigned int i = 0; i < snpStatBundle->numOfCPUThreads; ++i) {
    unsigned int overflowSize, overflowLimit;
    fread(&overflowSize, sizeof(unsigned int), 1, snapshot);
    fread(&overflowLimit, sizeof(unsigned int), 1, snapshot);
    snpStatBundle->snpOverflowCounterArray[i].limit = overflowLimit;
    snpStatBundle->snpOverflowCounterArray[i].size = overflowSize;
    snpStatBundle->snpOverflowCounterArray[i].counters = (SnpOverflowCounter *) malloc(
        overflowLimit * sizeof(SnpOverflowCounter));
    if (snpStatBundle->snpOverflowCounterArray[i].counters == NULL) {
      fprintf(stderr, "Memory Allocation Failed!\n\n");
      exit(1);
    }
    if (overflowLimit == 0) continue;
    fread(snpStatBundle->snpOverflowCounterArray[i].counters, sizeof(SnpOverflowCounter), overflowSize, snapshot);
  }

  fseek(snapshot, subfile_start[3], SEEK_SET);
  snpStatBundle->invalidSnpCounterPos = (unsigned int *) malloc(((textLength + 31) >> 5) * sizeof(unsigned int));
  if (subfile_length[3]) {
    fread(snpStatBundle->invalidSnpCounterPos, sizeof(unsigned int), (textLength + 31) >> 5, snapshot);
  } else {
    memset(snpStatBundle->invalidSnpCounterPos, 0, sizeof(unsigned int) * ((textLength + 31) >> 5));
  }

  fseek(snapshot, subfile_start[0], SEEK_SET);
  unsigned int curPtr;
  size_t capacity;
  fread(&(curPtr), sizeof(unsigned int), 1, snapshot);
  fread(&(capacity), sizeof(size_t), 1, snapshot);
  
  snpStatBundle->snpMemoryPool = createPool(capacity);
  MemoryPool *&pool = snpStatBundle->snpMemoryPool;
  fread(pool->address, sizeof(char), pool->capacity, snapshot);
  pool->curPtr = curPtr;

  fseek(snapshot, subfile_start[4], SEEK_SET);
  unsigned int snpDscSize;
  snpDscSize = subfile_length[4] / sizeof(SnpDirectionalSoftclipCounter);
  snpStatBundle->snpDirectionalSoftclipCounterCapacity = snpDscSize;
  SnpDirectionalSoftclipCounter *dscCounter = (SnpDirectionalSoftclipCounter *) malloc(
      snpDscSize * sizeof(SnpDirectionalSoftclipCounter));
  fread(dscCounter, sizeof(SnpDirectionalSoftclipCounter), snpDscSize, snapshot);
  snpStatBundle->snpDirectionalSoftclipCounterSize = snpDscSize;
  snpStatBundle->snpDirectionalSoftclipCounter = dscCounter;
  fclose(snapshot);

}

