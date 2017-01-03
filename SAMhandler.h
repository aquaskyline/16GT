

#ifndef __SAM_H__
#define __SAM_H__

#include <stdio.h>
#include <stdlib.h>
#include "HSP.h"
#include "definitions.h"

#include "sam.h"

#define OCC_CACHE_SIZE               81920

typedef struct OCCPositionCache
{
  unsigned short ChromId;
  unsigned char ReadStrand;
  unsigned long long tp;
  int occMismatch;
  int resultSource;
  char *cigarString;
  unsigned int len;
} OCCPositionCache;

typedef struct OCCPositionCacheToDisk
{
  unsigned char cell[11];
} OCCPositionCacheToDisk;

typedef struct OCC
{
  unsigned int occPositionCacheCount;
  OCCPositionCache occPositionCache[OCC_CACHE_SIZE];
  OCCPositionCacheToDisk occPositionCacheToDisk[OCC_CACHE_SIZE];

  bam1_t SAMOutBuffer;
} OCC;

#define SAM_MDATA_SIZE                  2048







int SAMIUint8ConcatUint8(uint8_t *data, int *curSize,
                         uint8_t key);

int SAMIUint8ConcatUint32(uint8_t *data, int *curSize,
                          uint32_t key);

int SAMIUint8ConcatString(uint8_t *data, int *curSize,
                          char *key, int len);

void SAMOutputHeaderConstruct(bam_header_t *sheader, HSP *hsp);

void SAMOutputHeaderDestruct(bam_header_t *sheader);

void SAMOccurrenceConstruct(OCC *occ);

void SAMOccurrenceDestruct(OCC *occ);

#endif
