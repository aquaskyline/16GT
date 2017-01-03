

#ifndef __BWT_H__
#define __BWT_H__

#include <stdio.h>
#include <stdlib.h>
#include "HSP.h"

const size_t BWT_HOLLOW_BETWEEN_METADATA_PAYLOAD = 3;
const size_t OCC_HOLLOW_BETWEEN_METADATA_PAYLOAD = 3;

typedef struct BWT
{
  unsigned int textLength;
  unsigned int saInterval;
  unsigned int inverseSaInterval;
  unsigned int inverseSa0;
  unsigned int *cumulativeFreq;
  unsigned int *bwtCode;
  unsigned int *occValue;
  unsigned int *occValueMajor;
  unsigned int *saValue;
  unsigned int *inverseSa;
  unsigned int *cachedSaIndex;
  unsigned int cachedSaIndexNumOfChar;
  unsigned int *saValueOnBoundary;
  unsigned int *decodeTable;
  unsigned int decodeTableGenerated;
  unsigned int bwtSizeInWord;
  unsigned int occSizeInWord;
  unsigned int occMajorSizeInWord;
  unsigned int saValueSizeInWord;
  unsigned int inverseSaSizeInWord;
  unsigned int cachedSaIndexSizeInWord;
} BWT;

void BWTFree(BWT *bwt, char isShareIndex);

#endif
