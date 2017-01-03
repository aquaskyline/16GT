#ifndef __COUNTER_READER_H__
#define __COUNTER_READER_H__

#include <cstdio>
#include <cstring>

#include "coreStruct.h"

typedef struct InsertSeqs
{
  char **insertSeq;
  unsigned int *hqCount;
  unsigned int *lqCount;
  int size;
} InsertSeqs;

typedef struct DeleteSeqs
{
  char *deleteLength;
  unsigned int *hqCount;
  unsigned int *lqCount;
  int size;
} DeleteSeqs;

typedef struct SnpCounterInfo
{
  unsigned int pos;
  char overflow;
  SnpOverflowCounter counter;

  InsertSeqs *insSeqs;
  DeleteSeqs *delSeqs;

  unsigned int leftSoftclip;
  unsigned int rightSoftclip;
} SnpCounterInfo;

#endif
