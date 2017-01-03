#ifndef _ASC_CORESTRUCT_H_
#define _ASC_CORESTRUCT_H_

#include "definitions.h"
#include "SimpleMemoryPool.h"



#pragma pack(1)
typedef struct SnpCounter
{
  unsigned short weightedCount1;
  unsigned short weightedCount2;

  union
  {
    struct
    {
      unsigned char posStrandCount1;
      unsigned char negStrandCount1;
      unsigned char posStrandCount2;
      unsigned char negStrandCount2;
    } strandCounts;

    unsigned int arrayIndex;
  } data;
#ifdef ENABLE_SOFTCLIP_COUNTER
  unsigned char softClipCount;
#endif
} SnpCounter;
#pragma pack()

#pragma pack(4)
typedef struct SnpOverflowCounter
{
#ifndef AMPLICON
  unsigned short weightedCount[ALPHABET_SIZE];
  unsigned short posStrandCount[ALPHABET_SIZE];
  unsigned short negStrandCount[ALPHABET_SIZE];
#else
  unsigned int weightedCount[ALPHABET_SIZE];
  unsigned int posStrandCount[ALPHABET_SIZE];
  unsigned int negStrandCount[ALPHABET_SIZE];
#endif
  union
  {
    struct
    {
      unsigned char insertSeq[3];
      unsigned char lqCount;
    } counters;

    unsigned int ptr;
  } insertion;

  union
  {
    struct
    {
      unsigned char deletionLength;
      unsigned char hqCount;
      unsigned char lqCount;
    } counters;

    unsigned int ptr;
  } deletion;
} SnpOverflowCounter;
#pragma pack()

#pragma pack(4)
typedef struct SnpOverflowCounterArray
{
  SnpOverflowCounter *counters;
  unsigned int size;
  unsigned int limit;
} SnpOverflowCounterArray;
#pragma pack()

#pragma pack(1)
typedef struct SnpDirectionalSoftclip
{
  unsigned int ambPosition;
  char leftOrRight;
} SnpDirectionalSoftClip;
#pragma pack()

#pragma pack(2)
typedef struct SnpDirectionalSoftclipCounter
{
  unsigned int ambPosition;
  unsigned int leftSoftclipCount;
  unsigned int rightSoftclipCount;
} SnpDirectionalSoftclipCounter;
#pragma pack()

#pragma pack(4)
typedef struct LongSoftClipPositionArray
{
  unsigned int *positions;
  unsigned int size;
  unsigned int limit;
} LongSoftClipPositionArray;
#pragma pack()

#pragma pack(1)
typedef struct LongSoftClipCounter
{
  unsigned int position;
  unsigned char count;
} LongSoftClipCounter;
#pragma pack()

#pragma pack(4)
typedef struct SnpBundle
{
  SnpCounter *snpCounter;
  SnpOverflowCounterArray *snpOverflowCounterArray;
  int numOfCPUThreads;
  unsigned int *invalidSnpCounterPos;
  MemoryPool *snpMemoryPool;
  SnpDirectionalSoftclipCounter *snpDirectionalSoftclipCounter;
  unsigned int textLength;
  unsigned int snpDirectionalSoftclipCounterSize;
  unsigned int snpDirectionalSoftclipCounterCapacity;
} SnpBundle;
#pragma pack()

#pragma pack(1)
typedef struct SnpSetting
{
  unsigned int dnaLength;
  unsigned int mapqThreshold;
  unsigned int trimSize;
  unsigned int softClipThreshold;
  unsigned char *weightMap;
  unsigned char indelWeightThreshold;
  char enableQualityCorrection;
} VcSetting;
#pragma pack()

__inline__ unsigned int findRegionByPosition(unsigned int position, unsigned int textLength, int numOfCPUThreads)
{
  int i = 0;

  for (i = 0; i < numOfCPUThreads - 1; ++i)
    {
      if (position < (textLength / numOfCPUThreads) * (i + 1))
        {
          return i;
        }
    }

  return i;
}

#define SORT_DIGIT_BY_VAR(array, var, shift_bits, num_elements, buckets, buffer, record_type) \
  do { memset(buckets, 0, NUM_BUCKETS * sizeof(unsigned int));          \
      for (unsigned int i=0; i<num_elements; ++i) buckets[((array)[i].var>>shift_bits)&0xFFFF]++; \
      for (unsigned int b=0, acc=0; b<NUM_BUCKETS; ++b) {unsigned int t=acc; acc+=buckets[b]; buckets[b]=t;} \
      for (unsigned int i=0; i<num_elements; ++i) buffer[buckets[((array)[i].var>>shift_bits)&0xFFFF]++] = array[i]; \
      { record_type *t = buffer; buffer = array; array = t; } } while (0)


#endif
