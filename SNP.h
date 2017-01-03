#ifndef _SNP_H_
#define _SNP_H_

#include <stdio.h>
#include <stdlib.h>
#include "coreStruct.h"

#define RA_INDEL_PATTERN_SEQ_LENGTH 32

#pragma pack(1)
typedef struct SnpOverflowBuffer
{
  unsigned int position;
  char strand;
  unsigned char info;



  unsigned char weight;
  unsigned char length;
  unsigned char insertSeq[100];
} SnpOverflowBuffer;
#pragma pack()

#pragma pack(4)
typedef struct SnpOverflowBufferArray
{
  SnpOverflowBuffer *buffer;
  unsigned int size;
  unsigned int limit;
} SnpOverflowBufferArray;
#pragma pack()

#pragma pack(4)
typedef struct SnpUpdateOverflowWrapperObj
{
  SnpBundle snpBundle;
  VcSetting *vcSetting;
  unsigned int threadId;
} SnpUpdateOverflowWrapperObj;
#pragma pack()

#pragma pack(1)
typedef struct RA_Window
{
  unsigned int startPos;
  unsigned int endPos;
  unsigned int raInfoStart;
  unsigned char condFlag;

} RA_Window;
#pragma pack()

#pragma pack(1)
typedef struct RA_AlignmentInfo
{
  unsigned int readID;
  unsigned int windowID;
  unsigned int startPos;
  unsigned int readsOffset;
  unsigned char indelPattern[RA_INDEL_PATTERN_SEQ_LENGTH];
  unsigned char cigar[RA_CIGAR_LENGTH];
  unsigned char strand;
  unsigned char dpScore;
} RA_AlignmentInfo;
#pragma pack()

#pragma pack(1)
typedef struct SC_DP_AlignmentInfo
{
  unsigned int windowID;
  unsigned int readID;
  unsigned int maxIndelLength;
  unsigned char indelPattern[RA_INDEL_PATTERN_SEQ_LENGTH];
} SC_DP_AlignmentInfo;
#pragma pack()

#pragma pack(1)
typedef struct RA_DP_AlignmentInfo
{
  unsigned int readID;
  unsigned int windowID;
  unsigned int modifiedRefID;
  unsigned int startPos;
  unsigned char cigar[RA_CIGAR_LENGTH];
  unsigned char dpScore;
} RA_DP_AlignmentInfo;
#pragma pack()



void destroySnpCounter(SnpCounter *snpCounter, unsigned int dnaLength);

void constructDirectionalSoftclipCounterFromFiles(int numOfCPUThreads,
    FILE **snpDirectionalSoftclipFilePtr,
    FILE *&snpDirectionalSoftclipFileDpPtr,
    FILE *&snpDirectionalSoftclipFileUnpairPtr,
    FILE **snpDirectionalSoftclipDedupFilePtr,
    FILE *&snpDirectionalSoftclipDedupFileDpPtr,
    FILE *&snpDirectionalSoftclipDedupFileUnpairPtr,
    SnpDirectionalSoftclipCounter *&snpDirectionalSoftclipCounter,
    unsigned int &snpDirectionalSoftclipCounterSize,
    unsigned int &snpDirectionalSoftclipCounterCapacity);


void constructSnpOverflowCounter(SnpOverflowCounterArray **snpOverflowCounterArrayPtr, int numOfCPUThreads);

unsigned int addSnpOverflowCounterToArray(SnpOverflowCounterArray &snpOverflowCounterArray, int numOfCPUThreads);

void destroySnpOverflowCounter(SnpOverflowCounterArray *snpOverflowCounterArray, int numOfCPUThreads);

void destroyMemoryPool(MemoryPool *pool);

void constructSnpOverflowBufferArray(int numOfCPUThreads);

void emptySnpOverflowBuffer(int numOfCPUThreads);

void destroySnpOverflowBufferArray(int numOfCPUThreads);


void updateSnpCounterForReads(SnpBundle snpBundle, VcSetting *vcSetting, int trimHeadSize, int trimTailSize,
                              unsigned char *recalScores,
                              LongSoftClipPositionArray *longSoftClipArray, unsigned int longSoftClipThreshold,
                              unsigned int position, const char strand,
                              unsigned char *query, char *qualities, unsigned int readLength,
                              char *cigar, FILE *overflowFilePtr, FILE *softclipFilePtr);

char updateSnpCounter(SnpBundle snpBundle, VcSetting *vcSetting,
                      unsigned int position, char base, char weight, char strand, char strandCount);

void readFilesForOverflow(SnpBundle snpBundle, VcSetting *vcSetting,
                          const char *filename);

void startUpdateOverflowCounter(SnpBundle snpBundle, VcSetting *vcSetting);


void updateDupSnpCounter(SnpBundle snpBundle, VcSetting *vcSetting,
                         unsigned int position, char base, char weight,
                         char strand, char strandCount);

void updateDupSnpCounterForInsert(SnpBundle snpBundle, VcSetting *vcSetting,
                                  unsigned int position, unsigned char insertLength,
                                  unsigned char avgWeight, unsigned char *insertSeq);

void updateDupSnpCounterForDelete(SnpBundle snpBundle, VcSetting *vcSetting,
                                  unsigned int position, unsigned char deleteLength, char weight);


#endif
