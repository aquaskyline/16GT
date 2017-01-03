

#ifndef __HSP_H__
#define __HSP_H__

#include "SNP.h"
#include "definitions.h"
#include <pthread.h>

#define MAX_SEQ_NAME_LENGTH         256

typedef struct occRec
{
  unsigned int readID;
  unsigned int ambPosition;
  unsigned char strand;
  unsigned char source;
  char score;

} occRec;


typedef struct AlgnResult
{

  occRec *occ_list;
  unsigned int occTotalNum;
  unsigned int availableSize;
} AlgnResult;


typedef struct AlgnResultArrays
{
  AlgnResult **algnArrays;
  unsigned int numArrays;
} AlgnResultArrays;

#define CIGAR_MATCHMISMATCH 'M'
#define CIGAR_INSERT 'I'
#define CIGAR_DELETE 'D'
#define CIGAR_SOFT_CLIP 'S'

typedef struct RecalQ
{

  unsigned int mismatch[Q_SCORE_ARRAY_SIZE];
  unsigned int nucleotides[Q_SCORE_ARRAY_SIZE];
} RecalQ;

#define NUM_BUCKETS 65536

#pragma pack(1)
typedef struct PairEndAlignmentInfo
{
  unsigned int read_id;

  unsigned int pos;
  unsigned short diff;

  unsigned char mapq;
  unsigned char avgBaseQuality;
} PEEntry;
#pragma pack()

#pragma pack(1)
typedef struct SingleEndAlignmentInfo
{
  unsigned int read_id;

  unsigned int pos;

  unsigned char mapq;
  unsigned char avgBaseQuality;
} SEEntry;
#pragma pack()

#include "SimpleMemoryPool.h"
#include "readIndex.h"

#define GRID_SAMPLING_FACTOR_2_POWER 18

typedef struct SeqActualOffset
{
  unsigned int startPos;
  unsigned int endPos;
} SeqActualOffset;

typedef struct Ambiguity
{
  unsigned int startPos;
  unsigned int rightOfEndPos;
  int symbol;
} Ambiguity;

typedef struct CIResult
{
  char *cigarString;
  unsigned int hitPos;
  int strand;
  unsigned int readID;
  int DPScore;
  char *altPattern;
} CIResult;

typedef struct HSP
{
  int numOfSeq;
  int numOfAmbiguity;
  unsigned int dnaLength;
  unsigned int minSeqLength;
  unsigned int *packedDNA;
  Annotation *annotation;
  SeqOffset *seqOffset;
  Ambiguity *ambiguity;
  unsigned int numOfRemovedSegment;
  unsigned int numOfGridEntry;
  unsigned int *ambiguityMap;
  Translate *translate;
  SeqActualOffset *seqActualOffset;


  AlgnResultArrays *algnResultArrays;
  char isFastq;
  int dpMatchScore;
  int dpMisMatchScore;
  int alignmentType;
  int readType;
  int peMaxOutputPerRead;
  int *x0_array;
  int *x1_array;
  int *mismatch_array;
  int minMAPQ;
  int maxMAPQ;
  int bwaLikeScore;
  int g_log_n[256];
  char *readGroup;
  char *sampleName;
  char *readGrpOption;
  int maxLenReadName;
  int ProceedDPForTooManyHits;
  int singleDPcutoffThreshold;


  void *bwt;
  void *soap3AnsArray;
  unsigned int *sa_start;
  unsigned int *occ_start;
  unsigned int *sa_num;
  unsigned int *occ_num;



  void *readsIDForSingleDP;
  void *allHits;






  char snpFlag;

  VcSetting *vcSetting;




  SnpBundle snpBundle;
  LongSoftClipPositionArray *longSoftClipArray;
  LongSoftClipPositionArray *dupLongSoftClipArray;
  LongSoftClipCounter *longSoftClipCounter;

  unsigned int mapqThreshold;


  unsigned int longSoftClipThreshold;
  unsigned int longIndelThreshold;
  unsigned int mBestLongIndelThreshold;
  unsigned int mBestLongIndelSupport;
  unsigned int mBestShortIndelThreshold;
  unsigned int mBestShortIndelSupport;
  unsigned int windowGapFactor;

  unsigned char *nonDupIDs;


  unsigned int maxBucketFiles;
  unsigned int *readsOffsetMap;

  FILE **dedupOutputFiles;
  FILE *dedupOutputDPFile;
  FILE *dedupOutputUnpairFile;

  int isBalsaOutputBam;


  unsigned int accumReadNum;







  char scoreRecalStage;
  RecalQ *recalQ;

  unsigned char **scoreMap;
  unsigned char currLane;

  unsigned int numAlignment[MAX_NUM_CPU_THREADS + 2];

  FILE *snpOverflowFilePtr[MAX_NUM_CPU_THREADS];
  FILE *snpOverflowFileDpPtr;
  FILE *snpOverflowFileUnpairPtr;

  FILE *snpDscFilePtr[MAX_NUM_CPU_THREADS];
  FILE *snpDscFileDpPtr;
  FILE *snpDscFileUnpairPtr;

  FILE *snpDscDedupFilePtr[MAX_NUM_CPU_THREADS];
  FILE *snpDscDedupFileDpPtr;
  FILE *snpDscDedupFileUnpairPtr;

  FILE *snpDscIncFilePtr;
  FILE *snpDscDecFilePtr;


  RA_Window *ra_window;


  unsigned int *window_id_map;
  unsigned int window_id_map_size;

  unsigned int *insert_high;
  unsigned int *insert_low;


  RA_AlignmentInfo *ra_alignmentInfo;
  unsigned int ra_alignmentInfo_size;


  SC_DP_AlignmentInfo *sc_dp_alignmentInfo;
  unsigned int sc_dp_alignmentInfo_size;


  RA_DP_AlignmentInfo *ra_dp_alignmentInfo;
  unsigned int ra_dp_alignmentInfo_size;


  PEEntry *peResult;
  unsigned int pe_size;
  SEEntry *seResult;
  unsigned int se_size;

  unsigned int *bucketSize;
  unsigned int bucketWidth;


  PrimerHash *primerHash;


  CIResult *ciResult;
  int curNumOfCiResult;
} HSP;

#define MAX_SEQ_NAME_LENGTH             256

static const char dnaChar[16] = {'A', 'C', 'G', 'T', 'M', 'R', 'S', 'V', 'W', 'Y', 'H', 'K', 'D', 'B', 'N', 'L'};

void HSPFree(HSP *hsp, const unsigned int trailerBufferInWord, char isShareIndex);

void getChrAndPos(unsigned int *occAmbiguityMap, Translate *occTranslate,
                  unsigned long long ambPos,
                  unsigned long long *tp, unsigned short *chr_id);

#endif

