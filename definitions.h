

#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_

#define MAX_FILEEXT_LEN 6
#define ALPHABET_SIZE 4
#define MAX_READ_LENGTH 1024
#define MAX_FIELD_LEN 256







#define PROGRAM_INFO_SAM "@PG\tID:SOAP3-dp\tPN:SOAP3-dp\tVN:v2.2\n"

#define Q_SCORE_ARRAY_SIZE 41
#define SNP_POOL_SIZE 1073741824u
#define NUM_BUCKETS 65536
#define SNP_DIRECTIONAL_SOFTCLIP
#define SNP_STAT_FLAG 2
#define RA_CIGAR_LENGTH 16
#define MAX_NUM_CPU_THREADS 99
#define MAX_PRIMER_PACKED_LENGTH 5
#define MAX_PRIMER_LENGTH 40

typedef struct Primer
{
  int length[2];
  unsigned int packedPattern[2][MAX_PRIMER_PACKED_LENGTH];
  char pattern[2][MAX_PRIMER_LENGTH];
  unsigned int primerPos[2];
  char strand;
} Primer;

#define POPCOUNT_MAP_SIZE 256

typedef struct PrimerInfo
{
  short id;
  unsigned char pos;
} PrimerInfo;


typedef struct PrimerList
{
  int size;
  PrimerInfo *primerInfo;
} PrimerList;



typedef struct PrimerGroup
{
  unsigned short *bitMap;
  unsigned char listSize;
  PrimerList *primerList;
} PrimerGroup;

typedef struct PrimerHashTable
{
  int groupSize;
  int numOfGroup;
  int bitMapLength;
  PrimerGroup *group;

  unsigned char popCountMap[POPCOUNT_MAP_SIZE];
} PrimerHashTable;

typedef struct PrimerHash
{
  PrimerHashTable *primerHashTable;
  PrimerHashTable **revPrimerHashTable;
  Primer *primerList;
  unsigned int numOfPrimerPair;
  int dnaMap[256];
  char revDnaMap[4];
} PrimerHash;


Primer *getTrimSize(PrimerHash *primerHash, unsigned char *read, int readLength, char isPosStrand,
                    char *cigarString, unsigned int alignPos,
                    int *trimHead, int *trimTail, const int minTrimSize);


#endif
