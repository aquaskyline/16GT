#ifndef _SNP_STAT_H_
#define _SNP_STAT_H_

#include <string.h>

#include "SNP_Meta.h"


#define DEPTH_THRES 65536
#define STEPS 100

#define UG UNIDENTIFIED_GENOTYPE
static const int Base2GTMap[] =
{
  0, 1, 2, 3, 10, UG,
  1, 4, 5, 6, 11, UG,
  2, 5, 7, 8, 12, UG,
  3, 6, 8, 9, 13, UG,
  10, 11, 12, 13, 14, 15,
  UG, UG, UG, UG, 15, UG
};
#undef UG

__inline__ int getGenotype(unsigned char b1, unsigned char b2)
{
  return Base2GTMap[b1 * COUNTER_NUM + b2];
}

typedef struct WholeGenomeMeta
{
  unsigned int zeroCount;
  unsigned int depthTooHigh;
  unsigned int depthTooLow;
  unsigned int hom;
  unsigned int het;
  unsigned int tri;
  unsigned int quad;
} WholeGenomeMeta;

typedef struct SnpMeta
{
  unsigned int zeroCount;
  unsigned int depthTooHigh;
  unsigned int depthTooLow;
  unsigned int hom;
  unsigned int het;
  unsigned int tri;
  unsigned int quad;
  unsigned int inconsistentDBnConsensusHom;
  unsigned int notConsensusButRefHom;
  unsigned int inconsistentDBnConsensusHet;
  unsigned int valid;
  unsigned int invalid;
} SnpMeta;

typedef struct SnpStats
{
  unsigned int totalDepth[DEPTH_THRES + 1];
  unsigned int baseDepth[DEPTH_THRES + 1];
  unsigned int lh[STEPS + 1];
  unsigned int sb[STEPS + 1];
  unsigned int bqb[STEPS + 1];
} SnpStats;

typedef struct SnpProbs
{
  double totalDepth[DEPTH_THRES + 1];
  double baseDepth[DEPTH_THRES + 1];
  double lh[STEPS + 1];
  double sb[STEPS + 1];
  double bqb[STEPS + 1];
} SnpProbs;


__inline__ int pVal2Index(double p)
{
  return int(p * STEPS + 0.4999);
}

__inline__ int pValf2Index(float p)
{
  return int(p * STEPS + 0.4999f);
}

#pragma pack(push)
#pragma pack(4)
typedef struct FourBasesDepth
{
  unsigned short depth;
  short id;
} FourBasesDepth;
#pragma pack(pop)

static __inline__ void sort4(FourBasesDepth *d)
{
#define mind(x, y) ((x.depth<y.depth)?x:y)
#define maxd(x, y) ((x.depth<y.depth)?y:x)
#define SWAP(x, y, _d) { FourBasesDepth tmp = maxd(_d[x], _d[y]); _d[y] = mind(_d[x], _d[y]); _d[x] = tmp; }

  SWAP(0, 1, d);
  SWAP(2, 3, d);
  SWAP(0, 2, d);
  SWAP(1, 3, d);
  SWAP(1, 2, d);
#undef SWAP
#undef maxd
#undef mind
}

__inline__ void updateToStats(const FourBasesDepth *const d, int homOrHet,
                              SnpStats *const snpStats,
                              const GenotypeLikelihood *const genotypeLikelihood,
                              const StrandBias *const strandBias)
{
  if (homOrHet == 0)
    {
      ++snpStats->lh[pVal2Index(genotypeLikelihood->likelihood[getGenotype(d[0].id, d[0].id)])];
      ++snpStats->sb[pValf2Index(strandBias->bias[d[0].id])];
    }
  else if (homOrHet == 1)
    {
      ++snpStats->lh[pVal2Index(genotypeLikelihood->likelihood[getGenotype(d[0].id, d[1].id)])];
      ++snpStats->sb[pValf2Index((strandBias->bias[d[0].id] + strandBias->bias[d[1].id]) * 0.5)];
    }
  else
    {
      fprintf(stderr, "[%s][%s][%d] Should never reach here.\n", __FILE__, __FUNCTION__, __LINE__);
      exit(1);
    }
}

__inline__ int myStrChr(const char *const str, const char c)
{
  int x = strlen(str);

  for (int i = 0; i < x; ++i)
    {
      if (str[i] == c)
        {
          return 1;
        }
    }

  return 0;
}


#endif
