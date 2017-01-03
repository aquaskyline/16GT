#ifndef _DEPENDENCIES_H_
#define _DEPENDENCIES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "definitions.h"

typedef unsigned int uint;

typedef struct InputOptions
{

  char *queryFileName;
  char *queryFileName2;
  int maxReadLength;

  int outputFormat;
  int isOutputBinary;





  int isBalsaOutputBam;





  char *outputPrefix;
  char *bamOutputPrefix;

  char readType;
  int numMismatch;
  int alignmentType;

  int insert_low;
  int insert_high;

  char enableDP;

  int maxHitNum;
  int maxHitNum2;

  char fileType;

  char isIlluminaQual;

  char isBAM;

  int GPUDeviceID;

  char *readGroup;

  char *sampleName;

  char *readGrpOption;

  char *dbSnpIndexFileName;

  char *indelDBIndexFileName;

  char *geneRegionFileName;

  char *resultPrefix;

  char isExome;
  char *exomeRegionFileName;

  char *tempPrefix;

  char outputSnapshot;

  char enableRandomForest;

  char enableQualityCorrection;

  char *primerListFileName;

  char verbose;

} InputOptions;


typedef struct IniParams
{
  char Ini_SaValueFileExt[MAX_FILEEXT_LEN];
  int Ini_NumOfCpuThreads;
  char Ini_HostAlignmentModelStr[4];
  int Ini_HostAlignmentModel;
  int Ini_GPUMemory;
  int Ini_PEStrandLeftLeg;
  int Ini_PEStrandRightLeg;
  unsigned int Ini_MaxOutputPerRead;
  unsigned int Ini_PEMaxOutputPerPair;
  unsigned int Ini_MaxHitsEachEndForPairing;
  int Ini_MatchScore;
  int Ini_MismatchScore;
  int Ini_GapOpenScore;
  int Ini_GapExtendScore;
  int Ini_DPScoreThreshold;
  int Ini_DiagonalWidth;
  int Ini_isDefaultThreshold;
  int Ini_Soap3MisMatchAllow;
  int Ini_maxMAPQ;
  int Ini_minMAPQ;
  int Ini_shareIndex;
  int Ini_maxReadNameLen;
  int Ini_maxFrontLenClipped;
  int Ini_maxEndLenClipped;
  int Ini_proceedDPForTooManyHits;
  int Ini_skipSOAP3Alignment;
  int Ini_bwaLikeScore;

  unsigned char Ini_weightBound0;
  unsigned char Ini_weightBound1;
  unsigned char Ini_weightBound2;
  unsigned char Ini_weightBound3;
  unsigned char Ini_IndelQualityThreshold;

  unsigned int Ini_memoryPoolSize;
  unsigned int Ini_statMapqThreshold;
  unsigned int Ini_statTrimSize;
  unsigned int Ini_statSoftClipThreshold;
  unsigned int Ini_LongSoftClipThreshold;
  unsigned int Ini_LongIndelThreshold;

  unsigned int Ini_MaxIndelPattern;
  unsigned int Ini_RASupportThreshold;

  int Ini_ScMatchScore;
  int Ini_ScMismatchScore;
  int Ini_ScGapOpenScore;
  int Ini_ScGapExtendScore;

  int Ini_IndelCallHqThreshold;
  int Ini_IndelCallLqThreshold;

  double Ini_BalanceSubError;
  double Ini_UnbalanceSubError;

  int Ini_RAW_P[23];
  unsigned int geneRegionLongIndel;
  unsigned int geneRegionLongIndelSupportThreshold;
  unsigned int nonGeneRegionLongIndel;
  unsigned int nonGeneRegionLongIndelSupportThreshold;
  unsigned int tandemRepeatIndel;
  unsigned int tandemRepeatIndelSupportThreshold;
  unsigned int tandemRepeatMmSupportThreshold;
  unsigned int tandemRepeatSearchRange;


  char Ini_skipScoreRecalibration;
  char Ini_skipDeduplication;







} IniParams;

typedef struct IndexFileNames
{
  char *iniFileName;
  char *bwtCodeFileName;
  char *occValueFileName;
  char *gpuOccValueFileName;
  char *lookupTableFileName;
  char *revBwtCodeFileName;
  char *revOccValueFileName;
  char *revGpuOccValueFileName;
  char *revLookupTableFileName;
  char *saCodeFileName;
  char *memControlFileName;


  char *packedDnaFileName;
  char *annotationFileName;
  char *ambiguityFileName;
  char *translateFileName;


  char *mmapOccValueFileName;
  char *mmapRevOccValueFileName;
  char *mmapPackedDnaFileName;

} IndexFileNames;

int ParseIniFile(char *iniFileName, IniParams &ini_params);

bool fileExists(const char *filePath);

bool dirOfPrefixExists(const char *prefix);

void *xmalloc(unsigned long long size);

void xfree(void *p);

double setStartTime();

double getElapsedTime(double startTime);

#endif
