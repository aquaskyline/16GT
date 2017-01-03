

#ifndef _CPUFUNCTIONS_H_
#define _CPUFUNCTIONS_H_

#include "dependencies.h"
#include "bwt.h"


typedef struct MultiInputItem
{
  char queryFile1[MAX_FIELD_LEN];
  char queryFile2[MAX_FIELD_LEN];
  int insert_low;
  int insert_high;
  char bamOutputPrefix[MAX_FIELD_LEN];
  char readGrpID[MAX_FIELD_LEN];
  char sampleName[MAX_FIELD_LEN];
  char readGrpOpt[MAX_FIELD_LEN];
} MultiInputItem;



MultiInputItem *loadMultiInputFile(char *fileName, int isPair, int isBAM, int &lineNum);



void updateInputOption(InputOptions *input_options, MultiInputItem *multiInputArray, int i);


void processIndexFileName(IndexFileNames &index, char *indexName, IniParams ini_params);

void loadIndex(IndexFileNames index, BWT **bwt, BWT **revBwt, HSP **hsp, uint **lkt, uint **revLkt,
               uint **revOccValue, uint **occValue, uint &numOfOccValue, uint &lookupWordSize, char isShareIndex);


#endif
