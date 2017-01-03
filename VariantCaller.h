#ifndef _VARIANT_CALLER_H_
#define _VARIANT_CALLER_H_

#include <stdio.h>
#include <stdlib.h>

#include "HSP.h"
#include "dependencies.h"
#include "coreStruct.h"
#include "SNP.h"
#include "SNP_Meta.h"

#include "likelihood_cache.h"
#include "definitions.h"
#include "struct.h"

char getSnpMetrics(unsigned int *packedSeq,
                   SnpBundle snpBundle, unsigned int position, unsigned int dnaLength,
                   InputStat &iStat, unsigned int &weightedCount,
                   unsigned int &posStrandCount, unsigned int &negStrandCount,
                   unsigned char &gcCount, unsigned char &gCount,
                   unsigned int &indelHqCount, unsigned int &indelLqCount,
                   IndelInfo *indelInfo);


void selectPossibleSNPs(SnpBundle snpBundle, unsigned int *packedSeq, unsigned int dnaLength, Annotation *annotation,
                        unsigned int *attriSize, char isExome, ExomeRegion *exomeRegion,
                        unsigned int numExomeRegion, LikelihoodCache *likelihood_cache, InputOptions &input_options,
                        IniParams ini_params, unsigned int *ambMap, Translate *translate, const char *snp_noRF_filename,
                        double startTime, double &lastEventTime);

void calRFpedictProb(const char *snp_noRF_filename,
                     const char *snp_filename,
                     unsigned int attriSize, char verbose);

#endif
