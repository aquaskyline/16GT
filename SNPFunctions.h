#ifndef _SNP_FUNCTIONS_H_
#define _SNP_FUNCTIONS_H_

#include "SNP_Meta.h"

double *createSnpPriorArray();

void freeSnpPriorArray(double *ary);

void computeStrandBias(MetaSnpCounter *snpCounter, StrandBias *strandBias);

void computeBaseQualityBias(MetaSnpCounter *snpCounter, BaseQualityBias *baseQualityBias);

#endif
