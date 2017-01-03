#ifndef _FISHER_EXACT_TEST_H_
#define _FISHER_EXACT_TEST_H_

#include "SNP_Meta.h"
#include "definitions.h"
#include "likelihood_cache.h"

#define USE_FISHER_LOOKUP_TABLE

void copyGenotypeTableToGPU ( unsigned int *genotype2int );

void copySnpPriorToGPU ( double *snpPrior );

void computeGenotypeLikelihoodWrapper ( InputStat *in, MetaReference *reference, GenotypeLikelihood *lh, IniParams ini_params,
                                        unsigned int input_size, LikelihoodCache *lc );

void computeGenotypeWrapper ( unsigned char *refChar,
                              GenotypeLikelihood *lh, SnpCallingInfo *info, IniParams ini_params,
                              unsigned int input_size );

void computeStrandBias ( InputStat *in, StrandBias *sb, IniParams ini_params, unsigned int input_size );

int prefill_likelihood_cache_with_p_err ( LikelihoodCache *lc, float b_err, float ub_err  );

#endif
