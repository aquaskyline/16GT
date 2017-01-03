#ifndef LIKELIHOOD_CACHE_H
#define LIKELIHOOD_CACHE_H

#include <stdlib.h>
#include <omp.h>

#include "fisher.h"
#include "ycsq.h"

#define FISHER_CACHE_SIZE 2000
#define YCSQ_CACHE_SIZE 15000

typedef struct FisherCache
{
  double cache1[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache2[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache3[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache4[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache5[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache6[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache7[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache8[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache9[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
  double cache10[FISHER_CACHE_SIZE * FISHER_CACHE_SIZE];
} FisherCache;

typedef struct LikelihoodCache
{
  FisherCache fc;
  double ycsq_cache[YCSQ_CACHE_SIZE];
} LikelihoodCache;

int likelihood_cache_prefill(LikelihoodCache *lc,
                             float homo_base_err, float homo_base_e_err,
                             float homo_indel_err, float homo_indel_e_err,
                             float het_base_err, float het_base_e_err,
                             float het_indel_1_err, float het_indel_e_1_err,
                             float het_indel_2_err, float het_indel_e_2_err);

double fisher_homo_base(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_homo_indel(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_homo_base_e(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_homo_indel_e(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_het_base(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_het_indel_1(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_het_indel_2(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_het_base_e(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_het_indel_e_1(LikelihoodCache *lc, int w_target, int err, int w_total);

double fisher_het_indel_e_2(LikelihoodCache *lc, int w_target, int err, int w_total);

#endif
