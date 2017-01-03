#include "likelihood_cache.h"

#define FISHER_CACHE_NOT_CACHED -1.0

int fisher_cache_prefill(FisherCache *fisher_cache,
                         float homo_base_err, float homo_base_e_err,
                         float homo_indel_err, float homo_indel_e_err,
                         float het_base_err, float het_base_e_err,
                         float het_indel_1_err, float het_indel_e_1_err,
                         float het_indel_2_err, float het_indel_e_2_err) {
  if (!fisher_cache) {
    return 0;
  }

  int i;
#pragma omp parallel for
  for (i = 0; i < FISHER_CACHE_SIZE; i++) {
    int j;
    for (j = 0; j < FISHER_CACHE_SIZE; j++) {
      int idx = i * FISHER_CACHE_SIZE + j;
      fisher_cache->cache1[idx] = fisher2(j, i, i * homo_base_err + 0.5f, i);
      fisher_cache->cache2[idx] = fisher2(j, i, i * homo_base_e_err + 0.5f, i);
      fisher_cache->cache3[idx] = fisher2(j, i, i * homo_indel_err + 0.5f, i);
      fisher_cache->cache4[idx] = fisher2(j, i, i * homo_indel_e_err + 0.5f, i);
      fisher_cache->cache5[idx] = fisher2(j, i, i * het_base_err + 0.5f, i);
      fisher_cache->cache6[idx] = fisher2(j, i, i * het_base_e_err + 0.5f, i);
      fisher_cache->cache7[idx] = fisher2(j, i, i * het_indel_1_err + 0.5f, i);
      fisher_cache->cache8[idx] = fisher2(j, i, i * het_indel_e_1_err + 0.5f, i);
      fisher_cache->cache9[idx] = fisher2(j, i, i * het_indel_2_err + 0.5f, i);
      fisher_cache->cache10[idx] = fisher2(j, i, i * het_indel_e_2_err + 0.5f, i);
    }
  }
  return 1;
}

int likelihood_cache_prefill(LikelihoodCache *lc,
                             float homo_base_err, float homo_base_e_err,
                             float homo_indel_err, float homo_indel_e_err,
                             float het_base_err, float het_base_e_err,
                             float het_indel_1_err, float het_indel_e_1_err,
                             float het_indel_2_err, float het_indel_e_2_err) {

  fisher_cache_prefill(&(lc->fc),
                       homo_base_err, homo_base_e_err,
                       homo_indel_err, homo_indel_e_err,
                       het_base_err, het_base_e_err,
                       het_indel_1_err, het_indel_e_1_err,
                       het_indel_2_err, het_indel_e_2_err);

#pragma omp parallel for
  for (int i = 0; i < YCSQ_CACHE_SIZE; i++) {
    lc->ycsq_cache[i] = chisqr(1, i / 10.0);
  }
}


static double ycsq_p(LikelihoodCache *lc, int a, int b, int c, int d) {

  double dist = chi2UniformDistance(a, b, c, d);

  int dist_int = (int) (dist * 10 + 0.5);
  if (dist_int >= YCSQ_CACHE_SIZE) {
    dist_int = YCSQ_CACHE_SIZE - 1;
  }
  return lc->ycsq_cache[dist_int];
}

double fisher_homo_base(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache1[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache1[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }
  return result;
}

double fisher_homo_indel(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache3[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache3[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_homo_base_e(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache2[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache2[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }
  return result;
}

double fisher_homo_indel_e(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache4[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache4[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_het_base(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache5[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache5[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_het_indel_1(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache7[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache7[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_het_indel_2(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache9[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache9[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_het_base_e(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache6[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache6[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_het_indel_e_1(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache8[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache8[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

double fisher_het_indel_e_2(LikelihoodCache *lc, int w_target, int err, int w_total) {
  FisherCache *fc = &(lc->fc);
  if (!fc || w_target >= FISHER_CACHE_SIZE || w_total >= FISHER_CACHE_SIZE) {
    return ycsq_p(lc, w_target, w_total, err, w_total);
  }


  int cache_idx = w_total * FISHER_CACHE_SIZE + w_target;
  double result = fc->cache10[cache_idx];
  if (result == FISHER_CACHE_NOT_CACHED) {
    result = fc->cache10[cache_idx] =
        fisher2(w_target, w_total, err, w_total);
  }

  return result;
}

