#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <omp.h>

#include "dependencies.h"
#include "FisherExactTest.h"

#define errorJitterStart 5
#define errorJitterEnd 6
#define geometricMeanRoot 0.5f
#define indelDifficultyCoefficient 1.375f

unsigned int gpuGenotype2int[16];

double gpuSnpPrior[256];

void copyGenotypeTableToGPU(unsigned int *genotype2int) {
  memcpy(gpuGenotype2int, genotype2int, 16 * sizeof(unsigned int));
}

void copySnpPriorToGPU(double *snpPrior) {
  memcpy(gpuSnpPrior, snpPrior, 256 * sizeof(double));
}

float fisher2tail(int a, int b, int c, int d);

float fisher2tail(int n11, int n1_, int n_1, int n) 
{
  n += n11 + n1_ + n_1;
  n1_ += n11;
  n_1 += n11;

  float poa;
  poa = hyper_323(n11, n1_, n_1, n); 

  int x;
  float pox = poa;
  float p_twotail = 0;
  p_twotail += poa;
  float poa_e = poa * 1.000000001f;

  int min = (int) fmaxf(0, n1_ + n_1 - n);
  for (x = n11; x > min; --x) {
    pox = pox * (float) (x * (n + x - n1_ - n_1)) / ((n1_ - x + 1) * (n_1 - x + 1));
    p_twotail += pox * (pox < poa_e);
  }

  int max = (int) fminf(n1_, n_1);
  pox = poa;
  for (x = n11 + 1; x <= max; ++x) {
    pox = pox * (float) (n1_ - (x - 1)) * (n_1 - (x - 1)) / ((n + x - n1_ - n_1) * x);
    p_twotail += pox * (pox < poa_e);
  }

  return p_twotail;

}

int prefill_likelihood_cache_with_p_err(LikelihoodCache *lc, float b_err, float ub_err) {

  float P_err_s = b_err / ALPHABET_SIZE;
  float P_err_i = ub_err / 2;

  float P_err_base = 3 * P_err_s + 2 * P_err_i;
  float hom_base_err = 1.0f - P_err_base;
  float hom_base_e_err = P_err_base;

  float P_err_indel = 4 * P_err_s + P_err_i;
  float hom_indel_err = 1.0f - P_err_indel;
  float hom_indel_e_err = P_err_indel;

  float P_err_base_base = 2 * P_err_s + 2 * P_err_i;
  float het_base_err = 0.5f - P_err_base_base;
  float het_base_e_err = P_err_base_base;

  float P_err_base_indel = 3 * P_err_s + P_err_i;
  float het_indel_1_err = 0.5f - P_err_base_indel;
  float het_indel_1_e_err = P_err_base_indel;

  float P_err_indel_indel = 4 * P_err_s;
  float het_indel_2_err = 0.5f - P_err_indel_indel;
  float het_indel_2_e_err = P_err_indel_indel;

  return likelihood_cache_prefill(lc,
                                  hom_base_err, hom_base_e_err,
                                  hom_indel_err, hom_indel_e_err,
                                  het_base_err, het_base_e_err,
                                  het_indel_1_err, het_indel_1_e_err,
                                  het_indel_2_err, het_indel_2_e_err);
}

#ifdef USE_FISHER_LOOKUP_TABLE

#define OVERFLOW_EXP 849
#define RET_EXP_LIMIT -1500

#define OVERFLOW_HANDLER 3.7537584144023501e+255

double checkOverflowAndReturn(double ret) {
  int ret_exp;
  frexp(ret, &ret_exp);
  if (ret == .0 || ret_exp - OVERFLOW_EXP < RET_EXP_LIMIT) {

    return 0;
  }

  double downRet = ret / OVERFLOW_HANDLER;
  if (downRet != .0) {
    return downRet;
  }
  return DBL_MIN;
}

static void GenerateIntegerArrayFlanking5(int num, int *ary) {
  ary[0] = num - 5;
  ary[1] = num - 4;
  ary[2] = num - 3;
  ary[3] = num - 2;
  ary[4] = num - 1;
  ary[5] = num;
  ary[6] = num + 1;
  ary[7] = num + 2;
  ary[8] = num + 3;
  ary[9] = num + 4;
  ary[10] = num + 5;
  for (int i = 0; i < 11; ++i) {
    if (ary[i] < 0) {
      ary[i] = 0;
    }
  }
}

double Likelihood_homo_base(unsigned int W[], int targetChar, int W_total,
                            double P_err_s, double P_err_i, LikelihoodCache *lc) {
  double P_err = 3 * P_err_s + 2 * P_err_i;
  int error0 = W_total * (1.0f - P_err) + 0.5f;
  int error1 = W_total * P_err + 0.5f;
  double homo_base = fisher_homo_base(lc, W[targetChar], error0, W_total);

  int ary[11];
  GenerateIntegerArrayFlanking5(W_total - W[targetChar], ary);
  double homo_base_e = fisher_homo_base_e(lc, ary[errorJitterStart], error1, W_total);
  for (int i = errorJitterStart + 1; i < errorJitterEnd; ++i) {
    double tmp = fisher_homo_base_e(lc, ary[i], error1, W_total);
    homo_base_e = homo_base_e > tmp ? homo_base_e : tmp;
  }
  double ret = homo_base * OVERFLOW_HANDLER * homo_base_e;

  return checkOverflowAndReturn(ret);
}

double Likelihood_homo_indel(unsigned int W[], int targetChar, int W_total,
                             double P_err_s, double P_err_i, LikelihoodCache *lc) {
  double P_err = 4 * P_err_s + P_err_i;
  int error0 = W_total * (1.0f - P_err) + 0.5f;
  int error1 = W_total * P_err + 0.5f;
  double homo_indel = fisher_homo_indel(lc, W[targetChar], error0, W_total);

  int ary[11];
  GenerateIntegerArrayFlanking5(W_total - W[targetChar], ary);
  double homo_indel_e = fisher_homo_indel_e(lc, ary[errorJitterStart], error1, W_total);
  for (int i = errorJitterStart + 1; i < errorJitterEnd; ++i) {
    double tmp = fisher_homo_indel_e(lc, ary[i], error1, W_total);
    homo_indel_e = homo_indel_e > tmp ? homo_indel_e : tmp;
  }
  double ret = homo_indel * OVERFLOW_HANDLER * homo_indel_e;

  return checkOverflowAndReturn(ret);
}

double Likelihood_het_base(unsigned int W[], int targetChar1, int targetChar2, int W_total,
                           double P_err_s, double P_err_i, LikelihoodCache *lc) {
  double P_err = 2 * P_err_s + 2 * P_err_i;
  double P_err_half = P_err * 0.5f;
  int error0 = W_total * (0.5f - P_err_half) + 0.5f;
  int error1 = W_total * P_err_half + 0.5f;
  int W_total_half = W_total * 0.5f + 0.5f;

  double het_base_1 = fisher_het_base(lc, W[targetChar1], error0, W_total_half);
  double het_base_2 = fisher_het_base(lc, W[targetChar2], error0, W_total_half);

  int ary[11];
  GenerateIntegerArrayFlanking5(W_total - W[targetChar1] - W[targetChar2], ary);
  double het_base_e = fisher_het_base_e(lc, ary[errorJitterStart], error1, W_total);
  for (int i = errorJitterStart + 1; i < errorJitterEnd; ++i) {
    double tmp = fisher_het_base_e(lc, ary[i], error1, W_total);
    het_base_e = het_base_e > tmp ? het_base_e : tmp;
  }

  double ret = pow(het_base_1, geometricMeanRoot) * pow(het_base_2, geometricMeanRoot) * OVERFLOW_HANDLER * het_base_e;

  return checkOverflowAndReturn(ret);
}

double Likelihood_het_indel_1(unsigned int W[], int targetChar1, int targetChar2, int W_total,
                              double P_err_s, double P_err_i, LikelihoodCache *lc) {
  double P_err = 3 * P_err_s + P_err_i;
  double P_err_half = P_err * 0.5f;
  int error0 = W_total * (0.5f - P_err_half) + 0.5f;
  int error1 = W_total * P_err_half + 0.5f;
  int W_total_half = W_total * 0.5f + 0.5f;

  double het_indel_1_char1 = fisher_het_indel_1(lc, W[targetChar1], error0, W_total_half);
  double het_indel_1_char2 = fisher_het_indel_2(lc, W[targetChar2], error0, W_total_half);

  int ary[11];
  GenerateIntegerArrayFlanking5(W_total - W[targetChar1] - W[targetChar2], ary);
  double het_indel_e_1 = fisher_het_indel_e_1(lc, ary[errorJitterStart], error1, W_total);
  for (int i = errorJitterStart + 1; i < errorJitterEnd; ++i) {
    double tmp = fisher_het_indel_e_1(lc, ary[i], error1, W_total);
    het_indel_e_1 = het_indel_e_1 > tmp ? het_indel_e_1 : tmp;
  }

  double ret =
      pow(het_indel_1_char1, geometricMeanRoot) * pow(het_indel_1_char2, geometricMeanRoot) * OVERFLOW_HANDLER *
      het_indel_e_1;

  return checkOverflowAndReturn(ret);
}

double Likelihood_het_indel_2(unsigned int W[], int targetChar1, int targetChar2, int W_total, double P_err_s,
                              LikelihoodCache *lc) {
  double P_err = 4 * P_err_s;
  double P_err_half = P_err * 0.5f;
  int error0 = W_total * (0.5f - P_err_half) + 0.5f;
  int error1 = W_total * P_err_half + 0.5f;
  int W_total_half = W_total * 0.5f + 0.5f;

  double het_indel_2_char1 = fisher_het_indel_2(lc, W[targetChar1], error0, W_total_half);
  double het_indel_2_char2 = fisher_het_indel_2(lc, W[targetChar2], error0, W_total_half);

  int ary[11];
  GenerateIntegerArrayFlanking5(W_total - W[targetChar1] - W[targetChar2], ary);
  double het_indel_e_2 = fisher_het_indel_e_2(lc, ary[errorJitterStart], error1, W_total);
  for (int i = errorJitterStart + 1; i < errorJitterEnd; ++i) {
    double tmp = fisher_het_indel_e_2(lc, ary[i], error1, W_total);
    het_indel_e_2 = het_indel_e_2 > tmp ? het_indel_e_2 : tmp;
  }

  double ret =
      pow(het_indel_2_char1, geometricMeanRoot) * pow(het_indel_2_char2, geometricMeanRoot) * OVERFLOW_HANDLER *
      het_indel_e_2;;

  return checkOverflowAndReturn(ret);
}

void computeGenotypeLikelihood(unsigned int W[], char refChar, GenotypeLikelihood *lh,
                               double p_err_s, double p_err_i, LikelihoodCache *lc) {
  unsigned int indelCount = W[4] + W[5];
  int idx = (int) refChar;
  unsigned int refWeightBackup;
  if (indelCount) {
    refWeightBackup = W[idx];
    double tmp = (double) indelCount * (indelDifficultyCoefficient) + 0.5f;
    tmp = W[idx] - tmp;
    if (tmp <= 0) { tmp = 0; }
    W[idx] = tmp;
  }

  int W_total = 0;
  int i, j;
  for (i = 0; i < COUNTER_NUM; ++i) {
    W_total += W[i];
  }
  int index = 0;
  for (i = 0; i < ALPHABET_SIZE; i++) {
    for (j = i; j < ALPHABET_SIZE; j++) {
      if (i == j) {
        lh->likelihood[index] = fmin(1.0, Likelihood_homo_base(W, i, W_total, p_err_s, p_err_i, lc));
      } else {
        lh->likelihood[index] = fmin(1.0, Likelihood_het_base(W, i, j, W_total, p_err_s, p_err_i, lc));
      }
      index++;
    }
  }
#ifdef GENOTYPE_16
  lh->likelihood[index++] = fmin(1.0, Likelihood_het_indel_1(W, 0, 4, W_total, p_err_s, p_err_i, lc));
  lh->likelihood[index++] = fmin(1.0, Likelihood_het_indel_1(W, 1, 4, W_total, p_err_s, p_err_i, lc));
  lh->likelihood[index++] = fmin(1.0, Likelihood_het_indel_1(W, 2, 4, W_total, p_err_s, p_err_i, lc));
  lh->likelihood[index++] = fmin(1.0, Likelihood_het_indel_1(W, 3, 4, W_total, p_err_s, p_err_i, lc));
  lh->likelihood[index++] = fmin(1.0, Likelihood_homo_indel(W, 4, W_total, p_err_s, p_err_i, lc));
  lh->likelihood[index++] = fmin(1.0, Likelihood_het_indel_2(W, 5, 4, W_total, p_err_s, lc));
#endif
  if (indelCount) {
    W[idx] = refWeightBackup;
  }
}

void
computeGenotypeLikelihoodWrapper(InputStat *in, MetaReference *reference, GenotypeLikelihood *lh, IniParams ini_params,
                                 unsigned int input_size, LikelihoodCache *lc) {
  omp_set_num_threads(ini_params.Ini_NumOfCpuThreads);
#pragma omp parallel for
  for (int i = 0; i < input_size; i++) {
    
    computeGenotypeLikelihood(in[i].W, reference[i].refChar, &(lh[i]), ini_params.Ini_BalanceSubError / ALPHABET_SIZE,
                              ini_params.Ini_UnbalanceSubError / 2, lc);
  }
}

#else

#endif

float Strand_bias(int F, int R, int S) {
  return fisher2tail(F, R, S * 0.5, S * 0.5);
}

void computeStrandBiasG(InputStat *in, StrandBias *sb) {
  for (int a = 0; a < ALPHABET_SIZE; ++a) {
    sb->bias[a] = fminf(1.0f, Strand_bias(in->F[a], in->R[a], in->F[a] + in->R[a]));
  }
}

void selectGenotype(unsigned char refChar,
                    GenotypeLikelihood *lh,
                    SnpCallingInfo *info) {
#ifdef GENOTYPE_10
  info->pD = lh->likelihood[0]
    + lh->likelihood[1]
    + lh->likelihood[2]
    + lh->likelihood[3]
    + lh->likelihood[4]
    + lh->likelihood[5]
    + lh->likelihood[6]
    + lh->likelihood[7]
    + lh->likelihood[8]
    + lh->likelihood[9];
#endif
#ifdef GENOTYPE_16
  info->pD = lh->likelihood[0]
             + lh->likelihood[1]
             + lh->likelihood[2]
             + lh->likelihood[3]
             + lh->likelihood[4]
             + lh->likelihood[5]
             + lh->likelihood[6]
             + lh->likelihood[7]
             + lh->likelihood[8]
             + lh->likelihood[9]
             + lh->likelihood[10]
             + lh->likelihood[11]
             + lh->likelihood[12]
             + lh->likelihood[13]
             + lh->likelihood[14]
             + lh->likelihood[15];
#endif

  if (info->pD == 0.0) {
    const int genotypeMap[][ALPHABET_SIZE] = {{0, 1, 2, 3},
                                              {1, 4, 5, 6},
                                              {2, 5, 7, 8},
                                              {3, 6, 8, 9}};
    info->genotype = genotypeMap[refChar][refChar];
    for (int i = 0; i < 16; ++i) { info->flh[i] = 0; }
    info->bestD = 0;
    info->secondD = 0;
    return;
  }

  info->flh[0] = lh->likelihood[0] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[0])];
  info->flh[1] = lh->likelihood[1] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[1])];
  info->flh[2] = lh->likelihood[2] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[2])];
  info->flh[3] = lh->likelihood[3] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[3])];
  info->flh[4] = lh->likelihood[4] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[4])];
  info->flh[5] = lh->likelihood[5] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[5])];
  info->flh[6] = lh->likelihood[6] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[6])];
  info->flh[7] = lh->likelihood[7] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[7])];
  info->flh[8] = lh->likelihood[8] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[8])];
  info->flh[9] = lh->likelihood[9] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[9])];
#ifdef GENOTYPE_16
  info->flh[10] = lh->likelihood[10] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[10])];
  info->flh[11] = lh->likelihood[11] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[11])];
  info->flh[12] = lh->likelihood[12] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[12])];
  info->flh[13] = lh->likelihood[13] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[13])];
  info->flh[14] = lh->likelihood[14] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[14])];
  info->flh[15] = lh->likelihood[15] / info->pD * gpuSnpPrior[((refChar) << 6) | (gpuGenotype2int[15])];
#endif
  register int bestI;
  register double flh0, flh1;
  if (info->flh[0] >= info->flh[1]) {
    bestI = 0;
    flh0 = info->flh[0];
    flh1 = info->flh[1];
  } else {
    bestI = 1;
    flh0 = info->flh[1];
    flh1 = info->flh[0];
  }

  if (info->flh[2] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[2];
    bestI = 2;
  } else if (info->flh[2] > flh1) {
    flh1 = info->flh[2];
  }

  if (info->flh[3] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[3];
    bestI = 3;
  } else if (info->flh[3] > flh1) {
    flh1 = info->flh[3];
  }

  if (info->flh[4] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[4];
    bestI = 4;
  } else if (info->flh[4] > flh1) {
    flh1 = info->flh[4];
  }

  if (info->flh[5] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[5];
    bestI = 5;
  } else if (info->flh[5] > flh1) {
    flh1 = info->flh[5];
  }

  if (info->flh[6] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[6];
    bestI = 6;
  } else if (info->flh[6] > flh1) {
    flh1 = info->flh[6];
  }

  if (info->flh[7] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[7];
    bestI = 7;
  } else if (info->flh[7] > flh1) {
    flh1 = info->flh[7];
  }

  if (info->flh[8] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[8];
    bestI = 8;
  } else if (info->flh[8] > flh1) {
    flh1 = info->flh[8];
  }

  if (info->flh[9] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[9];
    bestI = 9;
  } else if (info->flh[9] > flh1) {
    flh1 = info->flh[9];
  }
#ifdef GENOTYPE_16
  if (info->flh[10] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[10];
    bestI = 10;
  } else if (info->flh[10] > flh1) {
    flh1 = info->flh[10];
  }

  if (info->flh[11] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[11];
    bestI = 11;
  } else if (info->flh[11] > flh1) {
    flh1 = info->flh[11];
  }

  if (info->flh[12] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[12];
    bestI = 12;
  } else if (info->flh[12] > flh1) {
    flh1 = info->flh[12];
  }

  if (info->flh[13] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[13];
    bestI = 13;
  } else if (info->flh[13] > flh1) {
    flh1 = info->flh[13];
  }

  if (info->flh[14] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[14];
    bestI = 14;
  } else if (info->flh[14] > flh1) {
    flh1 = info->flh[14];
  }

  if (info->flh[15] > flh0) {
    flh1 = flh0;
    flh0 = info->flh[15];
    bestI = 15;
  } else if (info->flh[15] > flh1) {
    flh1 = info->flh[15];
  }
#endif

  info->bestD = flh0;
  info->secondD = flh1;
  info->genotype = bestI;
}

void computeGenotypeWrapper(unsigned char *refChar,
                            GenotypeLikelihood *lh, SnpCallingInfo *info, IniParams ini_params,
                            unsigned int input_size) {
  omp_set_num_threads(ini_params.Ini_NumOfCpuThreads);
#pragma omp parallel for
  for (int i = 0; i < input_size; ++i) {
    selectGenotype(refChar[i], &(lh[i]), &(info[i]));
    
    
  }
}

void computeStrandBias(InputStat *in, StrandBias *sb, IniParams ini_params, unsigned int input_size) {
  omp_set_num_threads(ini_params.Ini_NumOfCpuThreads);
#pragma omp parallel for
  for (int i = 0; i < input_size; ++i) {
    computeStrandBiasG(&(in[i]), &(sb[i]));
  }
}

