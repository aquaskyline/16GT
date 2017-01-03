#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "SNPFunctions.h"

#define snpAltHomProb 0.0005
#define snpRefHetProb 0.001
#define indelAltHomProb 0.00005
#define indelRefHetProb 0.0001
#define titvRatio 2.1

double *createSnpPriorArray() {
  double *ary = (double *) malloc(256 * sizeof(double));
  memset(ary, '\0', sizeof(double) * 256);

  char t_base, allele1, allele2;
  for (t_base = 0; t_base < 4; ++t_base) {
    for (allele1 = 0; allele1 < 6; ++allele1) {
      for (allele2 = 0; allele2 < 6; ++allele2) {
        if (allele1 == t_base && allele2 == t_base) {
          
          ary[t_base << 6 | allele1 << 3 | allele2] = 1;
        } else if (allele1 == t_base || allele2 == t_base) {
          
          if (allele1 > 3 || allele2 > 3) 
          {
            ary[t_base << 6 | allele1 << 3 | allele2] = indelRefHetProb;
          } else 
          {
            ary[t_base << 6 | allele1 << 3 | allele2] = snpRefHetProb;
          }
        } else if (allele1 == allele2) {
          
          if (allele1 > 3) 
          {
            ary[t_base << 6 | allele1 << 3 | allele2] = indelAltHomProb;
          } else 
          {
            ary[t_base << 6 | allele1 << 3 | allele2] = snpAltHomProb;
          }
        } else {
          if (allele1 > 3 && allele2 > 3) 
          {
            ary[t_base << 6 | allele1 << 3 | allele2] = indelAltHomProb * indelRefHetProb;
          } else if (allele1 > 3 || allele2 > 3) 
          {
            ary[t_base << 6 | allele1 << 3 | allele2] = snpAltHomProb * indelRefHetProb;
          } else 
          {
            
            ary[t_base << 6 | allele1 << 3 | allele2] = snpAltHomProb * snpRefHetProb;
          }
        }
        if (abs(allele1 - t_base) == 2) {
          
          ary[t_base << 6 | allele1 << 3 | allele2] *= titvRatio;
        }
        if (abs(allele2 - t_base) == 2) {
          
          ary[t_base << 6 | allele1 << 3 | allele2] *= titvRatio;
        }
      }
    }
  }
  


  return ary;
}

void freeSnpPriorArray(double *ary) {
  free(ary);
}





float fisher2tailC(int a, int b, int c, int d);

__inline__ static float lnfactC(int n) {
  return lgamma((float) (n + 1));
}

__inline__ float lnbicoC(int n, int k) {
  return (lnfactC(n) - lnfactC(k) - lnfactC(n - k));
}

__inline__ float hyper_323C(int n11, int n1_, int n_1, int n) {
  return (exp(lnbicoC(n1_, n11) + lnbicoC(n - n1_, n_1 - n11) - lnbicoC(n, n_1)));
}

__inline__ float fisher2tailC(int n11, int n1_, int n_1, int n) 
{
  
  
  
  
  
  n += n11 + n1_ + n_1; 
  n1_ += n11; 
  n_1 += n11;

  float poa; 
  poa = hyper_323C(n11, n1_, n_1, n); 

  int x;
  float pox = poa; 
  float p_twotail = 0;
  p_twotail += poa;
  float poa_e = poa * 1.000000001f;

  int min = 0;
  if (n1_ + n_1 - n > 0) {
    min = n1_ + n_1 - n;
  }

  for (x = n11; x > min; --x) {
    pox = pox * (float) (x * (n + x - n1_ - n_1)) / ((n1_ - x + 1) * (n_1 - x + 1));
    p_twotail += pox * (pox < poa_e);
  }

  int max = n_1;
  if (n1_ < n_1) {
    max = n1_;
  }
  pox = poa;
  for (x = n11 + 1; x <= max; ++x) {
    pox = pox * (float) (n1_ - (x - 1)) * (n_1 - (x - 1)) / ((n + x - n1_ - n_1) * x);
    p_twotail += pox * (pox < poa_e);
  }

  return p_twotail;
}

float getStrandBiasC(int F, int R, int S) {
  return fisher2tailC(F, R, S * 0.5, S * 0.5);
}

float getBaseQualityBias(int W, int S) {
  return fisher2tailC(W, S, S, S);
}

void computeStrandBias(MetaSnpCounter *snpCounter, StrandBias *strandBias) {
  unsigned int i;
  for (i = 0; i < ALPHABET_SIZE; ++i) {
    strandBias->bias[i] = fmin(1.0f,
                               getStrandBiasC(snpCounter->F[i], snpCounter->R[i], snpCounter->F[i] + snpCounter->R[i]));
  }
}

void computeBaseQualityBias(MetaSnpCounter *snpCounter, BaseQualityBias *baseQualityBias) {
  unsigned int i;
  for (i = 0; i < ALPHABET_SIZE; ++i) {
    baseQualityBias->bias[i] = fmin(1.0f,
                                    getBaseQualityBias(snpCounter->W[i], (snpCounter->F[i] + snpCounter->R[i]) * 4));
  }
}

