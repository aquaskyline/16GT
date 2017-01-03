#include "fisher.h"




static double lngamm(int z)




{
  double x = 0;
  x += 0.1659470187408462e-06 / (z + 7);
  x += 0.9934937113930748e-05 / (z + 6);
  x -= 0.1385710331296526 / (z + 5);
  x += 12.50734324009056 / (z + 4);
  x -= 176.6150291498386 / (z + 3);
  x += 771.3234287757674 / (z + 2);
  x -= 1259.139216722289 / (z + 1);
  x += 676.5203681218835 / (z);
  x += 0.9999999999995183;
  return (log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5));
}



static double lnfact(int n) {
  if (n <= 1) return (0);
  return (lngamm(n + 1));
}


double lnbico(int n, int k) {
  return (lnfact(n) - lnfact(k) - lnfact(n - k));
}


double hyper_323(int n11, int n1_, int n_1, int n) {
  return (exp(lnbico(n1_, n11) + lnbico(n - n1_, n_1 - n11) - lnbico(n, n_1)));
}

void fisher(int a, int b, int c, int d, double *p_left, double *p_right, double *p_twotail) {
  if (a < 0) a *= -1;
  if (b < 0) b *= -1;
  if (c < 0) c *= -1;
  if (d < 0) d *= -1;
  if (!(a || b || c || d)) {
    *p_left = 1.0;
    *p_right = 1.0;
    *p_twotail = 1.0;
  }

  int n11, n1_, n_1, n;
  n11 = a;            
  n1_ = a + b;        
  n_1 = a + c;        
  n = a + b + c + d;  

  int min, max;       
  min = (0 > n1_ + n_1 - n) ? 0 : n1_ + n_1 - n;
  max = (n1_ < n_1) ? n1_ : n_1;

  double poa; 
  poa = hyper_323(n11, n1_, n_1, n);
  *p_twotail = 0;
  *p_left = 0;
  *p_right = 0;

  int x;
  double pox; 
  int s = n11;
  double pos = poa; 

  for (x = n11 - 1; x >= min; x--) {
    if (x == s - 1) { 

      pox = pos * (double) (s * (n + s - n1_ - n_1)) / ((n1_ - x) * (n_1 - x));
    } else
      pox = hyper_323(x, n1_, n_1, n);

    s = x;
    pos = pox;

    *p_left += pox;
    if (pox < poa * 1.000000001) *p_twotail += pox;
  }

  s = n11;
  pos = poa;
  for (x = n11 + 1; x <= max; x++) {
    if (x == s + 1) { 
      pox = pos * (double) (n1_ - s) * (n_1 - s) / ((n + x - n1_ - n_1) * x);
    } else
      pox = hyper_323(x, n1_, n_1, n);
    s = x;
    pos = pox;
    *p_right += pox;
    if (pox < poa * 1.000000001) *p_twotail += pox;
  }

  *p_left += poa;
  *p_right += poa;
  *p_twotail += poa;
}

double fisher2(int a, int b, int c, int d) {
  double p_left = 0;
  double p_right = 0;
  double p_twotail = 0;
  fisher(a, b, c, d, &p_left, &p_right, &p_twotail);
  return p_twotail;
}
