#include "ycsq.h"

double chi2UniformDistance(const double a, const double b, const double c, const double d) {
  
  const double ab = a + b;
  const double cd = c + d;
  const double ac = a + c;
  const double bd = b + d;
  const double t = a + b + c + d;
  const double ea = ab * ac / t;
  const double eb = ab * bd / t;
  const double ec = ac * cd / t;
  const double ed = cd * bd / t;
  double dist = 0.;
  if (ea > 0) dist += (pow((abs(a - ea) - 0.5), 2) / ea);
  if (eb > 0) dist += (pow((abs(b - eb) - 0.5), 2) / eb);
  if (ec > 0) dist += (pow((abs(c - ec) - 0.5), 2) / ec);
  if (ed > 0) dist += (pow((abs(d - ed) - 0.5), 2) / ed);
  return dist;
}

double chisqr(int Dof, double Cv) {
  if (Cv < 0 || Dof < 1) {
    return 0.0;
  }

  double K = ((double) Dof) * 0.5;


  double PValue = pow(0.5, K) / tgamma(K) * pow(Cv, K - 1) * exp((-0.5) * Cv);

  return PValue > 1.0 ? 1.0 : PValue;
}

