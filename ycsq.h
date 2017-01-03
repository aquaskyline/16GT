#ifndef YCSQ_H
#define YCSQ_H




#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double chi2UniformDistance(const double a, const double b, const double c, const double d);

double chisqr(int Dof, double Cv);


#endif
