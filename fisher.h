
#ifndef FISHER_H
#define FISHER_H




#include <stdlib.h>
#include <math.h>
#include <float.h>


void fisher(int a, int b, int c, int d, double *p_left, double *p_right, double *p_twotail);


double fisher2(int a, int b, int c, int d);

double lnbico(int n, int k);

double hyper_323(int n11, int n1_, int n_1, int n);

#endif

