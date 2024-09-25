#include "zeta.h"

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>

void zetaC (double *v, double *z) {
  double *b2;
  int a, k;
  int m, n, m2;
  double p, a2, fred, s, sum;
  s = (*v);
  b2 = (double *)R_Calloc( 12 * sizeof(double), double);
  b2[0] = 1.0 / 6.0;
  b2[1] = -1.0 / 30.0;
  b2[2] = 1.0 / 42.0;
  b2[3] = -1.0 / 30.0;
  b2[4] = 5.0 / 66.0;
  b2[5] = -691.0 / 2730.0;
  b2[6] = 7.0 / 6.0;
  b2[7] = -3617.0 / 510.0;
  b2[8] = 4386.7 / 79.8;
  b2[9] = -1746.11 / 3.30;
  b2[10] = 8545.13 / 1.38;
  b2[11] = -2363.64091 / 0.02730;
  a = 12;
  k = 8;
  a2 = a * a;
  p = s / 2.0 / a2;
  sum = 1.0 / (s - 1.0) + 0.5 / a + b2[0] * p;
  for (m = 1; m < k; m++){
   m2 = m + m + 2;
   p=p*(s+m2-3.0)*(s+m2-2.0)/(m2-1.0)/m2/a2;
   sum = sum + p * b2[m];
  }
  fred = exp((s - 1.0) * log(1.0 * a));
  sum = 1.0 + sum / fred;
  for (n = 2; n < a; n++){
   sum = sum + 1.0 * exp(-s * log(1.0 * n));
  }
  *z = sum;
  R_Free(b2);
}
