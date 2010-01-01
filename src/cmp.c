/*******************************************************************/
/* Computation of the Conway Maxwell Poisson                       */
/*******************************************************************/

#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

double cmp(int x, double llambda, double nu, double lzcmp, int give_log)
  {
     double dev;
     dev=x*llambda-nu*lgamma(x+1.0)-lzcmp;
     if(give_log){
      return( dev );
     }else{
      return( exp(dev) );
     }
  }

double zcmp(double lambda, double nu, double err, int give_log)
  {
     double mj, z, aj;
     int j;
     z = 1.0 + lambda;
     aj = lambda;
     for (j = 2; j < 100; j++){
       mj=lambda/pow((double)j,nu);
       aj*=mj;
       z+=aj;
     }
     if(give_log){
      return( log(z) );
     }else{
      return( z );
     }
  }

void dcmp (int *x, double *lambda, double *nu, int *n, double *err, int *give_log, double *val) {
  int i, give_log1=1;
  double lzcmp;
  lzcmp = zcmp(*lambda, *nu, *err, give_log1);
  for (i = 0; i < *n; i++){
    val[i] = cmp(x[i], log(*lambda), *nu, lzcmp, *give_log);
  }
}

void rcmp (int *x, double *lambda, double *nu, int *n, int *K, double *err) {
  int i, Ki, ni, popi, give_log0=0, give_log1=1;
  double gb, lzcmp, llambda;
  double *pi = (double *) malloc(sizeof(double) * (*K));
  lzcmp = zcmp(*lambda, *nu, *err, give_log1);
  llambda = log(*lambda);
  Ki = (*K);
  ni = (*n);

  GetRNGstate();  /* R function enabling uniform RNG */

  for (i = 0; i < Ki; i++){
    pi[i] = cmp(i, llambda, *nu, lzcmp, give_log0);
  }
  for (i=1; i< Ki; i++){
    pi[i]=pi[i-1]+pi[i];
  }
  for (i = 0; i < ni; i++){
   gb = pi[Ki-1] * unif_rand();
   popi = 0;
   while(gb > pi[popi]){popi++;}
   x[i] = popi;
  }

  PutRNGstate();  /* Disable RNG before returning */

  free(pi);
}
