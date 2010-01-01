#ifndef CMP_H
#define CMP_H

void dcmp (int *x, double *lambda, double *nu, int *n, double *err, int *give_log, double *val);
double cmp(int x, double llambda, double nu, double lzcmp, int give_log);
double zcmp(double lambda, double nu, double err, int give_log);
void rcmp (int *x, double *lambda, double *nu, int *n, int *K, double *err);
#endif /* CMP_H */
