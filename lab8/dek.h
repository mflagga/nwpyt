#ifndef DEK_H
#define DEK_H

#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef complex double cmp;

void matmul(cmp *A, cmp *B, cmp *C, int N);
void dagger(cmp *A, cmp *Adagger, int N);
void initH(cmp *H, int N, double t, double eps);
void zapiszModulMacierzy(cmp *M, int N, char *nazwa);
void zapiszReMacierzy(cmp *M, int N, char *nazwa);
void zapiszImMacierzy(cmp *M, int N, char *nazwa);

#endif
