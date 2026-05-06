#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef complex double cmp;

void initV(double *V, int nx, int ny, double *x, double *y, double y0, double x0, double V0, double alpha){
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            V[i*(ny+1)+j]=pow(y[j]-y0,2)*V0/(pow(cosh(alpha*(x[i]-x0)),2));
        }
    }
}

void diagYslice(int i){
    
}
