#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef double complex cmp;

void initPsi(cmp *psi, int n, int nt, double sigmak, double *x, double x0, double k0){
    double sigma0=0.5*sigmak;
    for (int i=0;i<=n;i++){
        psi[i*(nt+1)]=pow(1.0/(2*sigma0*sigma0*M_PI),0.25)*exp(-pow(x[i]-x0,2)/(4.0*sigma0*sigma0))*cexp(I*k0*(x[i]-x0));
    }
}

void thomas(cmp *a, cmp *b, cmp *c, cmp *d, int n, cmp *x){
    cmp *cprim=calloc((n+1),sizeof(cmp));
    cmp *dprim=calloc((n+1),sizeof(cmp));
    cprim[1]=c[1]/b[1];
    for (int i=2;i<n;i++) cprim[i]=c[i]/(b[i]-a[i]*cprim[i-1]);
    dprim[1]=d[1]/b[1];
    for (int i=2;i<=n;i++) dprim[i]=(d[i]-a[i]*dprim[i-1])/(b[i]-a[i]*cprim[i-1]);
    x[n-1]=dprim[n];
    for (int i=n-1; i>=1; i--) x[i-1]=dprim[i]-cprim[i]*x[i];
    free(cprim);
    free(dprim);
}
