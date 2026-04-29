#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef double complex cmp;

void thomas(cmp *b, cmp *ac, int n, cmp *d, cmp *x){
    cmp *cprim = calloc(n,sizeof(cmp));
    cmp *dprim = calloc(n,sizeof(cmp));
    cprim[0]=ac[0]/b[0];
    for (int i=1;i<=n-2;i++){
        cprim[i]=ac[i]/(b[i]-ac[i]*cprim[i-1]);
    }
    dprim[0]=d[0]/b[0];
    for (int i=1;i<n;i++){
        dprim[i]=(d[i]-ac[i]*dprim[i-1])/(b[i]-ac[i]*cprim[i-1]);
    }
    x[n-1]=dprim[n-1];
    for (int i=n-2;i>=0;i--){
        x[i]=dprim[i]-cprim[i]*x[i+1];
    }
    free(cprim);
    free(dprim);
}

void solve(cmp *psi, cmp *V, int N, double t, double E, double k, double a){
    cmp *b = calloc(N,sizeof(cmp));
    cmp *ac = calloc(N,sizeof(cmp));
    for (int i=0;i<N;i++){
        b[i]=2*t+V[i]-E;
    }
    b[0] += -t*cexp(I*k*a);
    b[N-1] += -t*cexp(I*k*a);
    for (int i=0;i<N;i++) ac[i]=-t;
    cmp *s = calloc(N,sizeof(cmp));
    s[0]=t*(1.0-cexp(2.0*I*k*a));
    thomas(b,ac,N,s,psi);
    free(s);
    free(ac);
    free(b);
}

void initBarrier(cmp *V, double *x, int n, double start, double szerokosc, cmp V0){
    for (int i=0;i<=n;i++){
        if (x[i] >= start && x[i] <= start + szerokosc) V[i] = V0;
        else V[i]=0.0;
    }
}
