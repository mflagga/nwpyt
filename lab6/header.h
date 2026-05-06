#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>

typedef complex double cmp;

void initV(double *V, int nx, int ny, double *x, double *y, double y0, double x0, double V0, double alpha){
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            V[i*(ny+1)+j]=pow(y[j]-y0,2)*V0/(pow(cosh(alpha*(x[i]-x0)),2));
        }
    }
}

void diagYslice(size_t i, double *V, size_t ny, double *ev, double ty){
    size_t n=ny-1;
    gsl_matrix *H=gsl_matrix_calloc(n,n);
    for (size_t j=0;j<n;j++) gsl_matrix_set(H,j,j,2.0*ty+V[i*(ny+1)+j]);
    for (size_t j=0;j<n-1;j++){
        gsl_matrix_set(H,j+1,j,-ty);
        gsl_matrix_set(H,j,j+1,-ty);
    }
    gsl_vector *eval=gsl_vector_alloc(n);
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(n);
    gsl_eigen_symm(H,eval,w);
    gsl_eigen_symm_free(w);
    gsl_sort_vector(eval);
    for (size_t j=0;j<n;j++){
        ev[j]=gsl_vector_get(eval,j);
    }
    gsl_matrix_free(H);
    gsl_vector_free(eval);
}

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

void solve(cmp *psi, double *V, int N, double t, double E, double k, double a){
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

double TodE(double E, cmp *psi, double *V, int N, double t, double a, double m, double hbar){
    double k=sqrt(2.0*m*E)/hbar;
    solve(psi,V,N,t,E,k,a);
    return pow(cabs(psi[N-1]),2);
}
