#include "header.h"

int main(){
    int n=500;
    int N=n+1;
    int nt=400;
    int Nt=nt+1;
    double L=5.0;
    double tmax=5.0;
    double dt=tmax/nt;
    double dx=L/n;
    double hbar=1.0;
    double m=1.0;
    double E=1.0;
    double x0=L/5;
    double k0=sqrt(2.0*m*E)/hbar;
    double sigmak=0.5;
    double *x=malloc(N*sizeof(double));
    for (int i=0;i<=n;i++) x[i]=i*dx;
    double *t=malloc(Nt*sizeof(double));
    for (int p=0;p<=nt;p++) t[p]=dt*p;
    cmp *psi=calloc(N*Nt,sizeof(cmp));
    cmp *Gamma=calloc(N,sizeof(cmp));
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    free(psi);
    free(x);
    free(Gamma);
    free(t);
    return 0;
}
