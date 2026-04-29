#include "header.h"

int main(){
    // parametry jednostkowe
    double m=1.0;
    double hbar=1.0;

    // parametry ukladu
    double L=50.0;
    int n=500;
    cmp V0=0.02;
    double E=0.04;
    int nE=300;
    double LE=2*V0;

    // parametry wtórne
    int N=n+1;
    double a=L/n;
    double t=hbar*hbar/(2*m*a*a);
    double k=sqrtf(2.0*m*E)/hbar;
    double dE=LE/nE;

    // alokacja / inicjalizacja
    double *x = malloc(N*sizeof(double));
    for (int i=0;i<=n;i++) x[i] = i*a;
    cmp *psi = calloc(N,sizeof(cmp));
    cmp *V = calloc(N,sizeof(cmp));
    initBarrier(V,x,n,L/2,L/10,V0);

    // rozwiazanie
    solve(psi,V,N,t,E,k,a);

    // zapisanie
    FILE *p1psi=fopen("p1psi.csv","w");
    for (int i=0;i<=n;i++){
        fprintf(p1psi,"%lf,%lf\n",x[i],pow(cabs(psi[i]),2));
    }

    FILE *p1E=fopen("p1E.csv","w");
    double T;
    for (int i=0;i<nE;i++){
        E=(i+1)*dE;
        T=TodE(E,psi,V,N,t,a,m,hbar);
        fprintf(p1E,"%lf,%lf\n",E,T);
    }

    // czystki
    free(psi);
    free(V);
    free(x);
    fclose(p1psi);
    fclose(p1E);

    // return zero
    return 0;
}
