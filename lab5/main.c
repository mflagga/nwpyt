#include "header.h"

int main(){
    // parametry jednostkowe
    double m=1.0;
    double hbar=1.0;

    // parametry ukladu
    double L=50.0;
    int n=50;
    double E=0.04;

    // parametry wtórne
    int N=n+1;
    double a=L/n;
    double t=hbar*hbar/(2*m*a*a);
    double k=sqrtf(2.0*m*E)/hbar;

    // alokacja / inicjalizacja
    double *x = malloc(N*sizeof(double));
    for (int i=0;i<=n;i++) x[i] = i*a;
    cmp *psi = calloc(N,sizeof(cmp));
    cmp *V = calloc(N,sizeof(cmp));
    initBarrier(V,x,n,L/2,L/10,0.02);

    // rozwiazanie
    solve(psi,V,N,t,E,k,a);

    // zapisanie
    FILE *psifile=fopen("psifile.csv","w");
    for (int i=0;i<=n;i++){
        if (i!=0) fprintf(psifile,",");
        fprintf(psifile,"%lf",pow(cabs(psi[i]),2));
    }

    // czystki
    free(psi);
    free(V);
    fclose(psifile);

    // return zero
    return 0;
}
