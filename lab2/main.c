#include "header.h"

int main(){
    // parametry
    double L = 5;
    int N = 500;
    double dx = L/N;
    double delta = 1e-7;
    double m=0.0;
    double hbar=0.0;
    double dtau=0.5*m*dx*dx/hbar;
    // alokacja/inicjalizacja
    double *x = malloc((N+1)*sizeof(double));
    double *V = malloc((N+1)*sizeof(double));
    for (int i=0;i<=N;i++){
        x[i] = i*dx;
        V[i] = 0.0;
    }
    double *psi = malloc((N+1)*sizeof(double));
    double E=met_czas_urojony(L,N,V,delta,psi,x,hbar,m,dx,dtau);
    double *psi_analityczne = malloc((N+1)*sizeof(double));    

    // czystki
    free(x);
    free(V);
    free(psi);
    free(psi_analityczne);
    // return zero
    return 0;
}
