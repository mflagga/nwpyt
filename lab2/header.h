#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

double f(double x){
    return 0.4+0.0*x;
}

double psi0odx(double x, double L){
    return x*(L-x)*f(x);
}

void normuj_ffalowa(double *psi, int N, double dx){
    double norm=0.0;
    for (int i=0;i<=N;i++){
        norm += pow(fabs(psi[i]),2)*dx;
    }
    double c = 1.0/sqrt(norm);
    for (int i=0;i<=N;i++) psi[i] *= c; 
    psi[0]=0.0; psi[N]=0.0;
}

double licz_E(double *psi, double *V, double hbar, double m, double dx, int N){
    double E=0.0;
    double c=-hbar*hbar/(2.0*m*dx*dx);
    for (int i=1;i<N;i++){
        E += psi[i]*(c*(psi[i-1]-2*psi[i]+psi[i+1])+V[i]*psi[i])*dx;
    }
    return E;
}

double met_czas_urojony(double L, int N, double *V, double delta, double *psi, double *x, double hbar, double m, double dx, double dtau){
    int p=0;
    for (int i=0;i<=N;i++) psi[i] = psi0odx(x[i],L);
    double *psin = malloc((N+1)*sizeof(double));
    double E;
    double En=licz_E(psi,V,hbar,m,dx,N);
    double c=hbar*dtau/(2*m*dx*dx);
    do{
        E=En;
        for (int i=1;i<N;i++){
            psin[i]=c*(psi[i-1]-2*psi[i]+psi[i+1])-psi[i]*(V[i]*dtau/hbar-1.0);
        }
        normuj_ffalowa(psin,N,dx);
        En=licz_E(psin,V,hbar,m,dx,N);
        memcpy(psi,psin,(N+1)*sizeof(double));
        p++;
    }while(fabs(En-E)>=delta);
    printf("Zbiegło się w %d iteracjach\n",p);
    free(psin);
    return E;
}
