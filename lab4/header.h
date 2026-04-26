#define M_PI 3.141592653589793115997963468544185161590576171875
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

typedef double complex cmp;

void initPsi(cmp *psi, int n, int nt, double sigmak, double *x, double x0, double k0){
    double sigma0=0.5*sigmak;
    for (int i=0;i<=n;i++){
        psi[i*(nt+1)]=pow(1.0/(2*sigma0*sigma0*M_PI),0.25)*exp(-pow(x[i]-x0,2)/(4.0*sigma0*sigma0))*cexp(I*k0*(x[i]-x0));
    }
    psi[0]=0.0;
    psi[n]=0.0;
}

void thomas(cmp *a, cmp *b, cmp *c, cmp *d, int n, cmp *x){
    cmp *cprim = calloc(n,sizeof(cmp));
    cmp *dprim = calloc(n,sizeof(cmp));
    cprim[0]=c[0]/b[0];
    for (int i=1;i<=n-2;i++){
        cprim[i]=c[i]/(b[i]-a[i]*cprim[i-1]);
    }
    dprim[0]=d[0]/b[0];
    for (int i=1;i<=n-1;i++){
        dprim[i]=(d[i]-a[i]*dprim[i-1])/(b[i]-a[i]*cprim[i]);
    }
    x[n-1]=dprim[n-1];
    for (int i=n-2; i>=0;i--){
        x[i]=dprim[i]-cprim[i]*x[i+1];
    }


    free(cprim);
    free(dprim);
}

void trimul(cmp *in, cmp *out, cmp *tri, int n){
    for (int i=0;i<n;i++){
        out[i]=0.0;
        for (int k=(i!=0?(i-1):0); k<=(i!=(n-1)?(i+1):(n-1)); k++){
            out[i] += tri[i*n+k]*in[k];
        }
    }
}

void initH(cmp *H, int n, cmp *Gamma, double hbar, double m, double dx){
    double coeff=hbar*hbar/(m*dx*dx);
    for (int i=0;i<n;i++){
        H[i*n+i]=coeff + Gamma[i];
        if (i!=0){ 
            H[(i-1)*n+i]=-0.5*coeff; // nad
            H[i*n+i-1]=-0.5*coeff; // pod
        }
    }
}

void initA(cmp *Aa, cmp* Ab, int n, double hbar, double m, double dx, double dt, cmp *Gamma){
    double coeff=hbar*hbar/(m*dx*dx);
    for (int i=0;i<n;i++){
        Aa[i]=0.5*I*dt*(-0.5*coeff);
        Ab[i]=1.0 + 0.5*I*dt*(coeff + Gamma[i+1]);
    }
}

void initB(cmp *B, int n, double hbar, double m, double dx, double dt, cmp *Gamma){
    double coeff=hbar*hbar/(m*dx*dx);
    for (int i=0;i<n;i++){
        B[i*n+i]=1.0 - 0.5*I*dt*(coeff + Gamma[i+1]);
        if (i!=0){ 
            B[(i-1)*n+i]=-0.5*I*dt*(-0.5*coeff); // nad
            B[i*n+i-1]=-0.5*I*dt*(-0.5*coeff); // pod
        }
    }
}

void initGamma(cmp *Gamma, int n){
    for (int i=0;i<n;i++){
        Gamma[i] = 0.0;
    }
}

void solve(cmp *psi, cmp *Gamma, int N, double hbar, double m ,double dx, double dt, int Nt){
    cmp *Aa=calloc((N-2),sizeof(cmp));
    cmp *Ab=calloc(N-2,sizeof(cmp));
    initA(Aa,Ab,N-2,hbar,m,dx,dt,Gamma);
    cmp *B=calloc((N-2)*(N-2),sizeof(cmp));
    initB(B,N-2,hbar,m,dx,dt,Gamma);
    cmp *c = calloc((N-2),sizeof(cmp));
    cmp *psi_uciete=calloc((N-2),sizeof(cmp));
    for (int p=0; p<Nt-1;p++){
        for (int i=0;i<N-2;i++){
            psi_uciete[i]=psi[(i+1)*Nt+p]; // kopia robocza wektroa psi
        }
        trimul(psi_uciete,c,B,N-2); // rozwiazanie na c
        thomas(Aa,Ab,Aa,c,N-2,psi_uciete); // roziwazanie na psi^{p+1}
        for (int i=0;i<N-2;i++){
            psi[(i+1)*Nt+p+1]=psi_uciete[i]; // zapisanie na nowe psi
        }
    }   

    free(Aa);
    free(Ab);
    free(B);
    free(c);
    free(psi_uciete);
}
