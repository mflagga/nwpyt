#include "header.h"

int main(){
    // parametry
    int n=500;
    int N=n+1;
    int nt=40;
    int Nt=nt+1;
    double L=10.0;
    double tmax=5.0;
    double dt=tmax/nt;
    double dx=L/n;
    double hbar=1.0;
    double m=1.0;
    double E=1.0;
    double x0=L/5;
    double k0=sqrt(2.0*m*E)/hbar;
    double sigmak=0.2;
    // alokacje inicjalizacje
    double *x=malloc(N*sizeof(double));
    for (int i=0;i<=n;i++) x[i]=i*dx;
    double *t=malloc(Nt*sizeof(double));
    for (int p=0;p<=nt;p++) t[p]=dt*p;
    cmp *psi=calloc(N*Nt,sizeof(cmp));
    cmp *Gamma=calloc(N,sizeof(cmp));

    // punkt pierwszy
    FILE *p1misc=fopen("p1misc.csv","w");
    fprintf(p1misc,"%d,%d,%lf,%lf,",n,nt,L,sigmak);
    initGamma(Gamma, N);
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    FILE *p1 = fopen("p1.csv","w");
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p++){
            if (p!=0 || i!=0) fprintf(p1,",");
            fprintf(p1,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }
    fprintf(p1,"\n");
    sigmak*=1.2;
    fprintf(p1misc,"%lf,",sigmak);
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p++){
            if (p!=0 || i!=0) fprintf(p1,",");
            fprintf(p1,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }
    fprintf(p1,"\n");
    sigmak*=1.2;
    fprintf(p1misc,"%lf",sigmak);
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p++){
            if (p!=0 || i!=0) fprintf(p1,",");
            fprintf(p1,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }


    // czystki
    free(psi);
    free(x);
    free(Gamma);
    free(t);
    fclose(p1);
    fclose(p1misc);
    // return zero
    return 0;
}
