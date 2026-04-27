#include "header.h"

int main(){
    // parametry
    int n=500;
    int N=n+1;
    int nt=400;
    int Nt=nt+1;
    double L=20.0;
    double tmax=1.0;
    double dt=tmax/nt;
    double dx=L/n;
    double hbar=1.0;
    double m=1.0;
    double E=40.0;
    double x0=5.0;
    double k0=sqrt(2.0*m*E)/hbar;
    double sigmak=0.2;
    // double a=1.0;
    double V0=50.0;
    int frames = 75;
    int fps=frames/5;
    int co_ktora=nt/frames;
    // alokacje inicjalizacje
    double *x=malloc(N*sizeof(double));
    for (int i=0;i<=n;i++) x[i]=i*dx;
    double *t=malloc(Nt*sizeof(double));
    for (int p=0;p<=nt;p++) t[p]=dt*p;
    cmp *psi=calloc(N*Nt,sizeof(cmp));
    cmp *Gamma=calloc(N,sizeof(cmp));

    FILE *fpsfile=fopen("fpsfile.csv","w");
    fprintf(fpsfile,"%d",fps);
    fclose(fpsfile);

    // punkt pierwszy
    FILE *p1misc=fopen("p1misc.csv","w");
    fprintf(p1misc,"%d,%d,%lf,%lf,",n,nt/co_ktora,L,sigmak);
    initGamma(Gamma, N);
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    FILE *p1 = fopen("p1.csv","w");
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p+=co_ktora){
            if (p!=0 || i!=0) fprintf(p1,",");
            fprintf(p1,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }
    fprintf(p1,"\n");
    sigmak*=1.5;
    fprintf(p1misc,"%lf,",sigmak);
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p+=co_ktora){
            if (p!=0 || i!=0) fprintf(p1,",");
            fprintf(p1,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }
    fprintf(p1,"\n");
    sigmak*=1.5;
    fprintf(p1misc,"%lf",sigmak);
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p+=co_ktora){
            if (p!=0 || i!=0) fprintf(p1,",");
            fprintf(p1,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }
    fprintf(p1misc,",%lf",dt);

    // punkt druig
    V0=0.3;
    initWch(Gamma,N,x,7.5,V0);
    sigmak=0.5;
    E=80.0;
    k0=sqrt(2.0*m*E)/hbar;
    initPsi(psi,n,nt,sigmak,x,x0,k0);
    solve(psi,Gamma,N,hbar,m,dx,dt,Nt);
    FILE *p2=fopen("p2.csv","w");
    for (int i=0;i<=n;i++){
        for (int p=0;p<=nt;p+=co_ktora){
            if (p!=0 || i!=0) fprintf(p2,",");
            fprintf(p2,"%lf",pow(cabs(psi[i*Nt+p]),2));
        }
    }
    FILE *p2gamma=fopen("p2gamma.csv","w");
    for (int i=0;i<=n;i++) fprintf(p2gamma,"%lf,%lf\n",creal(Gamma[i]),-cimag(Gamma[i]));




    // czystki
    free(psi);
    free(x);
    free(Gamma);
    free(t);
    fclose(p1);
    fclose(p1misc);
    fclose(p2);
    fclose(p2gamma);
    // return zero
    return 0;
}
