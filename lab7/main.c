#include "header.h"

int main(){
    
    // parametry jednostkowe
    double m=1.0;
    double hbar=1.0;

    // parametry układu
    double a=1.0;
    double V0=100.0;
    int N=200;
    int Nk=20;

    // parametry wtórne
    double w=a/2;
    // double xa=-2*a;
    // double xb=3*a;
    double dx=a/N;
    double t=hbar*hbar/(2.0*m*dx*dx);
    double dk=2.0*M_PI/(a*Nk);

    // alokacje inicjalizacje
    double *x=malloc(N*sizeof(double));
    for (int i=0;i<=N;i++) x[i]=-a/2+i*dx;
    double *V=calloc(N,sizeof(double));
    initV(V,N,V0,w,x);

    // rozwiazanie 1
    double *macierzPasmowa=malloc(N*(Nk+1)*sizeof(double));
    double *ev=malloc(N*sizeof(double));
    double k;
    for (int p=0;p<=Nk;p++){
        k=-M_PI/a + p*dk;
        diagHdlaK(k,N,dx,t,V,ev);
        for (int j=0;j<N;j++){
            macierzPasmowa[p*N+j]=ev[j];
        }
    }

    // zapisanie 
    FILE *mp=fopen("mp.csv","w");
    for (int p=0;p<=Nk;p++){
        for (int j=0;j<N;j++){
            if (p!=0 || j!=0) fprintf(mp,",");
            fprintf(mp,"%lf",macierzPasmowa[p*N+j]);
        }
    }

    FILE *rsV=fopen("rsV.csv","w");


    FILE *misc=fopen("misc.csv","w");
    fprintf(misc,"%d,%d,%lf",N,Nk,a);
    for (int i=0;i<N;i++) fprintf(rsV,"%lf,%lf\n",x[i],V[i]);

    // czystki
    free(x);
    free(V);
    free(macierzPasmowa);
    free(ev);
    fclose(mp);
    fclose(misc);
    fclose(rsV);

    // return zero
    return 0;
}
