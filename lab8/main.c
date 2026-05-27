#include "dek.h"

int main(){

    // parametry
    int N=100;
    double a=1.0;
    double eps=2.0;
    double t=-1.0;

    // alokacje / inicjalizacje
    cmp *F=malloc(N*N*sizeof(cmp));

    // defniniowanie F
    for (int l=0;l<N;l++){
        for (int m=0;m<N;m++){
            F[l*N+m]=cexp(2.0*M_PI*I*l*m/N)/sqrt(N);
        }
    }

    // sprawdzenie
    cmp *Fdagger=malloc(N*N*sizeof(cmp));
    cmp *jeden=malloc(N*N*sizeof(cmp));
    dagger(F,Fdagger,N);
    matmul(Fdagger,F,jeden,N);

    zapiszModulMacierzy(jeden,N,"jeden.csv");
    zapiszReMacierzy(F,N,"Fre.csv");
    zapiszImMacierzy(F,N,"Fim.csv");

    // hamiltonian
    cmp *H=calloc(N*N,sizeof(cmp));
    initH(H,N,t,eps);

    zapiszReMacierzy(H,N,"H.csv");

    // baza pedowa
    cmp *Hk=malloc(N*N*sizeof(cmp));
    matmul(Fdagger,H,jeden,N); // Fdagger * H -> temp
    matmul(jeden,F,Hk,N); // temp * F -> Hk

    zapiszReMacierzy(Hk,N,"Hk.csv");

    FILE *misc=fopen("misc.csv","w");
    fprintf(misc,"%d,%lf,%lf,%lf",N,a,eps,t);

    // czystki
    free(F);
    free(Fdagger);
    free(jeden);
    free(H);
    free(Hk);
    fclose(misc);

    // return zero
    return 0;
}
