#include "dek.h"

void matmul(cmp *A, cmp *B, cmp *C, int N){
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            C[i*N+j]=0.0+0.0*I;
            for (int k=0;k<N;k++){
                C[i*N+j] += A[i*N+k]*B[k*N+j];
            }
        }
    }
}

void dagger(cmp *A, cmp *Adagger, int N){
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            Adagger[i*N+j]=conj(A[j*N+i]);
        }
    }
}

void initH(cmp *H, int N, double t, double eps){ // zakladajac ze jest calloc
    for (int i=0;i<N;i++){
        H[i*N+i]=eps;
        if (i!=N-1){
            H[(i+1)*N+i]=t;
            H[i*N+i+1]=t;
        }
    }
    H[N-1]=t;
    H[(N-1)*N]=t;
}

void zapiszModulMacierzy(cmp *M, int N, char *nazwa){
    FILE *file=fopen(nazwa,"w");
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            if (i!=0 || j!=0) fprintf(file,",");
            fprintf(file,"%e",cabs(M[i*N+j]));
        }
    }
    fclose(file);
}

void zapiszReMacierzy(cmp *M, int N, char *nazwa){
    FILE *file=fopen(nazwa,"w");
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            if (i!=0 || j!=0) fprintf(file,",");
            fprintf(file,"%e",creal(M[i*N+j]));
        }
    }
    fclose(file);
}

void zapiszImMacierzy(cmp *M, int N, char *nazwa){
    FILE *file=fopen(nazwa,"w");
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            if (i!=0 || j!=0) fprintf(file,",");
            fprintf(file,"%e",cimag(M[i*N+j]));
        }
    }
    fclose(file);
}
