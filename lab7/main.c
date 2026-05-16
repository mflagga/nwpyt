#include "header.h"

int main(){
    
    // parametry jednostkowe
    double m=1.0;
    double hbar=1.0;

    // parametry układu
    double a=1.0;
    double V0=50.0;
    int N=500;
    int Nk=30;

    // parametry wtórne
    double w=a/2;
    double dx=a/N;
    double t=hbar*hbar/(2.0*m*dx*dx);
    double dk=2.0*M_PI/(a*Nk);

    // alokacje inicjalizacje
    double *x=malloc(N*sizeof(double));
    for (int i=0;i<N;i++) x[i]=-a/2+i*dx;
    double *V=malloc(N*sizeof(double));
    initV(V,N,V0,w,x);

    // rozwiazanie RS
    double *macierzPasmowa=malloc(N*(Nk+1)*sizeof(double));
    double *ev=malloc(N*sizeof(double));
    double k;
    for (int p=0;p<=Nk;p++){
        k=-M_PI/a + p*dk;
        diagHdlaK(k,N,dx,t,V,ev);
        for (int j=0;j<N;j++){
            macierzPasmowa[p*N+j]=ev[j];
            // tu trzeba dodac szukanie mina i maaxa dla danego poziomu E i wypelnieneie wektora W + zmienna max_poprzedni i wypelnienie wektora deltaE
        }
    }

    // szukanie szerokosci i przerw energetycznych w pasmach
    double min,max,maxp;
    double *W=malloc(N*sizeof(double));
    double *deltaE=malloc(N*sizeof(double));
    for (int j=0;j<N;j++){
        min=macierzPasmowa[j];
        max=min;
        for (int p=0;p<=Nk;p++){
            if (macierzPasmowa[p*N+j]>max) max=macierzPasmowa[p*N+j];
            if (macierzPasmowa[p*N+j]<min) min=macierzPasmowa[p*N+j];
        }
        if (j!=0) deltaE[j-1]=min-maxp;
        W[j]=max-min;
        maxp=max;
    }
    deltaE[N-1]=deltaE[N-2];

    // zapisanie 
    FILE *mp=fopen("mp.csv","w");
    for (int p=0;p<=Nk;p++){
        for (int j=0;j<N;j++){
            if (p!=0 || j!=0) fprintf(mp,",");
            fprintf(mp,"%lf",macierzPasmowa[p*N+j]);
        }
    }

    FILE *rsV=fopen("rsV.csv","w");
    for (int i=0;i<N;i++) fprintf(rsV,"%lf,%lf\n",x[i],V[i]);

    FILE *misc=fopen("misc.csv","w");
    fprintf(misc,"%d,%d,%lf",N,Nk,a);

    FILE *przerwy=fopen("przerwy.csv","w");
    for (int j=0;j<N;j++) fprintf(przerwy,"%lf,%lf\n",W[j],deltaE[j]);

    // metoda TBA
    double xa=-2*a;
    double xb=3*a;
    dx=(xb-xa)/(N-1);
    t=hbar*hbar/(2.0*m*dx*dx);
    for (int i=0;i<N;i++) x[i] = xa+i*dx;
    initV(V,N,V0,w,x);
    double E;
    double *phi=malloc(N*sizeof(double));
    diagStud(N,t,V,&E,phi,dx);

    int ia=(int)round(a/dx);
    double *Vtot=calloc(N,sizeof(double));
    for (int j=0;j<N;j++){
        for (int p=-1;p<=1;p++){
            if (j-p*ia>=0 && j-p*ia<N) Vtot[j] += V[j-p*ia];
        }
    }
    
    FILE *tbaV=fopen("tbaV.csv","w");
    for (int i=0;i<N;i++) fprintf(tbaV,"%lf,%lf,%lf,%lf\n",x[i],V[i],pow(fabs(phi[i]),2),Vtot[i]);

    double t01=0.0;
    for (int i=0;i<N;i++){
        if (i-ia>=0) t01 += phi[i]*(E+Vtot[i]-V[i-ia])*phi[i-ia]*dx;
    }

    double s01=0.0;
    for (int i=0;i<N;i++){
        if (i-ia>=0) s01 += phi[i]*phi[i-ia]*dx;
    }

    double eps;
    FILE *tbaE=fopen("tbaE.csv","w");
    for (int p=0;p<=Nk;p++){
        k=-M_PI/a + p*dk;
        eps=(E+2.0*t01*cos(k*a))/(1.0+2*s01*cos(k*a));
        if (p!=0) fprintf(tbaE,",");
        fprintf(tbaE,"%lf",eps);
    }
    
    // czystki
    free(x);
    free(V);
    free(macierzPasmowa);
    free(ev);
    fclose(mp);
    fclose(misc);
    fclose(rsV);
    fclose(tbaV);
    free(phi);
    free(Vtot);
    fclose(tbaE);
    free(W);
    free(deltaE);
    fclose(przerwy);

    // return zero
    return 0;
}
