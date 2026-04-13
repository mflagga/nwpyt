#include "header.h"

int main(){

    /* NIESKONCZONA STUDNIA POTENCJALU */

    // parametry
    double L = 5;
    int N = 500;
    double dx = L/N;
    double delta = 1e-10;
    double m=1.0;
    double hbar=1.0;
    double dtau=0.5*m*dx*dx/hbar;

    // alokacja/inicjalizacja
    double *x = malloc((N+1)*sizeof(double));
    double *V = calloc((N+1),sizeof(double));
    for (int i=0;i<=N;i++) x[i] = i*dx;
    double *psi = malloc((N+1)*sizeof(double));

    // metoda
    double E=met_czas_urojony(L,N,V,delta,psi,x,hbar,m,dx,dtau);

    // psi analityczne
    double *psi_analityczne = malloc((N+1)*sizeof(double));
    for (int i=0;i<=N;i++) psi_analityczne[i] = sqrt(2.0/L)*sin(M_PI*x[i]/L);

    // zapis do plikow
    FILE *Aa = fopen("Aa.csv","w");
    FILE *Amisc = fopen("Amisc.csv","w");
    for (int i=0;i<=N;i++) fprintf(Aa,"%lf,%lf,%lf\n",x[i],pow(fabs(psi[i]),2),pow(fabs(psi_analityczne[i]),2)); // zapis funkcji do pliku

    /* test dla dtau nie spełniającego warunku zbieżności */

    dtau=1.5*m*dx*dx/hbar;

    double Etest=met_czas_urojony(L,N,V,delta,psi,x,hbar,m,dx,dtau);

    fprintf(Amisc,"%lf,%d,%lf,%e,%lf",E,N,L,delta,Etest); // zapis parametrow do pliku
    FILE *Ab=fopen("Ab.csv","w");
    for (int i=0;i<=N;i++) fprintf(Ab,"%lf\n",pow(fabs(psi[i]),2));

    /* SKONCZONA STUDNIA POTENCJALU */

    // parametry
    dtau=0.5*m*dx*dx/hbar;
    double a=5.0;
    double L0=6.0;
    double V0=10.0;
    double xa, xb;
    int licznikmax=4;
    double plusrowne=2.0;
    L=L0;

    // alokacja
    double **ss_psi = malloc(licznikmax*sizeof(double*));
    for (int i=0;i<licznikmax;i++) ss_psi[i]=calloc(N+1,sizeof(double));
    double *ss_E=malloc(licznikmax*sizeof(double));

    FILE *BEL=fopen("BEL.csv","w");

    for (int licznik=0; licznik<licznikmax; licznik++){
        xa=L/2-a/2;
        xb=L/2+a/2;
        dx = L/N;
        dtau = 0.5*m*dx*dx/hbar;
        memset(V, 0.0, (N+1)*sizeof(double));
        for (int i=0;i<=N;i++) x[i] = dx*i;
        for (int i=0;i<=N;i++) if ((x[i]<=xa) || (x[i]>=xb)) V[i]=fabs(V0);
        ss_E[licznik]=met_czas_urojony(L,N,V,delta,ss_psi[licznik],x,hbar,m,dx,dtau);
        fprintf(BEL,"%lf,%lf\n",ss_E[licznik],L);
        L += plusrowne;
    }

    // zapis do plików
    FILE *Ba = fopen("Ba.csv","w");
    for (int i=0;i<=N;i++){
        for (int licznik=0;licznik<licznikmax;licznik++){
            if (licznik>0) fprintf(Ba,",");
            fprintf(Ba,"%lf",pow(fabs(ss_psi[licznik][i]),2));
        }
        fprintf(Ba,"\n");
    }
    FILE *Bmisc=fopen("Bmisc.csv","w");
    fprintf(Bmisc,"%lf,%lf,%d",a,fabs(V0),licznikmax);

    /* 2D nieskonczona studnia */

    // uwolnienie
    free(x);
    free(V);
    free(psi);

    // parametry
    N = 100;
    L = 5;
    dx = L/N;
    dtau = 0.5*(m*dx*dx)/(2*hbar);

    // alokacja
    x = malloc((N+1)*sizeof(double));
    double *y = malloc((N+1)*sizeof(double));
    V = calloc((N+1)*(N+1),sizeof(double));
    psi = calloc((N+1)*(N+1),sizeof(double));

    // inicjalizacja
    for (int i=0;i<=N;i++) x[i] = i*dx;
    for (int j=0;j<=N;j++) y[j] = j*dx;

    // metoda
    E=met_czas_urojony_2D(L,N,V,delta,psi,x,hbar,m,dx,dtau,y);

    // zapis
    FILE *Ca = fopen("Ca.csv","w");
    FILE *Cmisc = fopen("Cmisc.csv","w");
    for (int i=0;i<=N;i++){
        for (int j=0;j<=N;j++){
            fprintf(Ca,"%lf\n",pow(fabs(psi[i*(N+1)+j]),2));
        }
    }
    fprintf(Cmisc,"%d,%lf,%lf",N,L,E);

    /* 2D studnia z potencjalem gaussowskim */
    
    // pararmetry
    double sigma = L/8;
    V0=4.0;

    // inicjalizacja
    V_init_2D(V,V0,N,x,y,L,sigma);

    // metoda
    E=met_czas_urojony_2D(L,N,V,delta,psi,x,hbar,m,dx,dtau,y);

    // zapis
    FILE *Da=fopen("Da.csv","w");
    FILE *Dmisc=fopen("Dmisc.csv","w");
    for (int i=0;i<=N;i++){
        for (int j=0;j<=N;j++){
            fprintf(Da,"%lf\n",pow(fabs(psi[i*(N+1)+j]),2));
        }
    }
    fprintf(Dmisc,"%lf",E);

    // czystki
    free(x);
    free(y);
    free(V);
    free(psi);
    free(psi_analityczne);
    fclose(Aa);
    fclose(Amisc);
    fclose(Ab);
    for (int i=0;i<licznikmax;i++) free(ss_psi[i]);
    free(ss_psi);
    free(ss_E);
    fclose(Ba);
    fclose(BEL);
    fclose(Bmisc);
    fclose(Ca);
    fclose(Cmisc);
    fclose(Da);
    fclose(Dmisc);

    // return zero
    return 0;
}
