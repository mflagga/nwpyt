#include "header.h"

int main(){
    /* NIESKONCZONA STUDNIA POTENCJALU */
    // parametry
    double L = 5;
    int N = 1000;
    double dx = L/N;
    double delta = 1e-10;
    double m=1.0;
    double hbar=1.0;
    double dtau=0.5*m*dx*dx/hbar;
    // alokacja/inicjalizacja
    double *x = malloc((N+1)*sizeof(double));
    double *V = calloc((N+1),sizeof(double));
    for (int i=0;i<=N;i++){
        x[i] = i*dx;
    }
    double *psi = malloc((N+1)*sizeof(double));
    double E=met_czas_urojony(L,N,V,delta,psi,x,hbar,m,dx,dtau); // metoda
    double *psi_analityczne = malloc((N+1)*sizeof(double));
    for (int i=0;i<=N;i++) psi_analityczne[i] = sqrt(2.0/L)*sin(M_PI*x[i]/L);
    FILE *Aa = fopen("Aa.csv","w");
    FILE *Amisc = fopen("Amisc.csv","w");
    for (int i=0;i<=N;i++) fprintf(Aa,"%lf,%lf,%lf\n",x[i],pow(fabs(psi[i]),2),pow(fabs(psi_analityczne[i]),2)); // zapis funkcji do pliku
    // test dla dtau nie spełniającego warunku zbieżności
    dtau=1.5*m*dx*dx/hbar;
    double Etest=met_czas_urojony(L,N,V,delta,psi,x,hbar,m,dx,dtau);
    fprintf(Amisc,"%lf,%d,%lf,%e,%lf",E,N,L,delta,Etest); // zapis parametrow do pliku
    FILE *Ab=fopen("Ab.csv","w");
    for (int i=0;i<=N;i++) fprintf(Ab,"%lf\n",pow(fabs(psi[i]),2));
    /* SKONCZONA STUDNIA POTENCJALU */
    dtau=0.5*m*dx*dx/hbar;
    double a=5.0;
    double L0=6.0;
    double V0=10.0;
    double xa, xb;
    int licznikmax=5;
    double plusrowne=2.0;
    L=L0;
    double **ss_psi = malloc(licznikmax*sizeof(double*));
    for (int i=0;i<licznikmax;i++) ss_psi[i]=calloc(N+1,sizeof(double));
    double *ss_E=malloc(licznikmax*sizeof(double));
    FILE *BEL=fopen("BEL.csv","w");
    for (int licznik=0; licznik<licznikmax; licznik++){
        xa=L/2-a/2;
        xb=L/2+a/2;
        dx = L/N;
        for (int i=0;i<=N;i++) x[i] = dx*i;
        for (int i=0;i<=N;i++) if ((x[i]<=xa) || (x[i]>=xb)) V[i]=fabs(V0);
        ss_E[licznik]=met_czas_urojony(L,N,V,delta,ss_psi[licznik],x,hbar,m,dx,dtau);
        fprintf(BEL,"%lf,%lf\n",ss_E[licznik],L);
        L += plusrowne;
    }
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
    // czystki
    free(x);
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
    // return zero
    return 0;
}
