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
    double *V = malloc((N+1)*sizeof(double));
    for (int i=0;i<=N;i++){
        x[i] = i*dx;
        V[i] = 0.0;
    }
    double *psi = malloc((N+1)*sizeof(double));
    double E=met_czas_urojony(L,N,V,delta,psi,x,hbar,m,dx,dtau); // metoda
    double *psi_analityczne = malloc((N+1)*sizeof(double));
    for (int i=0;i<=N;i++) psi_analityczne[i] = sqrt(2.0/L)*sin(M_PI*x[i]/L);
    FILE *Aa = fopen("Aa.csv","w");
    FILE *Amisc = fopen("Amisc.csv","w");
    for (int i=0;i<=N;i++) fprintf(Aa,"%lf,%lf,%lf\n",x[i],pow(fabs(psi[i]),2),pow(fabs(psi_analityczne[i]),2)); // zapis funkcji do pliku
    fprintf(Amisc,"%lf,%d,%lf,%e",E,N,L,delta); // zapis parametrow do pliku
    /* SKONCZONA STUDNIA POTENCJALU */
    
    // czystki
    free(x);
    free(V);
    free(psi);
    free(psi_analityczne);
    fclose(Aa);
    fclose(Amisc);
    // return zero
    return 0;
}
