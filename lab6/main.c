#include "header.h"

int main(){

    // parametry jednostkowe
    double m=1.0;
    double hbar=sqrt(2.0*m*0.0381); // dla [nm] i [eV]

    // parametry układu
    double Lx=120.0;
    double Ly=50.0;
    int nx=360;
    int ny=150;
    double V0=0.02;
    double alpha=0.25;

    // parametry wtórne
    double dx=Lx/nx;
    double dy=Ly/ny;
    double ty=hbar*hbar/(2.0*m*dy*dy);

    // alokacja / inicjalizacja
    double *x=malloc((nx+1)*sizeof(double));
    for (int i=0;i<=nx;i++) x[i]=i*dx;
    double *y=malloc((ny+1)*sizeof(double));
    for (int j=0;j<=ny;j++) y[j]=j*dy;
    double *V=malloc((nx+1)*(ny+1)*sizeof(double));
    initV(V,nx,ny,x,y,Ly/2,Lx/2,V0,alpha);

    // zapisanie potencjału
    FILE *Vmap=fopen("Vmap.csv","w");
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            if (i!=0 || j!=0) fprintf(Vmap,",");
            fprintf(Vmap,"%lf",V[i*(ny+1)+j]);
        }
    }

    // zapisanie parametrów
    FILE *misc=fopen("misc.csv","w");
    fprintf(misc,"%lf,%lf,%d,%d",Lx,Ly,nx,ny);

    // czystki
    free(x);
    free(y);
    free(V);
    fclose(Vmap);
    fclose(misc);

    // return zero
    return 0;
}
