#include "header.h"

int main(){

    // parametry jednostkowe
    double m=1.0;
    double hbar=sqrt(2.0*m*0.0381); // dla [nm] i [eV]
    double e=1.0;

    // parametry układu
    double Lx=60.0;
    double Ly=25.0;
    int nx=720;
    int ny=300;
    double V0=0.02;
    double alpha=0.5;
    double LE=10*V0;
    int nE=400;

    // parametry wtórne
    double dx=Lx/nx;
    double dy=Ly/ny;
    double ty=hbar*hbar/(2.0*m*dy*dy);
    double tx=hbar*hbar/(2.0*m*dx*dx);
    double dE=LE/nE;

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

    // obliczanie nowego potencjalu 
    double *Vn=malloc((nx+1)*(ny-1)*sizeof(double));
    double *evals=malloc((ny-1)*sizeof(double));
    for (int i=0;i<=nx;i++){
        diagYslice(i,V,ny,evals,ty);
        for (int j=0;j<ny-1;j++) Vn[i*(ny-1)+j]=evals[j];
    }

    // zapisanie potencjału nowego
    FILE *Vnmap=fopen("Vnmap.csv","w");
    for (int i=0;i<=nx;i++){
        for (int j=0;j<ny-1;j++){
            if (i!=0 || j!=0) fprintf(Vnmap,",");
            fprintf(Vnmap,"%lf",Vn[i*(ny-1)+j]);
        }
    }

    // rozwiazanie metoda qtbm
    double *Vr=malloc((nx+1)*sizeof(double));
    cmp *psi=malloc((nx+1)*sizeof(cmp));
    double *T=malloc(nE*(ny-1)*sizeof(double));
    double E;
    for (int j=0;j<ny-1;j++){
        for (int i=0;i<=nx;i++) Vr[i]=Vn[i*(ny-1)+j];
        for (int l=0;l<nE;l++){
            E=dE*(l+1);
            T[j*nE+l]=TodE(E,psi,Vr,nx+1,tx,dx,m,hbar);
        }
    }
    double *G=calloc(nE,sizeof(double));
    for (int l=0;l<nE;l++){
        for (int j=0;j<ny-1;j++){
            G[l] += T[j*nE+l];
        }
        G[l] *= e*e/(M_PI*hbar);
    }

    // zapis konduktacji od energii
    FILE *Gfile=fopen("Gfile.csv","w");
    for (int l=0;l<nE;l++){
        E=dE*(l+1);
        fprintf(Gfile,"%lf,%lf\n",E,G[l]);
    }

    // czystki
    free(x);
    free(y);
    free(V);
    fclose(Vmap);
    fclose(misc);
    free(Vn);
    free(evals);
    fclose(Vnmap);
    free(Vr);
    free(psi);
    free(T);
    free(G);
    fclose(Gfile);

    // return zero
    return 0;
}
