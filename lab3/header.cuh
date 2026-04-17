#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <string>

using namespace std;

#define typSrodek 0
#define typDirich1 1
#define typDirich2 2
#define typNeumann 3
#define typKarmann 4

__global__
void fill1D(double *x, int n, double d, double xa){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if (idx<=n) x[idx]=xa+idx*d;
}

__global__
void fill2DwithC(double *M, int nx, int ny, double C){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx<=nx && idy<=ny){
        M[idx*(ny+1)+idy]=C;
    }
}

__global__
void system1(int *sys_mask, double *rho, int nx, int ny, double *x, double *y, double R, double r){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    int p = idx*(ny+1)+idy; 
    if (idx<=nx && idy<=ny){
        if (idx==0 || idx==nx || idy==0 || idy==ny){
            sys_mask[p]=typDirich1;
        }
        else{
            sys_mask[p]=typSrodek;
        }
        if (sys_mask[p]==typSrodek){
            if (hypot(x[idx],y[idy])<R){
                rho[p]=r;
            }
            else{
                rho[p]=0.0;
            }
        }
    }
}

__global__
void system2(int *sys_mask, double *rho, int nx, int ny, double r){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx<=nx && idy<=ny){
        int p = idx*(ny+1)+idy; 
        sys_mask[p]=typSrodek;
        if (idx==0 || idx==nx) sys_mask[p]=typKarmann;
        if (idy==0) sys_mask[p]=typDirich1;
        if (idy==ny) sys_mask[p]=typDirich2;
        if (sys_mask[p]==typSrodek || sys_mask[p]==typKarmann) rho[p]=r;
    }
}

__global__
void system3(int *sys_mask, double *rho, int nx, int ny, double Lx, double Ly, double *x, double *y,double r){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx<=nx && idy<=ny){
        int p = idx*(ny+1)+idy;
        sys_mask[p]=typSrodek;
        if (idx==0 || idx==nx || idy==0 || idy==ny) sys_mask[p]=typNeumann;
        sys_mask[0]=typDirich1; sys_mask[ny]=typDirich1; sys_mask[nx*(ny+1)]=typDirich1; sys_mask[nx*(ny+1)+ny]=typDirich1;
        double a=Lx/10;
        double b=Lx/10;
        double c=Ly/4;
        double d=0.4*Ly;
        if ((x[idx]>(-b-3.0*a)) && (x[idx]<-a-b) && (y[idy]>Ly/2-d) && (y[idy]<=Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>-a) && (x[idx]<a) && (y[idy]>Ly/2-c) && (y[idy]<=Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>a+b) && (x[idx]<3.0*a+b) && (y[idy]>Ly/2-d) && (y[idy]<=Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>-b-3.0*a) && (x[idx]<-a-b) && (y[idy]>=-Ly/2) && (y[idy]<d-Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>-a) && (x[idx]<a) && (y[idy]>=-Ly/2) && (y[idy]<=c-Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>a+b) && (x[idx]<3.0*a+b) && (y[idy]>=-Ly/2) && (y[idy]<d-Ly/2)) sys_mask[p]=typDirich1;
        if (sys_mask[p]!=typDirich1){
            rho[p]=r;
        }
        else{
            rho[p]=0.0;
        }
    }
}

__global__
void SOrelaxHalf(int nx, int ny, int *sys_mask, double dx, double dy, double *v, double *rho, double V_brzegowe1, double V_brzegowe2, int zerolubjeden, double eps, double omega, double *blad, bool licz_norme){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx <= nx && idy <= ny){
        if ((idx+idy)%2==zerolubjeden){
            int p = idx*(ny+1)+idy;
            int typ = sys_mask[p];
            double inv_dx2=1.0/(dx*dx);
            double inv_dy2=1.0/(dy*dy);
            double diag=2.0*(inv_dx2+inv_dy2);
            double vGS;
            double vs=v[p];
            switch(typ){
                case typSrodek:
                    vGS=((v[p-ny-1]+v[p+ny+1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    v[p]=(1.0-omega)*v[p]+omega*vGS;
                    break;
                case typDirich1:
                    v[p] = V_brzegowe1;
                    break;
                case typDirich2:
                    v[p] = V_brzegowe2;
                    break;
                case typNeumann:
                    if (idx==0){
                        vGS=((2.0*v[p+ny+1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idx==nx){
                        vGS=((2.0*v[p-ny-1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idy==0){
                        vGS=((v[p+ny+1]+v[p-ny-1])*inv_dx2+(2.0*v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idy==ny){
                        vGS=((v[p+ny+1]+v[p-ny-1])*inv_dx2+(2.0*v[p-1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    v[p]=(1.0-omega)*v[p]+omega*vGS;
                    break;
                case typKarmann:
                    if (idx==0){
                        vGS=((v[nx*(ny+1)+idy]+v[p+ny+1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idx==nx){
                        vGS=((v[p-ny-1]+v[idy])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    v[p]=(1.0-omega)*v[p]+omega*vGS;
                    break;
            }
            if (licz_norme){
                double diff=v[p]-vs;
                atomicAdd(blad,diff*diff);
            }
        }
    }
}

void SO_relax_full_once(dim3 block2, dim3 grid2, int nx, int ny, int *sys_mask, double dx, double dy, double *v, double *rho, double V_brzegowe1, double V_brzegowe2, double eps, double omega, double *blad, bool licz_norme){
    SOrelaxHalf<<<grid2,block2>>>(nx,ny,sys_mask,dx,dy,v,rho,V_brzegowe1,V_brzegowe2,1,eps,omega,blad,licz_norme);
    cudaDeviceSynchronize();
    SOrelaxHalf<<<grid2,block2>>>(nx,ny,sys_mask,dx,dy,v,rho,V_brzegowe1,V_brzegowe2,0,eps,omega,blad,licz_norme);
    cudaDeviceSynchronize();
}

void petla_relaksacyjna(dim3 block2, dim3 grid2, int nx, int ny, int *sys_mask, double dx, double dy, double *v, double *rho, double V_brzegowe1, double V_brzegowe2, double eps, double omega, int co_ile, double tol, double Lx, double Ly){
    double *blad; cudaMalloc(&blad, sizeof(double));
    double bladhost=0.0;
    int itmax=100000;
    bool licz_norme;
    auto t0 = chrono::high_resolution_clock::now();
    printf("(Lx,Ly)=(%.2lf,%.2lf); (nx,ny)=(%d,%d); omega=%.2lf; tol=%.0e; ",Lx,Ly,nx,ny,omega,tol);
    for (int i=1;i<=itmax;i++){
        licz_norme = i%co_ile==0;
        if(licz_norme){
            bladhost=0.0;
            cudaMemcpy(blad,&bladhost,sizeof(double),cudaMemcpyHostToDevice);
        }
        SO_relax_full_once(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,V_brzegowe1,V_brzegowe2,eps,omega,blad,licz_norme);
        if (licz_norme){
            cudaMemcpy(&bladhost,blad,sizeof(double),cudaMemcpyDeviceToHost);
            bladhost=sqrt(bladhost);
            if (bladhost<tol){
                cout<<"Zbieglo sie w "<<i<<" iteracjach ";
                break;
            }
        }
        if (i==itmax) cout<<"Osiegnieto itmax="<<itmax<<" ";
    }
    cudaFree(blad);
    auto t1 = chrono::high_resolution_clock::now();
    cout<<"w czasie "<<chrono::duration_cast<chrono::milliseconds>(t1-t0).count()/1000.0<<" s\n";
}
