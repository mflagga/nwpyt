#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <string>

using namespace std;

#define typSrodek 0
#define typDirich1 1
#define typNeumann 2
#define typKarmann 3
#define typDirich2 4

__global__
void fill1D(double *x, int n, double d, double xa){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if (idx<=n) x[idx]=xa+idx*d;
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
void system3(int *sys_mask, double *rho, int nx, int ny, double Lx, double Ly, double *x, double *y){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx<=nx && idy<=ny){
        int p = idx*(ny+1)+idy;
        sys_mask[p]=typSrodek;
        if (idx==0 || idx==nx || idy==0 || idy==ny) sys_mask[p]=typNeumann;
        sys_mask[0]=typDirich1; sys_mask[ny]=typDirich1; sys_mask[nx*(ny+1)]=typDirich1; sys_mask[nx*(ny+1)+ny]=typDirich1;
        double a=Lx/10;
        double b=Lx/10;
        double c=Ly/5;
        double d=0.4*Ly;
        if ((x[idx]>-b-1.5*a) && (x[idx]<-a-b) && (y[idy]>Ly/2-d) && (y[idy]<=ny)) sys_mask[p]=typDirich1;
        if ((x[idx]>-a) && (x[idx]<a) && (y[idy]>Ly/2-c) && (y[idy]<=ny)) sys_mask[p]=typDirich1;
        if ((x[idx]>a+b) && (x[idx]<1.5*a+b) && (y[idy]>Ly/2-d) && (y[idy]<=ny)) sys_mask[p]=typDirich1;
        if ((x[idx]>-b-1.5*a) && (x[idx]<-a-b) && (y[idy]>=0) && (y[idy]<d-Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>-a) && (x[idx]<a) && (y[idy]>=0) && (y[idy]<=c-Ly/2)) sys_mask[p]=typDirich1;
        if ((x[idx]>a+b) && (x[idx]<1.5*a+b) && (y[idy]>=0) && (y[idy]<d-Ly/2)) sys_mask[p]=typDirich1;
    }
}

__global__
void GSrelaxHalf(int nx, int ny, int *sys_mask, double dx, double dy, double *v, double *rho, double V_brzegowe1, double V_brzegowe2, int zerolubjeden, double eps){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx <= nx && idy <= ny){
        if ((idx+idy)%2==zerolubjeden){
            int p = idx*(ny+1)+idy;
            int typ = sys_mask[p];
            double inv_dx2=1.0/(dx*dx);
            double inv_dy2=1.0/(dy*dy);
            double diag=2.0*(inv_dx2+inv_dy2);
            switch(typ){
                case typSrodek:
                    v[p]=((v[p-ny-1]+v[p+ny+1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    break;
                case typDirich1:
                    v[p] = V_brzegowe1;
                    break;
                case typDirich2:
                    v[p] = V_brzegowe2;
                    break;
                case typNeumann:
                    if (idx==0){
                        v[p]=((2.0*v[p+ny+1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idx==nx){
                        v[p]=((2.0*v[p-ny-1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idy==0){
                        v[p]=((v[p+ny+1]+v[p-ny-1])*inv_dx2+(2.0*v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idy==ny){
                        v[p]=((v[p+ny+1]+v[p-ny-1])*inv_dx2+(2.0*v[p-1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    break;
                case typKarmann:
                    if (idx==0){
                        v[p]=((v[nx*(ny+1)+idy]+v[p+ny+1])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
                    else if(idx==nx){
                        v[p]=((v[p-ny-1]+v[idy])*inv_dx2+(v[p-1]+v[p+1])*inv_dy2+rho[p]/eps)/diag;
                    }
            }
        }
    }
}

void GS_relax_full_once(dim3 block2, dim3 grid2, int nx, int ny, int *sys_mask, double dx, double dy, double *v, double *rho, double V_brzegowe1, double V_brzegowe2, double eps){
    GSrelaxHalf<<<grid2,block2>>>(nx,ny,sys_mask,dx,dy,v,rho,V_brzegowe1,V_brzegowe2,1,eps);
    cudaDeviceSynchronize();
    GSrelaxHalf<<<grid2,block2>>>(nx,ny,sys_mask,dx,dy,v,rho,V_brzegowe1,V_brzegowe2,0,eps);
    cudaDeviceSynchronize();
}
