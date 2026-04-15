#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <string>

using namespace std;

__global__
void fill1D(double *x, int n, double d, double xa){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if (idx<=n) x[idx]=xa+idx*d;
}

__global__
void fillRho2D(double *rho, int nx, int ny, double *x, double *y, double R, double r){
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int idy = blockIdx.y*blockDim.y+threadIdx.y;
    if (idx <= nx && idy <=ny){
        if (hypot(x[idx],y[idy])<R){
            rho[idx*(ny+1)+idy]=r;
        }
        else{
            rho[idx*(ny+1)+idy]=0.0;
        }
    }
}
