#include "header.cuh"

int main(){
    srand(time(NULL));
    // parametry
    int nx=128;
    int ny=128;
    double L=5.0;
    double eps=1.0;
    double Lx=L;
    double Ly=L;
    double dx=Lx/nx;
    double dy=Ly/ny;
    int Nx=nx+1;
    int Ny=ny+1;
    int tp1=256;
    int tp2=8;
    double R=1.0;
    double r=10.0;
    dim3 block2 = dim3(tp2,tp2);
    dim3 grid2 = dim3((nx+block2.x+1)/block2.x,(ny+block2.y+1)/block2.y);
    // alokacja
    double *x; cudaMalloc(&x,Nx*sizeof(double));
    double *y; cudaMalloc(&y,Ny*sizeof(double));
    double *sys_mask; cudaMalloc(&sys_mask,Nx*Ny*sizeof(double));
    // inicjalizacja
    fill1D<<<(nx+tp1)/tp1,tp1>>>(x,nx,dx,-Lx/2);
    fill1D<<<(ny+tp1)/tp1,tp1>>>(y,ny,dy,-Ly/2);
    cudaDeviceSynchronize();
    // czystki
    cudaFree(x);
    cudaFree(y);
    cudaFree(sys_mask);
    // return zero
    return 0;
}
