#include "header.cuh"

int main(){
    // parametry
    int n=pow(2,5)-1;
    int nx=n;
    int ny=n;
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
    double R=0.3;
    double r=10.0;
    double VB1=-10.0;
    double VB2=10.0;
    double tol=1e-9;
    int co_ile=1000;
    double omega=1.9;
    dim3 block2 = dim3(tp2,tp2);
    dim3 grid2 = dim3((nx+block2.x+1)/block2.x,(ny+block2.y+1)/block2.y);
    // alokacja
    double *x; cudaMalloc(&x,Nx*sizeof(double));
    double *y; cudaMalloc(&y,Ny*sizeof(double));
    int *sys_mask; cudaMalloc(&sys_mask,Nx*Ny*sizeof(int));
    double *rho; cudaMalloc(&rho,Nx*Ny*sizeof(double));
    double *v; cudaMalloc(&v,Nx*Ny*sizeof(double));
    // inicjalizacja
    fill1D<<<(nx+tp1)/tp1,tp1>>>(x,nx,dx,-Lx/2);
    fill1D<<<(ny+tp1)/tp1,tp1>>>(y,ny,dy,-Ly/2);
    cudaDeviceSynchronize();

    // układ pierwszy
    system1<<<grid2,block2>>>(sys_mask,rho,nx,ny,x,y,R,r);
    // system3<<<grid2,block2>>>(sys_mask,rho,nx,ny,Lx,Ly,x,y,r);
    fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
    cudaDeviceSynchronize();
    petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
    double *vc = new double[Nx*Ny];
    cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    int *sys_maskc=new int[Nx*Ny];
    cudaMemcpy(sys_maskc,sys_mask,Nx*Ny*sizeof(int),cudaMemcpyDeviceToHost);
    double *rhoc=new double[Nx*Ny];
    cudaMemcpy(rhoc,rho,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    FILE *S1misc=fopen("S1misc.csv","w");
    FILE *S1a=fopen("S1a.csv","w");
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            fprintf(S1a,"%lf,%d,%lf\n",vc[i*Ny+j],sys_maskc[i*Ny+j],rhoc[i*Ny+j]);
        }
    }
    fprintf(S1misc,"%d,%d,%lf,%lf",nx,ny,VB1,VB2);


    // czystki
    cudaFree(x);
    cudaFree(y);
    cudaFree(sys_mask);
    cudaFree(rho);
    fclose(S1a);
    delete [] vc;
    delete [] sys_maskc;
    delete [] rhoc;
    fclose(S1misc);

    // return zero
    return 0;
}
