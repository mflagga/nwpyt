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
    double R=0.5;
    double r=10.0;
    double VB1=0.0;
    double VB2=0.0;
    double tol=1e-9;
    int co_ile=1000;
    double omega=1.0;
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
    printf("UKLAD PIERWSZY:\n");
    system1<<<grid2,block2>>>(sys_mask,rho,nx,ny,x,y,R,r);
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
    fprintf(S1misc,"%d,%d,%lf,%lf,%lf,%lf",nx,ny,VB1,VB2,Lx,Ly);
    FILE *S1rL=fopen("S1rL.csv","w");
    double Lvals[] = {5.0, 10.0, 15.0, 20.0};
    omega=1.9;
    for (double L: Lvals){
        Lx=L; Ly=L;
        dx=Lx/nx; dy=Ly/ny;
        fill1D<<<(nx+tp1)/tp1,tp1>>>(x,nx,dx,-Lx/2);
        fill1D<<<(ny+tp1)/tp1,tp1>>>(y,ny,dy,-Ly/2);
        cudaDeviceSynchronize();
        system1<<<grid2,block2>>>(sys_mask,rho,nx,ny,x,y,R,r);
        fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
        cudaDeviceSynchronize();
        petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
        cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
        for (int i=0;i<=nx;i++){
            for (int j=0;j<=ny;j++){
                if (i!=0 || j!=0) fprintf(S1rL,",");
                fprintf(S1rL,"%lf",vc[i*Ny+j]);
            }
        }
        fprintf(S1rL,"\n");
    }

    // uklad trzeci
    printf("UKLAD TRZECI:\n");
    L=5.0;
    Lx=L;
    Ly=L;
    // r=7; vb=-10;
    r=7.0;
    VB1=-10.0;
    dx=Lx/nx; dy=Ly/ny;
    fill1D<<<(nx+tp1)/tp1,tp1>>>(x,nx,dx,-Lx/2);
    fill1D<<<(ny+tp1)/tp1,tp1>>>(y,ny,dy,-Ly/2);
    cudaDeviceSynchronize();
    system3<<<grid2,block2>>>(sys_mask,rho,nx,ny,Lx,Ly,x,y,r);
    fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
    cudaDeviceSynchronize();
    petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
    cudaMemcpy(sys_maskc,sys_mask,Nx*Ny*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(rhoc,rho,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    FILE *S3a = fopen("S3a.csv","w");
    FILE *S3misc = fopen("S3misc.csv","w");
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            fprintf(S3a,"%lf,%d,%lf\n",vc[i*Ny+j],sys_maskc[i*Ny+j],rhoc[i*Ny+j]);
        }
    }
    fprintf(S3misc,"%d,%d,%lf,%lf,%lf,%lf",nx,ny,VB1,VB2,Lx,Ly);

    FILE *S3r=fopen("S3r.csv","w");

    // r=14; vb=-10;
    VB1=-10.0;
    r=14.0;
    system3<<<grid2,block2>>>(sys_mask,rho,nx,ny,Lx,Ly,x,y,r);
    fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
    cudaDeviceSynchronize();
    petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
    cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            if (i!=0 || j!=0) fprintf(S3r,",");
            fprintf(S3r,"%lf",vc[i*Ny+j]);
        }
    }
    fprintf(S3r,"\n");

    // r=7; vb=-20;
    VB1=-20.0;
    r=7.0;
    system3<<<grid2,block2>>>(sys_mask,rho,nx,ny,Lx,Ly,x,y,r);
    fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
    cudaDeviceSynchronize();
    petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
    cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            if (i!=0 || j!=0) fprintf(S3r,",");
            fprintf(S3r,"%lf",vc[i*Ny+j]);
        }
    }
    fprintf(S3r,"\n");

    // r=14; vb=-20
    VB1=-20.0;
    r=14.0;
    system3<<<grid2,block2>>>(sys_mask,rho,nx,ny,Lx,Ly,x,y,r);
    fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
    cudaDeviceSynchronize();
    petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
    cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            if (i!=0 || j!=0) fprintf(S3r,",");
            fprintf(S3r,"%lf",vc[i*Ny+j]);
        }
    }
    fprintf(S3r,"\n");

    // uklad drugi
    printf("UKLAD DRUGI:\n");
    Lx=10.0;
    Ly=5.0;
    ny=n;
    nx = (ny+1)*2-1;
    VB1=-10.0;
    VB2=-5.0;
    r=5.0;
    dx=Lx/nx;
    dy=Ly/ny;
    Nx=nx+1;
    Ny=ny+1;
    grid2 = dim3((nx+block2.x+1)/block2.x,(ny+block2.y+1)/block2.y);
    // alokacja/inicjalizacja na nowo
    cudaFree(x); cudaMalloc(&x,Nx*sizeof(double));
    cudaFree(y); cudaMalloc(&y,Ny*sizeof(double));
    cudaFree(sys_mask); cudaMalloc(&sys_mask,Nx*Ny*sizeof(int));
    cudaFree(rho); cudaMalloc(&rho,Nx*Ny*sizeof(double));
    cudaFree(v); cudaMalloc(&v,Nx*Ny*sizeof(double));
    fill1D<<<(nx+tp1)/tp1,tp1>>>(x,nx,dx,-Lx/2);
    fill1D<<<(ny+tp1)/tp1,tp1>>>(y,ny,dy,-Ly/2);
    cudaDeviceSynchronize();
    system2<<<grid2,block2>>>(sys_mask,rho,nx,ny,r);
    fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
    cudaDeviceSynchronize();
    delete [] vc; vc = new double[Nx*Ny];
    delete [] sys_maskc; sys_maskc = new int[Nx*Ny];
    delete [] rhoc; rhoc = new double[Nx*Ny];
    // rozwiazanie
    petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
    cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(sys_maskc,sys_mask,Nx*Ny*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(rhoc,rho,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
    
    FILE *S2a = fopen("S2a.csv","w");
    FILE *S2misc = fopen("S2misc.csv","w");
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
            fprintf(S2a,"%lf,%d,%lf\n",vc[i*Ny+j],sys_maskc[i*Ny+j],rhoc[i*Ny+j]);
        }
    }
    fprintf(S2misc,"%d,%d,%lf,%lf,%lf,%lf",nx,ny,VB1,VB2,Lx,Ly);

    // rozne wartosci VB1

    FILE *S2rV = fopen("S2rV.csv","w");
    for (int j=0;j<=ny;j++){
        if (j!=0) fprintf(S2rV,",");
        fprintf(S2rV,"%lf",vc[(nx/2)*Ny+j]);
    }
    fprintf(S2rV,"\n");
    double vbvals[] = {-15.0, -20.0, -25.0};
    for (double V: vbvals){
        VB1=V;
        fill2DwithC<<<grid2,block2>>>(v,nx,ny,r/2);
        cudaDeviceSynchronize();
        petla_relaksacyjna(block2,grid2,nx,ny,sys_mask,dx,dy,v,rho,VB1,VB2,eps,omega,co_ile,tol,Lx,Ly);
        cudaMemcpy(vc,v,Nx*Ny*sizeof(double),cudaMemcpyDeviceToHost);
        for (int j=0;j<=ny;j++){
            if (j!=0) fprintf(S2rV,",");
            fprintf(S2rV,"%lf",vc[(nx/2)*Ny+j]);
        }
        fprintf(S2rV,"\n");
    }

    // czystki
    cudaFree(x);
    cudaFree(y);
    cudaFree(sys_mask);
    cudaFree(rho);
    cudaFree(v);
    fclose(S1a);
    delete [] vc;
    delete [] sys_maskc;
    delete [] rhoc;
    fclose(S1misc);
    fclose(S1rL);
    fclose(S3a);
    fclose(S3misc);    
    fclose(S3r);
    fclose(S2a);
    fclose(S2misc);
    fclose(S2rV);

    // return zero
    return 0;
}
