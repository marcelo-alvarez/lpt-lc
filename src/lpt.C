#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "allvars.h"

void Displace_1LPT(float *d, float *vx, float *vy, float *vz)
{

  int i,j,k,r,index;  

  // Local copies of parameters
  int N=Parameters.N;
  float BoxSize=Parameters.BoxSize;

  // Cell sizes and volumes in k space
  float dk=2.*M_PI/BoxSize;
  float d3k=pow(dk,3);
  float kmax=dk*N;

  // FFT of delta to k-space
  rfftwnd_mpi(plan, 1, d, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(d, N, size_fftw, 1.5);

  // Convert to displacements
  for(i=0;i<Nlocal;i++){
    float kx=(i+Nlocal*myid)*dk;
    if(i+Nlocal*myid>=N/2) kx=kx-kmax;
    float kx2=kx*kx;
    
    for(j=0;j<N;j++){
      float ky=j*dk;
      if(j>=N/2) ky=ky-kmax;
      float ky2=ky*ky;
      
      for(k=0;k<=N/2;k++){
	float kz=k*dk;
	if(k>=N/2) kz=kz-kmax;
	float kz2=kz*kz;
	
	int index=i*N*(N+2)+j*(N+2)+2*k;
	
	float rk2=kx2+ky2+kz2;
	float rk=sqrt(rk2);

	if(i+Nlocal*myid==N/2){
	  vx[index]   = 0;
	  vx[index+1] = 0;
	} else{
	  vx[index]   = -kx/rk2*d[index+1];
	  vx[index+1] =  kx/rk2*d[index];
	}

	if(j==N/2){
	  vy[index]   = 0;
	  vy[index+1] = 0;
	} else{
	  vy[index]   = -ky/rk2*d[index+1];
	  vy[index+1] =  ky/rk2*d[index];
	}

	if(k==N/2){
	  vz[index]   = 0;
	  vz[index+1] = 0;
	} else{
	  vz[index]   = -kz/rk2*d[index+1];
	  vz[index+1] =  kz/rk2*d[index];
	}

	if(rk2==0){
	  vx[index]=0; vx[index+1]=0;
	  vy[index]=0; vy[index+1]=0;
	  vz[index]=0; vz[index+1]=0;
	}

      }
    }
  }

  // Convert density back to real space
  rfftwnd_mpi(iplan, 1, d, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(d, N, size_fftw, 1.5);

  // Convert displacements back to real space
  rfftwnd_mpi(iplan, 1, vx, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vy, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vz, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(vx, N, size_fftw, 1.5);
  NormalizeArray(vy, N, size_fftw, 1.5);
  NormalizeArray(vz, N, size_fftw, 1.5);

  if(myid==0) printf("\n displacements calculated...");
  
  return;
  
}

void Displace_2LPT(float  *d1, float  *d2, 
			 float *vxx, float *vyy, float *vzz, 
			 float *vxy, float *vxz, float *vyz)
{

  int i,j,k,r,index;

  // Local copies of parameters
  int N=Parameters.N;
  float BoxSize=Parameters.BoxSize;

  // Cell sizes and volumes in k space
  float dk=2.*M_PI/BoxSize;
  float d3k=pow(dk,3);
  float kmax=dk*N;

  // FFT of delta to k-space
  rfftwnd_mpi(plan, 1, d1, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(d1, N, size_fftw, 1.5);

  // Convert to displacement
  for(i=0;i<Nlocal;i++){
    float kx=(i+Nlocal*myid)*dk;
    if(i+Nlocal*myid>=N/2) kx=kx-kmax;
    float kx2=kx*kx;
    
    for(j=0;j<N;j++){
      float ky=j*dk;
      if(j>=N/2) ky=ky-kmax;
      float ky2=ky*ky;
      
      for(k=0;k<=N/2;k++){
	float kz=k*dk;
	if(k>=N/2) kz=kz-kmax;
	float kz2=kz*kz;
	
	int index=i*N*(N+2)+j*(N+2)+2*k;
	
	float rk2=kx2+ky2+kz2;
	float rk=sqrt(rk2);

	if(i+Nlocal*myid==N/2){
	  vxx[index]   = 0;
	  vxx[index+1] = 0;
	  vxy[index]   = 0;
	  vxy[index+1] = 0;
	  vxz[index]   = 0;
	  vxz[index+1] = 0;
	} else{
	  vxx[index]   = kx*kx/rk2*d1[index];
	  vxx[index+1] = kx*kx/rk2*d1[index+1];
	  vxy[index]   = kx*ky/rk2*d1[index];
	  vxy[index+1] = kx*ky/rk2*d1[index+1];
	  vxz[index]   = kx*kz/rk2*d1[index];
	  vxz[index+1] = kx*kz/rk2*d1[index+1];
	}

	if(j==N/2){
	  vyy[index]   = 0;
	  vyy[index+1] = 0;
	  vxy[index]   = 0;
	  vxy[index+1] = 0;
	  vyz[index]   = 0;
	  vyz[index+1] = 0;
	} else{
	  vyy[index]   = ky*ky/rk2*d1[index];
	  vyy[index+1] = ky*ky/rk2*d1[index+1];
	  vyz[index]   = ky*kz/rk2*d1[index];
	  vyz[index+1] = ky*kz/rk2*d1[index+1];
	}

	if(k==N/2){
	  vzz[index]   = 0;
	  vzz[index+1] = 0;
	} else{
	  vzz[index]   = kz*kz/rk2*d1[index];
	  vzz[index+1] = kz*kz/rk2*d1[index+1];
	}

	if(rk2==0){
	  vxx[index]=0; vxx[index+1]=0;
	  vxy[index]=0; vxy[index+1]=0;
	  vxz[index]=0; vxz[index+1]=0;
	  vyy[index]=0; vyy[index+1]=0;
	  vyz[index]=0; vyz[index+1]=0;
	  vzz[index]=0; vzz[index+1]=0;
	}

      }
    }
  }

  // Convert displacements back to real space
  rfftwnd_mpi(iplan, 1, vxx, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vyy, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vzz, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vxy, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vxz, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, vyz, NULL, FFTW_NORMAL_ORDER);

  NormalizeArray(vxx, N, size_fftw, 1.5);
  NormalizeArray(vxy, N, size_fftw, 1.5);
  NormalizeArray(vxz, N, size_fftw, 1.5);
  NormalizeArray(vyy, N, size_fftw, 1.5);
  NormalizeArray(vyz, N, size_fftw, 1.5);
  NormalizeArray(vzz, N, size_fftw, 1.5);

  // Now apply 2LPT to obtain delta_2 in real space
  for(int ic=0;ic<Nlocal;ic++){
    for(int jc=0;jc<N;jc++){
      for(int kc=0;kc<N;kc++){

	index = ic*N*(N+2) + jc*(N+2) + kc;
	
	d2[index] = 0.0 
	  + vxx[index]*vyy[index]  
	  + vxx[index]*vzz[index] 
	  + vyy[index]*vzz[index] 
	  - vxy[index]*vxy[index]
	  - vyz[index]*vyz[index]
	  - vxz[index]*vxz[index];
	d2[index] = d2[index] * -1.0;
      }
    }
  }

  // reset xdv to zero
  for(int ic=0;ic<size_fftw;ic++){vxx[ic]=0;vyy[ic]=0;vzz[ic]=0;}

  if(myid==0) printf("\n 2LPT calculated...");
  
  return;
  
}

int ijk2index(int i, int j, int k){
   
   int index;
   long N=Parameters.N;

   index = i*N*(N+2)+j*(N+2)+k;

   return index;

}
