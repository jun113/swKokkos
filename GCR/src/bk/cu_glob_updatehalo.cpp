#include "hip/hip_runtime.h"

#include "gcr.h"
#include <hip/hip_runtime.h>
#include <hipblas.h>
#include <iostream>
#include <fstream>
#include <mpi.h>

//extern "C" void glob_updatehalo_x0_(double *x0, int *ims, int *ime, int *jms, int *jme, int *kts, int *kte);
//extern "C" void module_halo_mp_glob_updatehalo_real_3d_(float* h_p, int *,int *); 
extern "C" void glob_updatehalo_p_(float *p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);


#include "hip/hip_runtime.h"

#define index_buff1(i,k,j) (k-kts+1) + (j-jts)*(kte-kts+3) + (i) * (kte-kts+3)*(jte-jts+1)
#define index_buff2(i,k,j) i-its+1 +(k-kts+1)*(ite-its+3) + (j)*(kte-kts+3)*(ite-its+3) + 2*(kte-kts+3)*(jte-jts+1)

__global__ void gcr_memd2h(float* __restrict__ d_p, float* __restrict__ d_buff,int its, int ite, int jts,int jte,int kts, int kte){
        int NX = ite - its + 3;
	int NY = kte - kts + 3;
        int NZ = jte - jts + 3;
        //int jend=min(jde-1,jte);
        //jend=jte;
        int i= blockDim.x * blockIdx.x + threadIdx.x + its-1;
        int k= blockDim.y * blockIdx.y + threadIdx.y +kts-1;
        int j= blockDim.x * blockIdx.x + threadIdx.x + jts -((NX-1)/blockDim.x+1)*blockDim.x;

        if(j>=jts&&j<=jte){
            if(k>=kts-1&&k<=kte+1){
            d_buff[index_buff1(0,k,j)]=d_p[index4b(its,k,j,0)];
            d_buff[index_buff1(1,k,j)]=d_p[index4b(ite,k,j,0)];
        }}
        if(i>=its-1&&i<=ite+1){
            if(k>=kts-1&&k<=kte+1){
            d_buff[index_buff2(i,k,0)]=d_p[index4b(i,k,jts,0)];
            d_buff[index_buff2(i,k,1)]=d_p[index4b(i,k,jte,0)];
        }}

}
__global__ void gcr_memh2d(float* __restrict__ d_p, float* __restrict__ d_buff,int its, int ite, int jts,int jte,int kts, int kte){
	int NX = ite - its + 3;
	int NY = kte - kts + 3;
        int NZ = jte - jts + 3;
        //int jend=min(jde-1,jte);
        //jend=jte;
        int i= blockDim.x * blockIdx.x + threadIdx.x + its-1;
        int k= blockDim.y * blockIdx.y + threadIdx.y +kts-1;
        int j= blockDim.x * blockIdx.x + threadIdx.x + jts - ((NX-1)/blockDim.x+1)*blockDim.x;
        if(j>=jts&&j<=jte){
            if(k>=kts-1&&k<=kte+1){
                d_p[index4b(its-1,k,j,0)]=d_buff[index_buff1(0,k,j)];
                d_p[index4b(ite+1,k,j,0)]=d_buff[index_buff1(1,k,j)];
        }}
        if(i>=its-1&&i<=ite+1){
            if(k>=kts-1&&k<=kte+1){
            d_p[index4b(i,k,jts-1,0)]=d_buff[index_buff2(i,k,0)];
            d_p[index4b(i,k,jte+1,0)]=d_buff[index_buff2(i,k,1)];
        }}
}
//#include "cu_glob_updatehalo.hpp"

       





extern void glob_updatehalo_p_(float *h_p, int *m_tmp, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);

extern "C" void cu_glob_updatehalo_p(float* h_p, float* d_p, float* buff, float* d_buff,
    int its, int ite, int kts, int kte, int jts, int jte, 
        int m_tmp, int iter_max, int* mytid) {
  
  /*
  Size of buff: 2 * NY * (NX + NZ) 

  1. d_p -> d_buff
  2. d_buff -> h_buff
  3. h_buff -> h_p
  4. glob_updatehalo(h_p)
  5. h_p -> h_buff
  6. h_buff -> d_buff
  7. d_buff -> d_p
  */
  
  int NX = ite - its + 3;
  int NY = kte - kts + 3;
  int NZ = jte - jts + 3;
  int NG1=NX*NY*NG1;
  int zero=0;
  //int size_buff = 4 * (NX + NZ) * NY * sizeof(float);

  //int stride_1 = (NX + NZ) * NY;
  //int stride_2 = NX * NY;


    dim3 threadsPerBlock1(8, 4, 1);
    dim3 numBlocks1(ceil(((float)NX)/8)+ceil(((float)NZ)/8),ceil(((float) NY)/4),1);

        hipLaunchKernelGGL(gcr_memd2h, numBlocks1, threadsPerBlock1, 0, 0,  d_p+m_tmp*NG1,d_buff, its, ite, jts, jte, kts, kte);
    	CHECK(hipDeviceSynchronize()); 
        CHECK(hipMemcpy(buff,d_buff,2*(NX*NY+NY*NZ)*sizeof(float),hipMemcpyDeviceToHost));
        //boundary_init(h_p[index4b(i,k,j,0)]=0.0;)
        for(int j=jts;j<=jte;j++){
        for(int k=kts-1;k<=kte+1;k++){
                h_p[index4b(its,k,j,zero)]=buff[index_buff1(0,k,j)];
                h_p[index4b(ite,k,j,zero)]=buff[index_buff1(1,k,j)];
        }}

        for(int i=its-1;i<=ite+1;i++){
        for(int k=kts-1;k<=kte+1;k++){
            h_p[index4b(i,k,jts,zero)]=buff[index_buff2(i,k,0)];
            h_p[index4b(i,k,jte,zero)]=buff[index_buff2(i,k,1)];
        }}

//    #ifdef mpi_ed
        //module_halo_mp_glob_updatehalo_real_3d_(h_p,&zero,&iter_max);
        glob_updatehalo_p_(h_p,&zero,&iter_max, &its, &ite, &kts, &kte, &jts, &jte);
//    #endif  
        
        for(int j=jts;j<=jte;j++){
        for(int k=kts-1;k<=kte+1;k++){
            //buff[(k-kts+1)+(j-jts)*NY]=h_p[index4b(its,k,j,zero)];
            buff[index_buff1(0,k,j)]=h_p[index4b(its-1,k,j,zero)];
            buff[index_buff1(1,k,j)]=h_p[index4b(ite+1,k,j,zero)];
        }}
        for(int i=its-1;i<=ite+1;i++){
        for(int k=kts-1;k<=kte+1;k++){
            buff[index_buff2(i,k,0)]=h_p[index4b(i,k,jts-1,zero)];
            buff[index_buff2(i,k,1)]=h_p[index4b(i,k,jte+1,zero)];
        }}

        CHECK(hipMemcpy(d_buff,buff,2*(NX*NY+NY*NZ)*sizeof(float),hipMemcpyHostToDevice));
        hipLaunchKernelGGL(gcr_memh2d, numBlocks1, threadsPerBlock1, 0, 0,  d_p+m_tmp*NG1,d_buff, its, ite, jts, jte, kts, kte);
    	CHECK(hipDeviceSynchronize()); 

  return ;
}
