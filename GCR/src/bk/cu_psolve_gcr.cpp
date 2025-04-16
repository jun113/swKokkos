#define CUDA
//#define precondition
#include <math.h>
//#ifdef CUDA
#include <hip/hip_runtime.h>
#include <hipblas.h>
//#endif
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include "gptl.h"
//#include "pop_halo_c.hpp"
#include "norm.hpp"
#include "gcr.h"
#include "matrixpro.hpp"
#include "svrasr.hpp"

extern double norm2(const float *v, const  int n);
extern "C" void matrixpro_x(float *a,double *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend, int ims, int ime,  int jms, int jme);
extern "C" void update_a(float *a_helm, float *a0, int its, int ite, int jts, int jte, int kts, int kte);
//extern "C" void update_a(float *a_helm, int its, int ite, int jts, int jte, int kts, int kte);
extern "C" void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern "C" void svrasr( float*  yl,int m, float*  b, int ni, int nk, int nj);
extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);
extern "C" void cu_glob_updatehalo_p(float* h_p, float* d_p, float* h_buff, float* d_buff,
    int its, int ite, int kts, int kte, int jts, int jte, 
        int m, int iter_max, int* mytid);

__global__ void init_p(float *d_p, double *d_x0, int its, int ite, int jts, int jte, int kts, int kte,int jde) {
    #include "defindex.hpp"
    cupara_jte1(d_p[index4b(i,k,j,0)]=d_x0[index_x(i,k,j)];)
}
__global__ void trans_a(float* __restrict__ d_a_helm, float*  __restrict__ d_a0, int its, int ite, int jts, int jte, int kts, int kte,int jde) {
    #include "defindex.hpp"
    for(int m=1;m<=19;m++){
        cupara_jte(d_a0[index_a(m,i,k,j)]=d_a_helm[index_a0(m,i,k,j)];)
    }
}
__global__ void gcr_init( float* __restrict__ d_f0, float* __restrict__ d_p, float* __restrict__ d_r,int its, int ite, int jts,int jte,int kts, int kte, int jde){
	#include "defindex.hpp"
        cupara_jend(d_r[index3(i,k,j)]  = d_f0[index3(i,k,j)]-d_r[index3(i,k,j)];)
        cupara_jte(d_p[index4b(i,k,j,0)] =d_r[index3(i,k,j)];)
	}

__global__ void gcr_iteration(int iter_max,float* __restrict__ d_a0, float* __restrict__ d_p, float* __restrict__ d_ap,int its, int ite, int jts,int jte,int kts, int kte, int jde){
	#include "defindex.hpp"
        cupara_jte1(d_p[index4b(i,k,j,0)] = d_p[index4b(i,k,j,iter_max-1)];)
}

__global__ void gcr_iteration3(int m,double* __restrict__ c2, double* __restrict__ aps,float* __restrict__ d_p, int its, int ite, int jts,int jte,int kts, int kte, int jde){
		#include "defindex.hpp"
        for(int l=0;l<=m;l++){
            double bl=-c2[l]/aps[l];
            cupara_jte1( 
                    d_p[index4b(i,k,j,m+1)]=d_p[index4b(i,k,j,m+1)]+bl*d_p[index4b(i,k,j,l)];
                    )
            } 
}
__global__ void gcr_iteration2(int m,float* __restrict__ d_a0,double *d_x0, float* __restrict__ d_ar,double  d_ac, float* __restrict__ d_p, float* __restrict__ d_r, float* __restrict__ d_ap,int its, int ite, int jts,int jte,int kts, int kte, int jde){
	#include "defindex.hpp"
    cupara_jte1(
        d_x0[index_x(i,k,j)]= d_x0[index_x(i,k,j)]+d_ac*d_p[index4b(i,k,j,m)];
    )
    cupara_jend(
        d_r[index3(i,k,j)] = d_r[index3(i,k,j)]-d_ac*d_ap[index4(i,k,j,m)];
    )
    cupara_jte(d_p[index4b(i,k,j,m+1)] =d_r[index3(i,k,j)];)
}
/*
__global__ void gcr_iteration2b(int m,float* __restrict__ d_a0,double *d_x0, float* __restrict__ d_ar,double  d_ac, float* __restrict__ d_p, float* __restrict__ d_r, float* __restrict__ d_ap,int its, int ite, int jts,int jte,int kts, int kte, int jde){
	#include "defindex.hpp"
    cupara_jend(
        d_r[index3(i,k,j)] = d_r[index3(i,k,j)]-d_ac*d_ap[index4(i,k,j,m)];
    )
//    cupara_jte(d_p[index4b(i,k,j,m+1)] =d_r[index3(i,k,j)];)
}
*/
__global__ void gcr_iteration2b(int m,float* __restrict__ d_a0,double *d_x0, float* __restrict__ d_ar,double  d_ac, float* __restrict__ d_p, float* __restrict__ d_r, float* __restrict__ d_ap,int its, int ite, int jts,int jte,int kts, int kte, int jde){
	#include "defindex.hpp"
    cupara_jte(d_p[index4b(i,k,j,m+1)] =d_r[index3(i,k,j)];)
}
/*void cu_glob_updatehalo_p(float* h_p, float* d_p, float* h_buff, float* d_buff,
    int its, int ite, int kts, int kte, int jts, int jte, 
        int m, int iter_max, int* mytid);
*/
//#include "boundary_halo_fun.hpp"
//#define full_halo
#define precondition
//#define share_mem
extern "C" void module_gcr_mp_matrixpro_(float *, float *, float *, 
            int *ids, int *ide, int *jds, int *jde, int *kds, int *kde,
            int *ims, int *ime, int *jms, int *jme, int *kms, int *kme,
            int *its, int *ite, int *jts, int *jte, int *kts, int *kte);
extern  int  init_flag;
extern "C" void init_spmv(float *a_helm, float *d_a0, float *h_ap, float *d_ap, int its, int ite, int jts, int jte, int kts, int kte, int jend, int *mytid); 
//extern "C" void test_sme(int process_Rank);
//#ifdef trans
//extern "C" void  cu_psolve_main(float ep, float* __restrict__ a_helm, float* __restrict__ f0, float* __restrict__ cm_helm,int iter_max, double *x0, 
//#else
extern "C" void  cu_psolve_main(float ep, float* __restrict__ a0, float* __restrict__ f0, float* __restrict__ cm_helm,int iter_max, double *x0, 
//#endif
    int dep, int  jdep, int  ids, int  ide, int  jds, int  jde, int  kds, int  kde, 
    int ims, int  ime, int  jms, int  jme, int  kms, int  kme, 
    int its, int  ite, int  jts, int  jte, int  kts, int  kte, int *mytid){
    //FILE * fp=fopen("data.txt","w");
    //FILE * fp0=fopen("data0.txt","w");
    //fprintf(fp,"d=\n");
    GPTLstart ("defvar");
    #include "defvar.hpp"
    GPTLstop ("defvar");
    GPTLstart ("initialization");
/*#ifdef share_mem
    dim3 threadsPerBlock_mat(16, 4, 6);
    dim3 numBlocks_mat(ceil(((float)NX)/(threadsPerBlock_mat.x-2)),ceil(((float) NY-1)/(threadsPerBlock_mat.y)),ceil(((float)NZ)/(threadsPerBlock_mat.z-2)));
#else
    dim3 threadsPerBlock_mat(16, 4, 4);
    //dim3 threadsPerBlock(10, 6, 6);
    dim3 numBlocks_mat(ceil(((float)NX)/16),ceil(((float) NY)/4),ceil(((float)NZ)/4));
    //dim3 numBlocks(1,1,1);
#endif*/
    dim3 threadsPerBlock(64, 4, 1);
    dim3 numBlocks(ceil(((float)NX)/64),ceil(((float) NY)/4),NZ);
//hipLaunchKernelGGL(cu_matrixpro, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_p+m*NG1,d_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
    dim3 threadsPerBlock2(64, 1, 1); 
    dim3 numBlocks2(1,NZ,1); //halo=1

    dim3 block_size(64,4,1);
    dim3 num_blocks((NX + block_size.x - 1) / block_size.x,(NY + block_size.y - 1) / block_size.y,NZ);
    //GPTLstart ("cudamemcpy");
        //glob_updatehalo_x0_(x0, &ims, &ime, &jms, &jme, &kts, &kte);
	CHECK(hipMemcpy(d_x0,x0,NG2*sizeof(double),hipMemcpyHostToDevice));
	CHECK(hipMemset(d_p,0,NG1*(iter_max+1)*sizeof(float)));
       
/*
    GPTLstart ("trans_a");
    	CHECK(hipDeviceSynchronize()); //test time
        hipLaunchKernelGGL(trans_a, numBlocks_mat, threadsPerBlock_mat , 0, 0, d_a_helm, d_a0, its, ite, jts, jte, kts, kte, jde);
    	CHECK(hipDeviceSynchronize()); //test time
    GPTLstop ("trans_a");
*/
#ifdef trans
    GPTLstart ("trans_a");
    	//CHECK(hipDeviceSynchronize()); //test time
        hipLaunchKernelGGL(trans_a, num_blocks, block_size, 0, 0, d_a_helm, d_a0, its, ite, jts, jte, kts, kte, jde);
    	CHECK(hipDeviceSynchronize()); //test time // a_helm -> a0
        //initialize(h_p,NG1);
    GPTLstop ("trans_a");
#endif
 //hipLaunchKernelGGL(init_p, num_blocks, block_size, 0, 0, d_p, d_x0, its, ite, jts, jte, kts, kte, jde);
// hipLaunchKernelGGL(cu_matrixpro, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_p,d_r, its, ite, jts, jte, kts, kte, jend);
    hipLaunchKernelGGL(cu_matrixpro_x, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_x0,d_r, its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
   /* 
    cuprintr();
	CHECK(hipMemcpy(h_ap,d_r,NG*sizeof(float),hipMemcpyDeviceToHost));
        init_spmv(a0, d_a0, h_ap, d_p,  its,  ite,  jts,  jte,  kts,  kte,  jend,  mytid); 
        return;
*/
    hipLaunchKernelGGL(gcr_init, numBlocks, threadsPerBlock, 0, 0, d_f0, d_p,d_r, its, ite, jts, jte, kts, kte, jde);
    int zero=0;
    int m_tmp=0;
    //CHECK(hipDeviceSynchronize()); //must sync d_p
    cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
    hipLaunchKernelGGL(cu_svrasr, numBlocks2, threadsPerBlock2, 0, 0, d_p, 0,d_cm_helm, NX, NY, NZ);
    //CHECK(hipDeviceSynchronize()); //must sync d_p
    cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;


    GPTLstop ("initialization");
//cuprintr();
    int s;
    for(s=0;s<max_iteration;s++){
        GPTLstart ("cu_iteration");  
        int m = s%(iter_max-1);
        m_tmp=m+1;
        if ((m == 0) &&( s!=0) ) {
            hipLaunchKernelGGL(gcr_iteration, numBlocks, threadsPerBlock, 0, 0, iter_max,d_a0, d_p,d_ap, its, ite, jts, jte, kts, kte, jde);
        //    cuinitialize(d_r,NG);
        //hipLaunchKernelGGL(init_p, num_blocks, block_size, 0, 0, d_p, d_x0, its, ite, jts, jte, kts, kte, jde);
        //hipLaunchKernelGGL(cu_matrixpro, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_p,d_r, its, ite, jts, jte, kts, kte, jend);
        //hipLaunchKernelGGL(gcr_init, numBlocks, threadsPerBlock, 0, 0, d_f0, d_p,d_r, its, ite, jts, jte, kts, kte, jde);
        //cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;

        }
    	CHECK(hipDeviceSynchronize()); //sync d_p
            GPTLstart ("matrixpro");
        hipLaunchKernelGGL(cu_matrixpro, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_p+m*NG1,d_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
    	CHECK(hipDeviceSynchronize()); //test time
            GPTLstop ("matrixpro");
        hipblasSdot(h, NG, d_r, 1, d_r, 1, dd);
        CHECK(hipMemcpy(&d_r4, dd, sizeof(float), hipMemcpyDeviceToHost));
            GPTLstart ("reduction");
        d=d_r4;
        hipblasSdot(h, NG, d_r, 1, d_ap+m*NG, 1, result);
        hipblasSdot(h, NG, d_ap+m*NG, 1, d_ap+m*NG, 1, result+1);
        CHECK(hipMemcpy(c1_r4, result, 2*sizeof(float), hipMemcpyDeviceToHost));
        //CHECK(hipMemcpy(&c1_r4[1], result+1, sizeof(float), hipMemcpyDeviceToHost));
            GPTLstop ("reduction");
        c1[0]=c1_r4[0];
        c1[1]=c1_r4[1];
        MPI_Allreduce(c1,c2,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        double ac=c2[0]/c2[1];
        aps[m]=c2[1];
 
#define origin2
#ifdef origin2
        hipLaunchKernelGGL(gcr_iteration2, numBlocks, threadsPerBlock, 0, 0, m, d_a0, d_x0,d_ar, ac, d_p, d_r, d_ap, its, ite, jts, jte, kts, kte, jde);
//    hipLaunchKernelGGL(cu_matrixpro, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_p1,d_ap, its, ite, jts, jte, kts, kte, jend);
//    hipLaunchKernelGGL(gcr_final, numBlocks, threadsPerBlock, 0, 0,  d_x0, d_ap, d_r, d_p1,  its, ite, jts, jte, kts, kte, jend); 
#else
        hipLaunchKernelGGL(gcr_iteration2, numBlocks, threadsPerBlock, 0, 0, m, d_a0, d_x0,d_ar, ac, d_p, d_r, d_ap, its, ite, jts, jte, kts, kte, jde);
        hipLaunchKernelGGL(cu_matrixpro_x, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_x0,d_r, its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
    hipLaunchKernelGGL(gcr_init, numBlocks, threadsPerBlock, 0, 0, d_f0, d_p+(m+1)*NG1,d_r, its, ite, jts, jte, kts, kte, jde);
//        hipLaunchKernelGGL(gcr_iteration2b, numBlocks, threadsPerBlock, 0, 0, m, d_a0, d_x0,d_ar, ac, d_p, d_r, d_ap, its, ite, jts, jte, kts, kte, jde);
#endif

        //int two=2;
        //CHECK(hipMemcpy(h_p+2*NG1,d_p+(m+1)*NG1,NG1*sizeof(float),hipMemcpyDeviceToHost));
        //glob_updatehalo_p_(h_p,&two,&iter_max, &its, &ite, &kts, &kte, &jts, &jte);
        //svrasr(h_p,two, cm_helm, NX, NY, NZ); 
        CHECK(hipDeviceSynchronize()); //must sync d_r and d_p
        hipblasSdot(h, NG, d_r, 1, d_r, 1, dd);
        CHECK(hipDeviceSynchronize()); //test time


            GPTLstart ("cu_glob_updatehalo");
        cu_glob_updatehalo_p(h_p, d_p+(m+1)*NG1, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
            GPTLstop ("cu_glob_updatehalo");
            GPTLstart ("precondition");
        hipLaunchKernelGGL(cu_svrasr, numBlocks2, threadsPerBlock2, 0, 0, d_p+(m+1)*NG1, 0,d_cm_helm, NX, NY, NZ);
        CHECK(hipDeviceSynchronize()); //must sync d_p
            GPTLstop ("precondition");
    //    CHECK(hipMemcpy(h_p,d_p+(m+1)*NG1,NG1*sizeof(float),hipMemcpyDeviceToHost));
    //    print_diff();
        cu_glob_updatehalo_p(h_p, d_p+(m+1)*NG1, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
         CHECK(hipMemcpy(&d_r4, dd, sizeof(float), hipMemcpyDeviceToHost));
        d=d_r4;
 
//        CHECK(hipDeviceSynchronize()); //must sync d_p
            GPTLstart ("matrixpro");
        hipLaunchKernelGGL(cu_matrixpro, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_p+(m+1)*NG1,d_ar, its, ite, jts, jte, kts, kte,jend);
        CHECK(hipDeviceSynchronize()); //test time
            GPTLstop ("matrixpro");


            GPTLstart ("reduction");
        c1[m+1] = d;
        for(l=0;l<=m;l++){
            hipblasSdot(h, NG, d_ar, 1, d_ap+l*NG, 1, result+l);
        }
            CHECK(hipMemcpy(c1_r4, result, (m+1)*sizeof(float), hipMemcpyDeviceToHost));
        for(l=0;l<=m;l++){
            c1[l]=c1_r4[l];
        }  
            GPTLstop ("reduction");
   
        MPI_Allreduce(c1,c2,m+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        d = c2[m+1];
        CHECK(hipMemcpy(d_c1,c2,(m+1)* sizeof(double), hipMemcpyHostToDevice));
        CHECK(hipMemcpy(d_aps,aps,(m+1)* sizeof(double), hipMemcpyHostToDevice));
        hipLaunchKernelGGL(gcr_iteration3, numBlocks, threadsPerBlock, 0, 0, m,d_c1,d_aps, d_p, its, ite, jts, jte, kts, kte, jde);
        //if(mytid[0]==0)printf("its=%d,cu_res=%E\n",s,d);
        GPTLstop ("cu_iteration");  
        if (d <= ep || s == max_iteration-1 ) { 
            if(mytid[0]==0)printf(" RES of gcr  %E in %d iterations\n", sqrt(d) ,s);
            break;}
         if(mytid[0]==0)printf(" RES of gcr  %E in %d iterations\n", sqrt(d) ,s);
   //MPI_Barrier(MPI_COMM_WORLD);
        //return;
	//	if(m==1)break;
        /*
            GPTLstart ("cudamemcpy_x0");
	CHECK(hipMemcpy(x0,d_x0,NG2*sizeof(double),hipMemcpyDeviceToHost));
    GPTLstop ("cudamemcpy_x0");
        matrixpro_x(a0, x0,h_ap,its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
        initialize(h_r,NG);
        para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]- h_ap[index3(i,k,j)];)

    c1[0] = norm2(h_r, NG);
    MPI_Allreduce(c1,c2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    c2[0]=sqrt(c2[0]);
    if(mytid[0]==12)
        printf("RES of GCR is %e \n",c2[0]);
*/

   }
    GPTLstart ("cudamemcpy_x0");
	CHECK(hipMemcpy(x0,d_x0,NG2*sizeof(double),hipMemcpyDeviceToHost));
    GPTLstop ("cudamemcpy_x0");
    /*
        matrixpro_x(a0, x0,h_ap,its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
        initialize(h_r,NG);
        para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]- h_ap[index3(i,k,j)];)

        hipLaunchKernelGGL(cu_matrixpro_x, numBlocks_mat, threadsPerBlock_mat, 0, 0, d_a0, d_x0,d_r, its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
    hipLaunchKernelGGL(gcr_init, numBlocks, threadsPerBlock, 0, 0, d_f0, d_p,d_r, its, ite, jts, jte, kts, kte, jde);
	CHECK(hipMemcpy(h_r,d_r,NG*sizeof(float),hipMemcpyDeviceToHost));
    c1[0] = norm2(h_r, NG);
    MPI_Allreduce(c1,c2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    c2[0]=sqrt(c2[0]);
    if(mytid[0]==5|| mytid[0]==12)
        printf("RES of gmres is %e \n",c2[0]);
*/


}
