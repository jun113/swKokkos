#include "hip/hip_runtime.h"
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
    cupara_jend(
        d_x0[index_x(i,k,j)]= d_x0[index_x(i,k,j)]+d_ac*d_p[index4b(i,k,j,m)];
        d_r[index3(i,k,j)] = d_r[index3(i,k,j)]-d_ac*d_ap[index4(i,k,j,m)];
    )
    cupara_jte(d_p[index4b(i,k,j,m+1)] =d_r[index3(i,k,j)];)
}
void cu_glob_updatehalo_p(float* h_p, float* d_p, float* h_buff, float* d_buff,
    int its, int ite, int kts, int kte, int jts, int jte, 
        int m, int iter_max, int* mytid);
//#define precondition
//#define full_halo
#include "boundary_halo_fun.hpp"
void  cu_psolve_main(float ep, float* __restrict__ a0, float* __restrict__ f0, float* __restrict__ cm_helm,int iter_max, double *x0, 
    int dep, int  jdep, int  ids, int  ide, int  jds, int  jde, int  kds, int  kde, 
    int ims, int  ime, int  jms, int  jme, int  kms, int  kme, 
    int its, int  ite, int  jts, int  jte, int  kts, int  kte, int *mytid){

    #include "defvar.hpp"
    dim3 threadsPerBlock(8, 4, 4);
    dim3 numBlocks(ceil(((float)NX)/8),ceil(((float) NY)/4),ceil(((float)NZ)/4));
    dim3 threadsPerBlock1(8, 4, 1);
    dim3 numBlocks1(ceil(((float)NX)/8)+ceil(((float)NZ)/8),ceil(((float) NY)/4),1);
    dim3 threadsPerBlock2(32, 1, 1); //is fatter than (32,4,1)
    dim3 numBlocks2(1,NZ,1); //halo=1
    GPTLstart ("cudamemcpy");
	CHECK(hipMemcpy(d_x0,x0,NG2*sizeof(double),hipMemcpyHostToDevice));
	CHECK(hipMemset(d_p,0,NG1*(iter_max+1)*sizeof(float)));
        //initialize(h_p,NG1);
    GPTLstop ("cudamemcpy");

//cuprintp(0);
    hipLaunchKernelGGL(matrixpro_sme_x, numBlocks, threadsPerBlock, 0, 0, d_a0, d_x0,d_r, its, ite, jts, jte, kts, kte, jde);
//    hipLaunchKernelGGL(matrixpro_sme, numBlocks, threadsPerBlock, 0, 0, d_a0, d_p,d_r, its, ite, jts, jte, kts, kte, jde);

    hipLaunchKernelGGL(gcr_init, numBlocks, threadsPerBlock, 0, 0, d_f0, d_p,d_r, its, ite, jts, jte, kts, kte, jde);
    int zero=0;
#ifdef full_halo
    CHECK(hipMemcpy(h_p,d_p,NG1*sizeof(float),hipMemcpyDeviceToHost));
    glob_updatehalo_p_(h_p,&zero,&iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    #ifdef precondition
    svrasr(h_p,zero, cm_helm, NX, NY, NZ); 
    glob_updatehalo_p_(h_p,&zero,&iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    #endif
    CHECK(hipMemcpy(d_p,h_p,NG1*sizeof(float),hipMemcpyHostToDevice));
#else
    cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
    #ifdef precondition
    hipLaunchKernelGGL(cu_svrasr, numBlocks2, threadsPerBlock2, 0, 0, d_p, 0,d_cm_helm, NX, NY, NZ);
    cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
    #endif
#endif


//cuprintr();
    int s;
    for(s=0;s<max_iteration;s++){
        GPTLstart ("cu_iteration");  
        int m = s%(iter_max-1);
        if ((m == 0) &&( s!=0) ) {
            hipLaunchKernelGGL(gcr_iteration, numBlocks, threadsPerBlock, 0, 0, iter_max,d_a0, d_p,d_ap, its, ite, jts, jte, kts, kte, jde);
        }
    	CHECK(hipDeviceSynchronize()); 
        hipLaunchKernelGGL(matrixpro_sme, numBlocks, threadsPerBlock, 0, 0, d_a0, d_p+m*NG1,d_ap+m*NG, its, ite, jts, jte, kts, kte, jde);
    	CHECK(hipDeviceSynchronize()); 
//cuprintp(m);
//cuprintap(m);
        //initialize(d_c1,iter_max+4);
        //initialize(d_c2,iter_max+4);
        hipblasSdot(h, NG, d_r, 1, d_r, 1, dd);
        hipblasSdot(h, NG, d_r, 1, d_ap+m*NG, 1, result);
        hipblasSdot(h, NG, d_ap+m*NG, 1, d_ap+m*NG, 1, result+1);
        CHECK(hipMemcpy(&d_r4, dd, sizeof(float), hipMemcpyDeviceToHost));
        d=d_r4;
        CHECK(hipMemcpy(c1_r4, result, sizeof(float), hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(&c1_r4[1], result+1, sizeof(float), hipMemcpyDeviceToHost));
        c1[0]=c1_r4[0];
        c1[1]=c1_r4[1];
        //ac[m]=c1[0]/c1[1];
        //aps[m]=c1[1];
        MPI_Allreduce(c1,c2,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        double ac=c2[0]/c2[1];
        aps[m]=c2[1];
//        if(mytid[0]==0)printf("line=%d,c1=%E,%E,c2=%E,%E,ac=%E,r_ap=%f\n",__LINE__,c1[0],c1[1],c2[0],c2[1],ac[m],c12_tmp);
        //if(mytid[0]==0)printf("ac=%E,aps=%E\n",__LINE__,ac,aps[m]);
 
        hipLaunchKernelGGL(gcr_iteration2, numBlocks, threadsPerBlock, 0, 0, m, d_a0, d_x0,d_ar, ac, d_p, d_r, d_ap, its, ite, jts, jte, kts, kte, jde);
        hipblasSdot(h, NG, d_r, 1, d_r, 1, dd);
        //CHECK(hipMemcpy(h_p,d_p,(m_tmp+1)*NG1*sizeof(float),hipMemcpyDeviceToHost));
        //module_halo_mp_glob_updatehalo_real_3d_(h_p,&m_tmp,&iter_max);
        
#ifdef full_halo
        CHECK(hipMemcpy(h_p,d_p+(m+1)*NG1,NG1*sizeof(float),hipMemcpyDeviceToHost));
        glob_updatehalo_p_(h_p,&zero,&iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    #ifdef precondition
        svrasr(h_p,zero, cm_helm, NX, NY, NZ); 
        glob_updatehalo_p_(h_p,&zero,&iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    #endif
        CHECK(hipMemcpy(d_p+(m+1)*NG1,h_p,NG1*sizeof(float),hipMemcpyHostToDevice));
#else
        //cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte, m+1, iter_max,  mytid) ;
        CHECK(hipDeviceSynchronize()); //must sync d_p
    GPTLstart ("cu_glob_updatehalo");
#include "boundary_halo.hpp"
//        cu_glob_updatehalo_p(h_p, d_p+(m+1)*NG1, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
    GPTLstop ("cu_glob_updatehalo");
    #ifdef precondition
    GPTLstart ("cu_svrasr");
        hipLaunchKernelGGL(cu_svrasr, numBlocks2, threadsPerBlock2, 0, 0, d_p, m+1,d_cm_helm, NX, NY, NZ);
        CHECK(hipDeviceSynchronize()); //must sync d_p
    GPTLstop ("cu_svrasr");
        //cu_glob_updatehalo_p(h_p, d_p, buff, d_buff,its, ite, kts, kte, jts, jte,  m+1, iter_max,  mytid) ;
    GPTLstart ("cu_glob_updatehalo");
        cu_glob_updatehalo_p(h_p, d_p+(m+1)*NG1, buff, d_buff,its, ite, kts, kte, jts, jte, zero, iter_max,  mytid) ;
    GPTLstop ("cu_glob_updatehalo");
    #endif
#endif
         CHECK(hipMemcpy(&d_r4, dd, sizeof(float), hipMemcpyDeviceToHost));
        d=d_r4;
 
        //CHECK(hipDeviceSynchronize()); //must sync d_p
        hipLaunchKernelGGL(matrixpro_sme, numBlocks, threadsPerBlock, 0, 0, d_a0, d_p+(m+1)*NG1,d_ar, its, ite, jts, jte, kts, kte,jde);
        c1[m+1] = d;
        for(l=0;l<=m;l++){
            hipblasSdot(h, NG, d_ar, 1, d_ap+l*NG, 1, result+l);
            CHECK(hipMemcpy(&c1_r4[l], result+l, sizeof(float), hipMemcpyDeviceToHost));
            c1[l]=c1_r4[l];
        }  
   
        MPI_Allreduce(c1,c2,m+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        d = c2[m+1];
    	//if(mytid[0]==0){printf("line=%d,its=%d\n",__LINE__,s);}
        if (d <= ep || s == max_iteration-1 ) { 
            if(mytid[0]==0)printf(" RES of gcr  %E in %d iterations\n", d ,s);
            break;}
        CHECK(hipMemcpy(d_c1,c2,(m+1)* sizeof(double), hipMemcpyHostToDevice));
        CHECK(hipMemcpy(d_aps,aps,(m+1)* sizeof(double), hipMemcpyHostToDevice));
        hipLaunchKernelGGL(gcr_iteration3, numBlocks, threadsPerBlock, 0, 0, m,d_c1,d_aps, d_p, its, ite, jts, jte, kts, kte, jde);
        //para_jte1( d_p[index4b(i,k,j,m+1)]=d_p[index4b(i,k,j,m+1)]+bl*d_p[index4b(i,k,j,l)];)
//cuprintp(m+1);
        //if(mytid[0]==0)printf("its=%d,cu_res=%E\n",s,d);
         if(mytid[0]==0)printf(" RES of gcr  %E in %d iterations\n", d ,s);
        GPTLstop ("cu_iteration");  
//		if(m==5)break;
   }
	CHECK(hipMemcpy(x0,d_x0,NG2*sizeof(double),hipMemcpyDeviceToHost));
	//CHECK(hipFree(d_p));CHECK(hipFree(d_a0));CHECK(hipFree(d_f0));
    //CHECK(hipFree(d_r));CHECK(hipFree(d_ap));CHECK(hipFree(d_ar));
	//CHECK(hipFree(d_x0));CHECK(hipFree(d_ac));CHECK(hipFree(d_aps));
}
