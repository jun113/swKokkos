#include "norm.hpp"
#include "gcr.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
//#include <mkl.h>
//#include <mkl_cblas.h>
//#include "gptl.h"
//#include "pop_halo_c.hpp"
//#include "matrixpro.hpp"

extern "C" void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern "C" void svrasr( float*  yl,int m, float*  b, int ni, int nk, int nj);
extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);
extern "C" void matrixpro_x(float *a,double *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend, int ims, int ime,  int jms, int jme);


double norm2(const float *v, const  int n)
{
    double tmp_sum = 0;
    #pragma omp parallel for  reduction(+:tmp_sum) schedule(runtime)
    for(int i=0; i<n; i++) 
        tmp_sum += v[i] * v[i];
    return tmp_sum;//sqrt(tmp);
}

extern "C" void  psolve_main(float ep, float *a0, float *f0, float *cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid){

        int  idep, jdep, ids,ide,jds,jde,kds,kde,ims,ime,
                jms,jme,kms,kme, its,ite,jts,jte,kts,kte;
        assign_indices(ijk_index, idep, jdep, ids, ide, jds, jde, kds, kde,
            ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte);
    #include "defvar.hpp"

    //#pragma omp simd
    initialize(h_p,NG1*(iter_max+1));
    //#pragma omp simd
    initialize(h_ap,NG*(iter_max+1));
    initialize(h_r,NG);

    //#pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
    para_jte1(h_p[index4b(i,k,j,0)]=x0[index_x(i,k,j)];)

    matrixpro_c(a0, h_p,h_r,its, ite, jts, jte, kts, kte, jend);
    m = 0;
    d = 0.0;
    // gcr_init
    //#pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
    para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]-h_r[index3(i,k,j)];)
    printr();
    //#pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
    para_jte(h_p[index4b(i,k,j,m)] = h_r[index3(i,k,j)];)

//TODO
    glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    
    svrasr(h_p,m, cm_helm, NX, NY, NZ);
    glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    d = 0.0;
    //#pragma omp parallel for private(j,k,i) reduction(+:d) schedule(runtime) 
    for(j=jbegin;j<=jend;j++){
        for(k=kts-1;k<=kte+1;k++){
            for(i=ibegin;i<=iend;i++){
                d = d + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
            }
        }
    }
    c1[0] = d;
    MPI_Allreduce(c1,c2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    d = c2[0];
    //printf("d=%.3e\n",d);
    int s;
    for(s=0;s<max_iteration;s++){
        //GPTLstart ("gcr_iteration");  
        m = s%(iter_max-1);
        if ((m == 0) &&( s!=0) ) {
            para_jte1(h_p[index4b(i,k,j,0)] = h_p[index4b(i,k,j,iter_max-1)];)
        }
        //GPTLstart ("maxtrixpro");  
        matrixpro_c(a0, h_p+m*NG1,h_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
        //GPTLstop ("maxtrixpro");  
        double c11=0.0,c12=0.0;
        
        //#pragma omp parallel for private(j,k,i) reduction(+:c11,c12) schedule(runtime) 
        para_jend(
            c11 = c11 +  h_r[index3(i,k,j)]* h_ap[index4(i,k,j,m)];
            c12 = c12 + h_ap[index4(i,k,j,m)]* h_ap[index4(i,k,j,m)];
        )
        //GPTLstop ("reduction");  
            
        c1[0]=c11;
        c1[1]=c12;
        MPI_Allreduce(c1,c2,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        double ac=c2[0]/c2[1];
        aps[m]=c2[1];
  
        //if(mytid[0]==0) printf("ac=%.3e\n",ac);
        //#pragma omp parallel for private(j,k,i) schedule(runtime) 
        para_jte1(
            x0[index_x(i,k,j)]= x0[index_x(i,k,j)]+ac * h_p[index4b(i,k,j,m)];
        )
        para_jend(
            h_r[index3(i,k,j)] = h_r[index3(i,k,j)]-ac * h_ap[index4(i,k,j,m)];
        )
        d = 0.0;
        //#pragma omp parallel for private(j,k,i) reduction(+:d) schedule(runtime) 
        para_jend(
            d = d  + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
        )
        //#pragma omp parallel for private(j,k,i) schedule(runtime) 
        para_jte(h_p[index4b(i,k,j,m+1)] = h_r[index3(i,k,j)];)

        int m_tmp=m+1;
        //GPTLstart ("glob_updatehalo");  
        glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
        //GPTLstop ("glob_updatehalo");  
        
        //GPTLstart ("precondition");  
        svrasr(h_p,m_tmp, cm_helm, NX, NY, NZ);
        //GPTLstop ("precondition");  
        //GPTLstart ("glob_updatehalo");  
        glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
        //GPTLstop ("glob_updatehalo");  
        //GPTLstart ("maxtrixpro");  
        matrixpro_c(a0, h_p+m_tmp*NG1,h_ar, its, ite, jts, jte, kts, kte, jend);
        //GPTLstop ("maxtrixpro");  
        //GPTLstart ("reduction_cl");  
        for(l=0;l<=m;l++){
            double cl=0.0;
            //#pragma omp parallel for private(j,k,i) reduction(+:cl) schedule(runtime) 
            para_jend(cl = cl + h_ar[index3(i,k,j)]*h_ap[index4(i,k,j,l)];)
            c1[l]=cl;
        }  
        //GPTLstop ("reduction_cl");  

        c1[m+1] = d;
        MPI_Allreduce(c1,c2,m+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        d = c2[m+1];

        for(l=0;l<=m;l++){
             b[l]=-c2[l]/aps[l];
        }
        //#pragma omp parallel for private(j,k,i,l) schedule(runtime)   
        for(j=jts-1;j<=jte+1;j++){
        for(l=0;l<=m;l++){
            for(k=kts-1;k<=kte+1;k++){
                for(i=its-1;i<=ite+1;i++){
                    h_p[index4b(i,k,j,m+1)]+=b[l]*h_p[index4b(i,k,j,l)];//)
                }}}
        } 
        //GPTLstop ("reductionp");  
        //GPTLstop ("gcr_iteration");  
        if (d <= ep || s == max_iteration-1 ) { 
            if(mytid[0]==0)printf("RES of gcr is %e in %d iterations\n", sqrt(d) ,s);
            break;
        }
        if(mytid[0]==0)printf("RES of gcr is %e in %d iterations\n", sqrt(d) ,s);
        if(s==10)return;
    }
        matrixpro_x(a0, x0,h_ap,its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
        initialize(h_r,NG);
        para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]- h_ap[index3(i,k,j)];)

    c1[0] = norm2(h_r, NG);
    MPI_Allreduce(c1,c2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    c2[0]=sqrt(c2[0]);
    if(mytid[0]==5|| mytid[0]==12)
        printf("RES of gmres is %e \n",c2[0]);

	//delete[] h_p;
    //delete[] h_ap; 
   // delete[] h_ar; delete[] h_r;
}

