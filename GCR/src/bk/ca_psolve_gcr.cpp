#include "norm.hpp"
#include "gcr.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include "gptl.h"
#define cagcr
extern "C" void update_a(float *a_helm, int its, int ite, int jts, int jte, int kts, int kte);
//extern "C" void update_a(float *a_helm,float *a_helm_new, int its, int ite, int jts, int jte, int kts, int kte);
extern "C" void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern "C" void matrixpro_ca(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern "C" void module_halo_mp_glob_updatehalo_real_3d_(float* h_p, int *,int *); 
extern "C" void svrasr( float*  yl,int m, float*  b, int ni, int nk, int nj);
extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);

#define index_gram(i,j) (j+(i)*(step+1))
#define index_pc(i,j) (i+(j)*(2*step+1))
extern "C" void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);

double gram(float *a,float *G,float *b,int step){
    double r2=0.0;
    for(int ki=0;ki<2*step+1;ki++){
        double r1=0.0;
        for(int kj=0;kj<2*step+1;kj++){
            r1+=a[kj]*G[index_gram((kj)%(step+1),ki%(step+1))];
        }
        r2+=r1*b[ki];
    }
    return r2;
}
void  init_gram(float* __restrict__ G, float * Gm, float* __restrict__ h_ap, float* __restrict__ h_r, int NG, int step, int *mytid){
    for(int gj=0;gj<=step;gj++){
        for(int gi=0;gi<=gj;gi++){
            G[index_gram(gi,gj)]=cblas_sdot(NG,h_ap+gi*NG,1,h_ap+gj*NG,1);
            G[index_gram(gj,gi)]=G[index_gram(gi,gj)];
        }
        Gm[index_gram(0,gj)]=cblas_sdot(NG,h_r,1,h_ap+gj*NG,1);//[pi(:)]' * [ahPR(:,gj)]    
    }

    for(int gi=1;gi<=step;gi++){
        for(int gj=0;gj<=step;gj++){
            Gm[index_gram(gi,gj)]=G[index_gram(gi-1,gj)];
        }
    } 
/*
if(mytid[0]==0){
    printf("G=\n");                                                    
        for(int gi=0;gi<=step;gi++){
            for(int gj=0;gj<=step;gj++){
                printf("%4.2E  ",G[index_gram(gi,gj)]);
            }
            printf("\n");
        }
            printf("Gm=\n");                                                    
        for(int gi=0;gi<=step;gi++){
            for(int gj=0;gj<=step;gj++){
                printf("%4.2E  ",Gm[index_gram(gi,gj)]);
            }
            printf("\n");
        }
    }*/
}

extern "C" void  ca_psolve_main(double ep, float *a0, float *f0, float *cm_helm, int iter_max, double *x0, 
    int dep, int  jdep, int  ids, int  ide, int  jds, int  jde, int  kds, int  kde, 
    int ims, int  ime, int  jms, int  jme, int  kms, int  kme, 
    int its, int  ite, int  jts, int  jte, int  kts, int  kte, int *mytid){

    #include "defvar.hpp"

    initialize(h_p,NG1);
    initialize(h_ap,NG*(iter_max+1));
    initialize(h_r,NG);

    
    para_jte1(h_p[index4b(i,k,j,0)]=x0[index_x(i,k,j)];)
    matrixpro_c(a0, h_p,h_r,its, ite, jts, jte, kts, kte, jend);
    m = 0;
    d = 0.0;
    para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]-h_r[index3(i,k,j)];)
    //initialize(h_p,NG1);
    para_jte(h_p[index4b(i,k,j,0)] = h_r[index3(i,k,j)];)
    int s=0,zero=0;

#define precondition

    int step=iter_max;
        
        float *p_c=new float[(step+1)*(2*step+1)];
        float *G=new float[(step+1)*(step+1)];//={0.0};
        float *Gm=new float[(step+1)*(step+1)];//={0.0};
        float *y=new float[NG];
    while(s<max_iteration){
    //while(s<30){
    //while(s<6){
        GPTLstart ("iteration");  
#ifdef precondition
glob_updatehalo_p_(h_p, &zero, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    svrasr(h_p,zero, cm_helm, NX, NY, NZ);
#endif
        //printp(0);
        for(m=0;m<=step;m++){
#ifdef precondition
            glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
#endif
            matrixpro_c(a0, h_p+m*NG1,h_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
            para_jend(h_p[index4b(i,k,j,m+1)] = h_ap[index4(i,k,j,m)];)
int m_tmp=m+1;

#ifdef precondition
            glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
            svrasr(h_p,m_tmp, cm_helm, NX, NY, NZ);
        //printp(m_tmp);
#endif

        }
        init_gram(G, Gm, h_ap, h_r, NG, step, mytid);
        initialize(p_c,(step+1)*(2*step+1));
        float r_c[2*step+1];
        float x_c[2*step+1];
        para_step(r_c[l]=0.0;x_c[l]=0.0;)
        p_c[index_pc(0,0)]=1.0;r_c[step+1]=1.0;
//#define mpi_reduce
#ifdef mpi_reduce
        float betap[step+1];
        for(m=0;m<step;m++){
            s++;
            float D[m+2];
            D[0]=gram(r_c,Gm,p_c+m*(2*step+1),step);
            D[1]=gram(p_c+m*(2*step+1),G,p_c+m*(2*step+1),step);
            MPI_Allreduce(D,&D[2],2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
            float alpha=D[2]/D[3];
            para_step(x_c[l]+=alpha*p_c[index_pc(l,m)];)
            para_step(if(l!=0&&l!=step+1)r_c[l]+=-alpha*p_c[index_pc(l-1,m)];)
            para_step(p_c[index_pc(l,m+1)]=r_c[l];)
            for(int j=0;j<=m;j++){
                D[j]=(float)gram(r_c,G,p_c+j*(2*step+1),step);
            }   
            D[m+1]=gram(p_c+m*(2*step+1),G,p_c+m*(2*step+1),step);
            float beta2[m+2];
            MPI_Allreduce(D,beta2,m+2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
            betap[m]=beta2[m+1];
            for(int j=0;j<=m;j++){
                float beta=-beta2[j]/betap[j];
                para_step(p_c[index_pc(l,m+1)]+=beta*p_c[index_pc(l,j)];)
            }
        }
#else
        for(int m=0;m<step;m++){
            s++;
            double D=gram(r_c,Gm,p_c+m*(2*step+1),step);
            double alpha=D/gram(p_c+m*(2*step+1),G,p_c+m*(2*step+1),step);
            para_step(x_c[l]+=alpha*p_c[index_pc(l,m)];)
            para_step(if(l!=0&&l!=step+1)r_c[l]+=-alpha*p_c[index_pc(l-1,m)];)
            para_step(p_c[index_pc(l,m+1)]=r_c[l];)
            double beta;
            for(int j=0;j<=m;j++){
                beta=-gram(r_c,G,p_c+j*(2*step+1),step)/gram(p_c+j*(2*step+1),G,p_c+j*(2*step+1),step);
                para_step(p_c[index_pc(l,m+1)]+=beta*p_c[index_pc(l,j)];)
            }
        }
#endif
        if(mytid[0]==12){
                printf("x_c=\n");
                    for(int l=0;l<2*step+1;l++) printf("%f ", x_c[l]);
                    printf("\n");
        }
        para_jend(y[index3(i,k,j)]=x_c[step]*h_p[index4b(i,k,j,step)];)
        //para_jend(h_p[index3b(i,k,j)]+=(x_c[0]+x_c[step+1])*h_r[index3(i,k,j)];)
        for(int l=0;l<step;l++){
            para_jte(y[index3(i,k,j)]+=(x_c[l]+x_c[l+step+1])*h_p[index4b(i,k,j,l)];)
        } 
for(int ii=0;ii<=step;ii++){
//    	CHECK(hipDeviceSynchronize()); 
        //printp(ii)
}
        para_jte(x0[index_x(i,k,j)]+=y[index3(i,k,j)];)
//printy()



        initialize(h_p,NG1);
        //initialize(h_r,NG);
        para_jte(h_p[index4b(i,k,j,0)]=y[index3(i,k,j)];)//boundary of d_p is 0 
#ifdef precondition
        glob_updatehalo_p_(h_p, &zero, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
#endif
        matrixpro_c(a0, h_p,h_ap,its, ite, jts, jte, kts, kte, jend);
        para_jend(h_r[index3(i,k,j)] -= h_ap[index3(i,k,j)];)
        para_jte(h_p[index4b(i,k,j,0)] = h_r[index3(i,k,j)];)
        printr();
        GPTLstop ("iteration");  
        //return;
    }
}


