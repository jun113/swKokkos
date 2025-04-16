#include "gcr.h"

#include "athread_psolve_gcr.h"

#include "athread.h"
#include <mpi.h>

#include <math.h>

static float cblas_sdot(int N, const float *X, int incX, const float *Y, int incY) {
    float result = 0.0f;
    int ix = 0;
    int iy = 0;
    for (int i = 0; i < N; i++) {
        result += X[ix] * Y[iy];
        ix += incX;
        iy += incY;
    }
    return result;
}

#define para_jte(status)\
        for(j=jts;j<=jte;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its;i<=ite;i++){\
                    status\
                }\
            }\
        }

#define para_jte1(status)\
        for(j=jts-1;j<=jte+1;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its-1;i<=ite+1;i++){\
                    status\
                }\
            }\
        }

#define para_jend(status)\
        for(j=jts;j<=jend;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its;i<=ite;i++){\
                    status\
                }\
            }\
        }

#define printr()\
                c11_tmp=cblas_sdot(NG,h_r,1,h_r,1);\
       MPI_Allreduce(&c11_tmp,&c12_tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);\
                if(mytid[0]==12) printf("mytid=%d,r_norm=%E,tol=%E\n",mytid[0],sqrt(c12_tmp),sqrt(c12_tmp));

extern void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern void svrasr( float*  yl,int m, float*  b, int ni, int nk, int nj);
extern void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);
extern void matrixpro_x(float *a,double *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend, int ims, int ime,  int jms, int jme);


// inline double norm2(const float *v, const  int n)
// {
//     double tmp_sum = 0;
//     #pragma omp parallel for  reduction(+:tmp_sum) schedule(runtime)
//     for(int i=0; i<n; i++) 
//         tmp_sum += v[i] * v[i];
//     return tmp_sum;//sqrt(tmp);
// }

static void assign_indices_c (int* ind, int* idep, int* jdep, int* ids, int* ide, int* jds, int* jde, int* kds, int* kde,
                    int* ims, int* ime, int* jms, int* jme, int* kms, int* kme,
                    int* its, int* ite, int* jts, int* jte, int* kts, int* kte) {
    *idep = ind[0]; 
    *jdep = ind[1];
    *ids = ind[2];
    *ide = ind[3];
    *jds = ind[4];
    *jde = ind[5];
    *kds = ind[6];
    *kde = ind[7];
    *ims = ind[8];
    *ime = ind[9];
    *jms = ind[10];
    *jme = ind[11];
    *kms = ind[12];
    *kme = ind[13];
    *its = ind[14];
    *ite = ind[15];
    *jts = ind[16];
    *jte = ind[17];
    *kts = ind[18];
    *kte = ind[19];
    
//if(mytid[0]==0){
//    printf("kms: %d, kme: %d, its: %d, ite: %d\n", kms, kme, its, ite);
//}
}

void athread_psolve_main (float ep, float *a0, float *f0, float *cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid) {

  int idep, jdep, ids,ide,jds,jde,kds,kde,ims,ime,
      jms,jme,kms,kme, its,ite,jts,jte,kts,kte;
  assign_indices_c (ijk_index, &idep, &jdep, &ids, &ide, &jds, &jde, &kds, &kde,
      &ims, &ime, &jms, &jme, &kms, &kme, &its, &ite, &jts, &jte, &kts, &kte);

  int i,j,k,m,l;
  double d;
  // float d_r4;
  float c11_tmp;
  float c12_tmp;
 
  int max_iteration = 50;
  int ibegin = its;
  int iend   = ite;
  int jbegin = jts;
  // int jend   = min (jde - 1, jte);
  int jend = (jde - 1) < jte ? (jde - 1) : jte;
 
  int NX = ite - its + 3;
  int NY = kte - kts + 3;
  int NZ = jte - jts + 3;
  int NG1 = NX * NY * NZ;
  int NG = (ite - its + 1 )  * NY * (jte - jts + 1);
  int NG2 = (ime - ims + 1 ) * NY * (jme - jms + 1);
	
  float* h_p  = (float*) malloc (sizeof(float) * NG1 * iter_max);
  float* h_ap = (float*) malloc (sizeof(float) * NG * (iter_max + 1));
  float* h_ar = (float*) malloc (sizeof(float) * NG);
  float* h_r  = (float*) malloc (sizeof(float) * NG);

  double* aps = (double*) malloc (sizeof(double) * (iter_max + 1));
  double* b   = (double*) malloc (sizeof(double) * (iter_max + 1));
  double* c1  = (double*) malloc (sizeof(double) * (iter_max + 2));
  double* c2  = (double*) malloc (sizeof(double) * (iter_max + 2));

  memset (h_p,  0, sizeof(float) * NG1 * iter_max);
  memset (h_ap, 0, sizeof(float) * NG * (iter_max+1));
  memset (h_ar, 0, sizeof(float) * NG);
  memset (h_r,  0, sizeof(float) * NG);

  // para_jte1(h_p[index4b(i,k,j,0)]=x0[index_x(i,k,j)];)
  athread_kernel1_host (its, kts, jts, ite, kte, jte, ims, ime, jms, 
      NX, NY, NZ, h_p, x0);

  // matrixpro_c(a0, h_p,h_r,its, ite, jts, jte, kts, kte, jend);
  athread_matrixpro_host (a0, h_p, h_r, its, ite, jts, jte, kts, kte, jend);
  d = 0.0;
  // gcr_init
  m = 0;
  // para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]-h_r[index3(i,k,j)];)
  athread_kernel2_host (its, kts, jts, ite, kte, jte, ims, ime, jms, 
      jend, NX, NY, NZ, h_r, f0);

  printr();
  // para_jte(h_p[index4b(i,k,j,m)] = h_r[index3(i,k,j)];)
  athread_kernel3_host (m, its, kts, jts, ite, kte, jte, ims, ime, jms, 
      jend, NX, NY, NZ, h_p, h_r);

  //TODO
  glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    
  // svrasr(h_p,m, cm_helm, NX, NY, NZ);
  athread_svrasr_host (h_p, cm_helm, m, NX, NY, NZ);

  glob_updatehalo_p_ (h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
  d = 0.0;
  //#pragma omp parallel for private(j,k,i) reduction(+:d) schedule(runtime) 
  // for(j=jbegin;j<=jend;j++){
  //   for(k=kts-1;k<=kte+1;k++){
  //     for(i=ibegin;i<=iend;i++){
  //       d = d + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
  //     }
  //   }
  // }
  athread_kernel_reduce_d_host (its, kts, jts, ite, kte, jte,
                                ims, ime, jms, jend, NX, NY, NZ, 
                                h_r, &d);
  c1[0] = d;
  MPI_Allreduce (c1, c2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  d = c2[0];
  int s;

  for (s = 0; s < max_iteration; s++) {
    //GPTLstart ("gcr_iteration");  
    m = s % (iter_max - 1);
    if ((m == 0) && (s != 0)) {
      // para_jte1(h_p[index4b(i,k,j,0)] = h_p[index4b(i,k,j,iter_max-1)];)
      athread_kernel4_host (iter_max, its, kts, jts, ite, kte, jte, NX, NY, NZ, h_p);
    }
    //GPTLstart ("maxtrixpro");  
    // matrixpro_c(a0, h_p+m*NG1,h_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
    athread_matrixpro_host (a0, h_p+m*NG1, h_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
    //GPTLstop ("maxtrixpro");  
    double c11 = 0.0;
    double c12 = 0.0;
        
    //#pragma omp parallel for private(j,k,i) reduction(+:c11,c12) schedule(runtime) 
    // para_jend(
    //   c11 = c11 +  h_r[index3(i,k,j)]* h_ap[index4(i,k,j,m)];
    //   c12 = c12 + h_ap[index4(i,k,j,m)]* h_ap[index4(i,k,j,m)];
    // )
    athread_kernel_reduce_c11_host (m, its, kts, jts, ite, kte, jte,
        ims, ime, jms, jend, NX, NY, NZ, h_r, h_ap, &c11);
    athread_kernel_reduce_c12_host (m, its, kts, jts, ite, kte, jte,
        ims, ime, jms, jend, NX, NY, NZ, h_ap, &c12);
    //GPTLstop ("reduction");  
            
    c1[0] = c11;
    c1[1] = c12;
    MPI_Allreduce (c1, c2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double ac = c2[0] / c2[1];
    aps[m] = c2[1];
  
        //if(mytid[0]==0) printf("ac=%.3e\n",ac);
        //#pragma omp parallel for private(j,k,i) schedule(runtime) 
    // para_jte1(
    //     x0[index_x(i,k,j)]= x0[index_x(i,k,j)]+ac * h_p[index4b(i,k,j,m)];
    // )
    athread_kernel5_host (m, its, kts, jts, ite, kte, jte,
        ims, ime, jms, jend, NX, NY, NZ, ac, h_p, x0);

    // para_jend(
    //     h_r[index3(i,k,j)] = h_r[index3(i,k,j)]-ac * h_ap[index4(i,k,j,m)];
    // )
    athread_kernel6_host (m, its, kts, jts, ite, kte, jte,
        ims, ime, jms, jend, NX, NY, NZ, ac, h_r, h_ap);
    d = 0.0;
    //#pragma omp parallel for private(j,k,i) reduction(+:d) schedule(runtime) 
    // para_jend(
    //     d = d  + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
    // )
    athread_kernel_reduce_d_host (its, kts, jts, ite, kte, jte,
        ims, ime, jms, jend, NX, NY, NZ, h_r, &d);

    //#pragma omp parallel for private(j,k,i) schedule(runtime) 
    // para_jte(h_p[index4b(i,k,j,m+1)] = h_r[index3(i,k,j)];)
    athread_kernel3_host (m+1, its, kts, jts, ite, kte, jte, ims, ime, jms, 
      jend, NX, NY, NZ, h_p, h_r);

    int m_tmp = m + 1;
    //GPTLstart ("glob_updatehalo");  
    glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    //GPTLstop ("glob_updatehalo");  
        
    //GPTLstart ("precondition");  
    // svrasr(h_p,m_tmp, cm_helm, NX, NY, NZ);
    athread_svrasr_host (h_p, cm_helm, m_tmp, NX, NY, NZ);
    //GPTLstop ("precondition");  
    //GPTLstart ("glob_updatehalo");  
    glob_updatehalo_p_ (h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    //GPTLstop ("glob_updatehalo");  
    //GPTLstart ("maxtrixpro");  
    // matrixpro_c(a0, h_p+m_tmp*NG1,h_ar, its, ite, jts, jte, kts, kte, jend);
    athread_matrixpro_host (a0, h_p+m_tmp*NG1, h_ar, its, ite, jts, jte, kts, kte, jend);
    //GPTLstop ("maxtrixpro");  
    //GPTLstart ("reduction_cl");  
    for (l = 0; l <= m; l++) {
      double cl = 0.0;
      //#pragma omp parallel for private(j,k,i) reduction(+:cl) schedule(runtime) 
      // para_jend(cl = cl + h_ar[index3(i,k,j)]*h_ap[index4(i,k,j,l)];)
      athread_kernel_reduce_cl_host (l, its, kts, jts, ite, kte, jte,
          ims, ime, jms, jend, NX, NY, NZ, h_ar, h_ap, &cl);
      c1[l] = cl;
    }
    //GPTLstop ("reduction_cl");  

    c1[m+1] = d;
    MPI_Allreduce (c1, c2, m + 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    d = c2[m+1];

    // for (l = 0; l <= m; l++){
    //   b[l] = - c2[l] / aps[l];
    // }
    athread_kernel7_host (m, b, c2, aps);
    //#pragma omp parallel for private(j,k,i,l) schedule(runtime)   
    // for(j=jts-1;j<=jte+1;j++){
    //   for(l=0;l<=m;l++){
    //     for(k=kts-1;k<=kte+1;k++){
    //       for(i=its-1;i<=ite+1;i++){
    //         h_p[index4b(i,k,j,m+1)]+=b[l]*h_p[index4b(i,k,j,l)];//)
    //       }
    //     }
    //   }
    // } 
    for (l = 0; l <= m; l++) {
      athread_kernel8_host (l, m, its, kts, jts, ite, kte, jte,
          ims, ime, jms, jend, NX, NY, NZ, h_p, b);
    }
    //GPTLstop ("reductionp");  
    //GPTLstop ("gcr_iteration");  
    if (d <= ep || s == max_iteration-1 ) { 
      if(mytid[0]==0)printf("RES of gcr is %e in %d iterations\n", sqrt(d) ,s);
      break;
    }
    if (mytid[0]==0) printf("RES of gcr is %e in %d iterations\n", sqrt(d) ,s);
    if (s == 10) return;
  }
  // End iter
  // matrixpro_x(a0, x0,h_ap,its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
  athread_matrixpro_double_host (a0, x0, h_ap, its, ite, jts, jte, kts, kte, jend);

  // initialize(h_r,NG);
  memset (h_r, 0, sizeof(float) * NG);

  // para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]- h_ap[index3(i,k,j)];)
  athread_kernel9_host (its, kts, jts, ite, kte, jte,
      ims, ime, jms, jend, NX, NY, NZ, h_r, f0, h_ap);

  // c1[0] = norm2(h_r, NG);
  athread_kernel_reduce_norm2 (NG, h_r, c1);

  MPI_Allreduce(c1, c2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  c2[0] = sqrt (c2[0]);

  if(mytid[0] == 5 || mytid[0] == 12) {
    printf("RES of gmres is %e \n",c2[0]);
  }

  free (h_p);
  free (h_ap);
  free (h_ar);
  free (h_r);

  free (aps);
  free (b);
  free (c1);
  free (c2);

  return ;
}

#undef para_jte
#undef para_jte1
#undef para_jend
#undef printr