#include "gcr.h"
#include <omp.h>
#define index_yl(i,k,j,m) (i-its+(k-kts)*ni+(j-1)*ni*nk +m*(ni*nk*nj))
#define index_cm(i,k,j,m) (i-its+1+(k-kts+1)*NX+(j-jts+1)*NX*NY +(m-1)*(NX*NY*NZ))
#define index_b(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk+(m-1)*ni*nk*nj)
/*
   void ILU_3( float*  __restrict__ a_helm, float* __restrict__ b, int ni, int nk, int nj){
   int ng=ni*nk*nj;
   int i,k,j;
   for(j=1;j<=nj;j++){
   for(i=2;i<=ni;i++){
   yl[index_yl(i,1,j,m)]=yl[index_yl(i,1,j,m)]-b[index_b(i,1,j,2)]*yl[index_yl(i-1,1,j,m)];
   }
   for(k=2;k<=nk;k++){
   for(i=1;i<=ni;i++){
   yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,4)]*yl[index_yl(i,k-1,j,m)];
   }

   for(i=2;i<=ni;i++){
   yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,2)]*yl[index_yl(i-1,k,j,m)];
   }
   }
   }
   */

extern "C" void __module_halo_MOD_glob_updatehalo_real_3d(float* h_p, int *m,int *iter_max); 
//extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte){
//__module_halo_MOD_glob_updatehalo_real_3d(h_p,m,iter_max);
//     __module_halo_MOD_glob_updatehalo_real_3d(h_p,m,iter_max);
//        }

//extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);

extern "C" void ILU_5(float *a, float *T, int jend,  int  its, int ite, int jts, int jte, int kts, int kte){

    int i,k,j;
    int NX = ite - its + 3;
    int NY = kte - kts + 3;
    int NZ = jte - jts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;

    size_t size_T = 7 * (ite-its+3) * (kte-kts+3) * (jte-jts+3);
    for(i=0;i<NX*NY*NZ;i++){
        T[i]  = 1.0;
    }
//#pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
    for(j=jts;j<=jte;j++){
        for(k=kts-1;k<=kte+1;k++){
            for(i=its;i<=ite;i++){
                T[index_cm(i,k,j, 1)]  = a[index_a(1,i,k,j)];
                T[index_cm(i,k,j, 2)]  = a[index_a(2,i,k,j)];
                T[index_cm(i,k,j, 3)]  = a[index_a(3,i,k,j)];
                T[index_cm(i,k,j, 4)]  = a[index_a(10,i,k,j)];
                T[index_cm(i,k,j, 5)]  = a[index_a(15,i,k,j)];
                T[index_cm(i,k,j, 6)]  = a[index_a(4,i,k,j)];
                T[index_cm(i,k,j, 7)]  = a[index_a(5,i,k,j)];
            }
        }
    }
    
    int iter_max=7;
    for(int l=2;l<=iter_max;l++){
        __module_halo_MOD_glob_updatehalo_real_3d(T, &l,&iter_max); 
    }
    /*
    int l1=1,l2=6;
    __module_halo_MOD_glob_updatehalo_real_3d(T, &l1,&iter_max); 
    __module_halo_MOD_glob_updatehalo_real_3d(T, &l2,&iter_max); 
    */
    /*
  for(j=jts-1;j<=jte+1;j++){
      for(k=kts-1;k<=kte+1;k++){
          for(i=its-1;i<=ite+1;i++){ 
              T[index_cm(i,k,j,1)] = 1.0 / T[index_cm(i,k,j,1)];
          }
      }
  }*/
//#pragma omp parallel for private(j,k) schedule(runtime) collapse(2)
  for(j=jts-1;j<=jte+1;j++){
      if(j>=jts){
          for(i=its-1;i<=ite+1;i++){
              for(k=kts-1;k<=kte+1;k++){
                  T[index_cm(i,k,j,6)] = T[index_cm(i,k,j,6)] / T[index_cm(i,k,j-1,1)];
                  T[index_cm(i,k,j,1)] = T[index_cm(i,k,j,1)] - T[index_cm(i,k,j,6)] * T[index_cm(i,k,j-1,7)];
              }
          }
      }
      for(k=kts-1;k<=kte+1;k++){
          if(k>=kts){
              for(i=its-1;i<=ite+1;i++){
                  T[index_cm(i,k,j,4)] = T[index_cm(i,k,j,4)] / T[index_cm(i,k-1,j,1)];
                  T[index_cm(i,k,j,1)] = T[index_cm(i,k,j,1)] - T[index_cm(i,k,j,4)] * T[index_cm(i,k-1,j,5)];
              }
          }

          T[index_cm(its-1,k,j,1)] = 1.0 / T[index_cm(its-1,k,j,1)];
          for(i=its;i<=ite+1;i++){
              T[index_cm(i,k,j,2)] = T[index_cm(i,k,j,2)] * T[index_cm(i-1,k,j,1)];
              T[index_cm(i,k,j,1)] = 1.0 / (T[index_cm(i,k,j,1)] - T[index_cm(i,k,j,2)] * T[index_cm(i-1,k,j,3)]);
              //T[index_cm(i,k,j,1)] =  T[index_cm(i,k,j,1)]  - T[index_cm(i,k,j,2)] * T[index_cm(i-1,k,j,3)] ;
          }
      }
  }
}
    
    /*    for(j=jts-1;j<=jte+1;j++){
         //   for(k=kts-1;k<=kte+1;k++){
    //k=kts-1; //j=jbegin-1;
        if(j>=jts){
            for(k=kts-1;k<=kte+1;k++){
                for(i=its-1;i<=ite+1;i++){
                    T[index_cm(i,k,j,6)] = T[index_cm(i,k,j,6)] / T[index_cm(i,k,j-1,1)];
                    T[index_cm(i,k,j,1)] = T[index_cm(i,k,j,1)] - T[index_cm(i,k,j,6)] * T[index_cm(i,k,j-1,7)] ;
                }
            }
        }
        for(k=kts-1;k<=kte+1;k++){
            if(k>=kts){
                for(i=its-1;i<=ite+1;i++){
                    T[index_cm(i,k,j,4)] = T[index_cm(i,k,j,4)] / T[index_cm(i,k-1,j,1)];
                    T[index_cm(i,k,j,1)] = T[index_cm(i,k,j,1)] - T[index_cm(i,k,j,4)] * T[index_cm(i,k-1,j,5)] ;
                }
            }
            for(i=its;i<=ite+1;i++){
                T[index_cm(i,k,j,2)] = T[index_cm(i,k,j,2)] / T[index_cm(i-1,k,j,1)];
                T[index_cm(i,k,j,1)] =  T[index_cm(i,k,j,1)]  - T[index_cm(i,k,j,2)] * T[index_cm(i-1,k,j,3)] ;
            }
        }
    }
}
*/
