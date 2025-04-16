#include <omp.h>
#define index_yl(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk +m*(ni*nk*nj))
#define index_b(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk+(m-1)*ni*nk*nj)

void svrasr( float*  __restrict__ yl,int m, float* __restrict__ b, int ni, int nk, int nj){
      int ng=ni*nk*nj;
      int i,k,j;   
  #pragma omp parallel for private(j,k,i) schedule(runtime) 
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

           yl[index_yl(ni,nk,nj,m)]=yl[index_yl(ni,nk,nj,m)]*b[index_b(ni,nk,nj,1)];
  #pragma omp parallel for private(j,k,i) schedule(runtime) 
         for(j=nj;j>=1;j-=1){
           for(i=ni-1;i>=1;i-=1){
             yl[index_yl(i,nk,j,m)]=(yl[index_yl(i,nk,j,m)]-b[index_b(i,nk,j,3)]*yl[index_yl(i+1,nk,j,m)])*b[index_b(i,nk,j,1)];
           }
           for(k=nk-1;k>=1;k-=1){
             for(i=1;i<=ni;i++){
               yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,5)]*yl[index_yl(i,k+1,j,m)];
             }
             yl[index_yl(ni,k,j,m)]=yl[index_yl(ni,k,j,m)]*b[index_b(ni,k,j,1)];
             for(i=ni-1;i>=1;i-=1){
                yl[index_yl(i,k,j,m)]=(yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,3)]*yl[index_yl(i+1,k,j,m)])*b[index_b(i,k,j,1)];
             }
           }
         }
        return;
      }

void svrasr_opt_c( float*  __restrict__ yl,int m, float* __restrict__  b, int ni, int nk, int nj){
      int ng=ni*nk*nj;
      int i,k,j;  
      float yl_bench,left;
//#pragma omp parallel for schedule(dynamic)
//        #pragma omp parallel for num_threads(2)
        for(j=1;j<=nj;j++){
              left=yl[index_yl(1,1,j,m)];
          for(i=2;i<=ni;i++){
              yl_bench=yl[index_yl(i,1,j,m)];
              yl_bench+= -b[index_b(i,1,j,2)]*left;
              left=yl_bench;
              yl[index_yl(i,1,j,m)]=left;
          }
          for(k=2;k<=nk;k++){
              
              //yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,2)]*yl[index_yl(i-1,k,j,m)];
//            #pragma omp simd
            for(i=1;i<=ni;i++){
              yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,4)]*yl[index_yl(i,k-1,j,m)];
            }

              left=yl[index_yl(1,k,j,m)];
            for(i=2;i<=ni;i++){
              yl_bench=yl[index_yl(i,k,j,m)];
              yl_bench+= -b[index_b(i,k,j,2)]*left;
              left=yl_bench;
              yl[index_yl(i,k,j,m)]=left;
              //yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,2)]*yl[index_yl(i-1,k,j,m)];
            }
          }
        }

           yl[index_yl(ni,nk,nj,m)]=yl[index_yl(ni,nk,nj,m)]*b[index_b(ni,nk,nj,1)];
//        #pragma omp parallel for num_threads(2)
         for(j=nj;j>=1;j-=1){
              left=yl[index_yl(ni,nk,j,m)];
           for(i=ni-1;i>=1;i-=1){
              yl_bench=yl[index_yl(i,nk,j,m)];
              yl_bench=(yl_bench -b[index_b(i,nk,j,3)]*left)*b[index_b(i,nk,j,1)];
              left=yl_bench;
              yl[index_yl(i,nk,j,m)]=left;
             //yl[index_yl(i,nk,j,m)]=(yl[index_yl(i,nk,j,m)]-b[index_b(i,nk,j,3)]*yl[index_yl(i+1,nk,j,m)])*b[index_b(i,nk,j,1)];
           }
           for(k=nk-1;k>=1;k-=1){
//            #pragma omp simd
             for(i=1;i<=ni;i++){
               yl[index_yl(i,k,j,m)]=yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,5)]*yl[index_yl(i,k+1,j,m)];
             }
              yl[index_yl(ni,k,j,m)]=yl[index_yl(ni,k,j,m)]*b[index_b(ni,k,j,1)];

              left=yl[index_yl(ni,k,j,m)];
             for(i=ni-1;i>=1;i-=1){
              yl_bench=yl[index_yl(i,k,j,m)];
              yl_bench=(yl_bench -b[index_b(i,k,j,3)]*left)*b[index_b(i,k,j,1)];
              left=yl_bench;
              yl[index_yl(i,k,j,m)]=left;
               // yl[index_yl(i,k,j,m)]=(yl[index_yl(i,k,j,m)]-b[index_b(i,k,j,3)]*yl[index_yl(i+1,k,j,m)])*b[index_b(i,k,j,1)];
             }
           }
         }
        return;
      }