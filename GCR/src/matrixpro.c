#include <stdio.h>
#include "gcr.h"
extern void matrixpro_ca(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend){
    int i,k,j;
    int NX = ite - its + 3;
    int NY = kte - kts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
  #pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
    for(k=kts;k<=kte;k++){
      for(i=ibegin;i<=iend;i++){
       c[index3b(i,k,j)]  =                                  
                   +a[index_a(1 ,i,k,j)]*b[index3b(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index3b(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index3b(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index3b(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index3b(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index3b(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index3b(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index3b(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index3b(i-1,k  ,j+1)]      
                   +a[index_a(10,i,k,j)]*b[index3b(i  ,k-1,j  )]      
                   +a[index_a(11,i,k,j)]*b[index3b(i-1,k-1,j  )]      
                   +a[index_a(12,i,k,j)]*b[index3b(i+1,k-1,j  )]      
                   +a[index_a(13,i,k,j)]*b[index3b(i  ,k-1,j-1)]      
                   +a[index_a(14,i,k,j)]*b[index3b(i  ,k-1,j+1)]      
                   +a[index_a(15,i,k,j)]*b[index3b(i  ,k+1,j  )]      
                   +a[index_a(16,i,k,j)]*b[index3b(i-1,k+1,j  )]      
                   +a[index_a(17,i,k,j)]*b[index3b(i+1,k+1,j  )]      
                   +a[index_a(18,i,k,j)]*b[index3b(i  ,k+1,j-1)]      
                   +a[index_a(19,i,k,j)]*b[index3b(i  ,k+1,j+1)];
      }
    }
  }
    k=kts-1;
  #pragma omp parallel for private(j,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
      for(i=ibegin;i<=iend;i++){
       c[index3b(i,k,j)]  =                                 
                   +a[index_a(1 ,i,k,j)]*b[index3b(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index3b(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index3b(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index3b(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index3b(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index3b(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index3b(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index3b(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index3b(i-1,k  ,j+1)]      
                   +a[index_a(15,i,k,j)]*b[index3b(i  ,k+1,j  )]      
                   +a[index_a(16,i,k,j)]*b[index3b(i-1,k+1,j  )]      
                   +a[index_a(17,i,k,j)]*b[index3b(i+1,k+1,j  )]      
                   +a[index_a(18,i,k,j)]*b[index3b(i  ,k+1,j-1)]      
                   +a[index_a(19,i,k,j)]*b[index3b(i  ,k+1,j+1)];
      }
  }
    k=kte+1;
  #pragma omp parallel for private(j,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
      for(i=ibegin;i<=iend;i++){
       c[index3b(i,k,j)]  =                                 
                   +a[index_a(1 ,i,k,j)]*b[index3b(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index3b(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index3b(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index3b(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index3b(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index3b(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index3b(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index3b(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index3b(i-1,k  ,j+1)]      
                   +a[index_a(10,i,k,j)]*b[index3b(i  ,k-1,j  )]      
                   +a[index_a(11,i,k,j)]*b[index3b(i-1,k-1,j  )]      
                   +a[index_a(12,i,k,j)]*b[index3b(i+1,k-1,j  )]      
                   +a[index_a(13,i,k,j)]*b[index3b(i  ,k-1,j-1)]      
                   +a[index_a(14,i,k,j)]*b[index3b(i  ,k-1,j+1)];
      }
  }
}

extern void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend){
    int i,k,j;
    int NX = ite - its + 3;
    int NY = kte - kts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
  #pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
    for(k=kts;k<=kte;k++){
      for(i=ibegin;i<=iend;i++){
       c[index3(i,k,j)]  =                                  
                   +a[index_a(1 ,i,k,j)]*b[index3b(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index3b(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index3b(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index3b(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index3b(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index3b(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index3b(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index3b(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index3b(i-1,k  ,j+1)]      
                   +a[index_a(10,i,k,j)]*b[index3b(i  ,k-1,j  )]      
                   +a[index_a(11,i,k,j)]*b[index3b(i-1,k-1,j  )]      
                   +a[index_a(12,i,k,j)]*b[index3b(i+1,k-1,j  )]      
                   +a[index_a(13,i,k,j)]*b[index3b(i  ,k-1,j-1)]      
                   +a[index_a(14,i,k,j)]*b[index3b(i  ,k-1,j+1)]      
                   +a[index_a(15,i,k,j)]*b[index3b(i  ,k+1,j  )]      
                   +a[index_a(16,i,k,j)]*b[index3b(i-1,k+1,j  )]      
                   +a[index_a(17,i,k,j)]*b[index3b(i+1,k+1,j  )]      
                   +a[index_a(18,i,k,j)]*b[index3b(i  ,k+1,j-1)]      
                   +a[index_a(19,i,k,j)]*b[index3b(i  ,k+1,j+1)];
      }
    }
  }
    k=kts-1;
  #pragma omp parallel for private(j,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
      for(i=ibegin;i<=iend;i++){
       c[index3(i,k,j)]  =                                 
                   +a[index_a(1 ,i,k,j)]*b[index3b(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index3b(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index3b(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index3b(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index3b(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index3b(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index3b(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index3b(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index3b(i-1,k  ,j+1)]      
                   +a[index_a(15,i,k,j)]*b[index3b(i  ,k+1,j  )]      
                   +a[index_a(16,i,k,j)]*b[index3b(i-1,k+1,j  )]      
                   +a[index_a(17,i,k,j)]*b[index3b(i+1,k+1,j  )]      
                   +a[index_a(18,i,k,j)]*b[index3b(i  ,k+1,j-1)]      
                   +a[index_a(19,i,k,j)]*b[index3b(i  ,k+1,j+1)];
      }
  }
    k=kte+1;
  #pragma omp parallel for private(j,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
      for(i=ibegin;i<=iend;i++){
       c[index3(i,k,j)]  =                                 
                   +a[index_a(1 ,i,k,j)]*b[index3b(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index3b(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index3b(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index3b(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index3b(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index3b(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index3b(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index3b(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index3b(i-1,k  ,j+1)]      
                   +a[index_a(10,i,k,j)]*b[index3b(i  ,k-1,j  )]      
                   +a[index_a(11,i,k,j)]*b[index3b(i-1,k-1,j  )]      
                   +a[index_a(12,i,k,j)]*b[index3b(i+1,k-1,j  )]      
                   +a[index_a(13,i,k,j)]*b[index3b(i  ,k-1,j-1)]      
                   +a[index_a(14,i,k,j)]*b[index3b(i  ,k-1,j+1)];
      }
  }
}

extern void matrixpro_x(float *a,double *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend, int ims, int ime,  int jms, int jme){
    int i,k,j;
    int NX = ite - its + 3;
    int NY = kte - kts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
  #pragma omp parallel for private(j,k,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
    for(k=kts;k<=kte;k++){
      for(i=ibegin;i<=iend;i++){
       c[index3(i,k,j)]  =                                  
                   +a[index_a(1 ,i,k,j)]*b[index_x(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index_x(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index_x(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index_x(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index_x(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index_x(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index_x(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index_x(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index_x(i-1,k  ,j+1)]      
                   +a[index_a(10,i,k,j)]*b[index_x(i  ,k-1,j  )]      
                   +a[index_a(11,i,k,j)]*b[index_x(i-1,k-1,j  )]      
                   +a[index_a(12,i,k,j)]*b[index_x(i+1,k-1,j  )]      
                   +a[index_a(13,i,k,j)]*b[index_x(i  ,k-1,j-1)]      
                   +a[index_a(14,i,k,j)]*b[index_x(i  ,k-1,j+1)]      
                   +a[index_a(15,i,k,j)]*b[index_x(i  ,k+1,j  )]      
                   +a[index_a(16,i,k,j)]*b[index_x(i-1,k+1,j  )]      
                   +a[index_a(17,i,k,j)]*b[index_x(i+1,k+1,j  )]      
                   +a[index_a(18,i,k,j)]*b[index_x(i  ,k+1,j-1)]      
                   +a[index_a(19,i,k,j)]*b[index_x(i  ,k+1,j+1)];
      }
    }
  }
    k=kts-1;
  #pragma omp parallel for private(j,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
      for(i=ibegin;i<=iend;i++){
       c[index3(i,k,j)]  =                                 
                   +a[index_a(1 ,i,k,j)]*b[index_x(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index_x(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index_x(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index_x(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index_x(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index_x(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index_x(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index_x(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index_x(i-1,k  ,j+1)]      
                   +a[index_a(15,i,k,j)]*b[index_x(i  ,k+1,j  )]      
                   +a[index_a(16,i,k,j)]*b[index_x(i-1,k+1,j  )]      
                   +a[index_a(17,i,k,j)]*b[index_x(i+1,k+1,j  )]      
                   +a[index_a(18,i,k,j)]*b[index_x(i  ,k+1,j-1)]      
                   +a[index_a(19,i,k,j)]*b[index_x(i  ,k+1,j+1)];
      }
  }
    k=kte+1;
  #pragma omp parallel for private(j,i) schedule(runtime) collapse(2)
  for(j=jbegin;j<=jend;j++){
      for(i=ibegin;i<=iend;i++){
       c[index3(i,k,j)]  =                                 
                   +a[index_a(1 ,i,k,j)]*b[index_x(i  ,k  ,j  )]      
                   +a[index_a(2 ,i,k,j)]*b[index_x(i-1,k  ,j  )]      
                   +a[index_a(3 ,i,k,j)]*b[index_x(i+1,k  ,j  )]      
                   +a[index_a(4 ,i,k,j)]*b[index_x(i  ,k  ,j-1)]      
                   +a[index_a(5 ,i,k,j)]*b[index_x(i  ,k  ,j+1)]      
                   +a[index_a(6 ,i,k,j)]*b[index_x(i+1,k  ,j+1)]      
                   +a[index_a(7 ,i,k,j)]*b[index_x(i+1,k  ,j-1)]      
                   +a[index_a(8 ,i,k,j)]*b[index_x(i-1,k  ,j-1)]      
                   +a[index_a(9 ,i,k,j)]*b[index_x(i-1,k  ,j+1)]      
                   +a[index_a(10,i,k,j)]*b[index_x(i  ,k-1,j  )]      
                   +a[index_a(11,i,k,j)]*b[index_x(i-1,k-1,j  )]      
                   +a[index_a(12,i,k,j)]*b[index_x(i+1,k-1,j  )]      
                   +a[index_a(13,i,k,j)]*b[index_x(i  ,k-1,j-1)]      
                   +a[index_a(14,i,k,j)]*b[index_x(i  ,k-1,j+1)];
      }
  }
}

