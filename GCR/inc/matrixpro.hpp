#include "hip/hip_runtime.h"
#ifdef CUDA

#define smeIndex(i,k,j) (i + (k)*(blockDim.x)+ (j)*((blockDim.x)*(blockDim.y+2)))
__global__ void matrixpro_sme(float *d_a,float *d_var,float *c, int its, int ite, int jts, int jte, int kts, int kte, int jend){
    //dim3 threadsPerBlock(8, 4, 4);
    //dim3 numBlocks(22,91,90);
    
    int NX = ite - its + 3;
    int NY = kte - kts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
    int idx= blockDim.x * blockIdx.x + threadIdx.x;
    int kdx= blockDim.y * blockIdx.y + threadIdx.y;
    int jdx= blockDim.z * blockIdx.z + threadIdx.z;

    int i,j,k,k1;

    int ii=idx-2*(idx/blockDim.x)+ibegin-1;
    int jj=jdx-2*(jdx/blockDim.z)+jbegin-1;
    //int kk=kdx-2*(kdx/blockDim.y)+kts-1;
    int kk=kdx+kts;

    __shared__ float smem[16*6*6];

    i=threadIdx.x;
    k=threadIdx.y+1;
    j=threadIdx.z;
    //global memory to shared memory
    if(jj>=jbegin-1 && jj<=jend+1){
    if(kk>=kts && kk<=kte){
	if(ii>=ibegin-1 && ii<=iend+1){
            smem[smeIndex(i,k,j)] = d_var[index3b(ii,kk,jj)];
        if(k==1){smem[smeIndex( i, k-1, j)]=d_var[index3b( ii, kk-1, jj)];}
        if((k==blockDim.y && kk<kte)|| kk==kte){ smem[smeIndex(i,k+1,j)]=d_var[index3b(ii,kk+1,jj)];}
	}
    }}
	__syncthreads();

    if(jj>=jbegin && jj<=jend && j>0 && j<blockDim.z-1){
        if(ii>=ibegin && ii<=iend && i>0 && i<blockDim.x-1){
        //if(kk>=kts && kk<=kte && k>0 && k< blockDim.y-1){
            if(kk>=kts && kk<=kte ){
                 c[index3(ii,kk,jj)] =d_a[index_a(1 ,ii,kk,jj)] * smem[smeIndex(i  ,k  ,j  )]
                    +d_a[index_a(2 ,ii,kk,jj)] * smem[smeIndex(i-1,k  ,j  )]
                    +d_a[index_a(3 ,ii,kk,jj)] * smem[smeIndex(i+1,k  ,j  )]
                    +d_a[index_a(4 ,ii,kk,jj)] * smem[smeIndex(i  ,k  ,j-1)]
                    +d_a[index_a(5 ,ii,kk,jj)] * smem[smeIndex(i  ,k  ,j+1)]
                    +d_a[index_a(6 ,ii,kk,jj)] * smem[smeIndex(i+1,k  ,j+1)]
                    +d_a[index_a(7 ,ii,kk,jj)] * smem[smeIndex(i+1,k  ,j-1)]
                    +d_a[index_a(8 ,ii,kk,jj)] * smem[smeIndex(i-1,k  ,j-1)]
                    +d_a[index_a(9 ,ii,kk,jj)] * smem[smeIndex(i-1,k  ,j+1)]
                    +d_a[index_a(10,ii,kk,jj)] * smem[smeIndex(i  ,k-1,j  )]
                    +d_a[index_a(11,ii,kk,jj)] * smem[smeIndex(i-1,k-1,j  )]
                    +d_a[index_a(12,ii,kk,jj)] * smem[smeIndex(i+1,k-1,j  )]
                    +d_a[index_a(13,ii,kk,jj)] * smem[smeIndex(i  ,k-1,j-1)]
                    +d_a[index_a(14,ii,kk,jj)] * smem[smeIndex(i  ,k-1,j+1)]
                    +d_a[index_a(15,ii,kk,jj)] * smem[smeIndex(i  ,k+1,j  )]
                    +d_a[index_a(16,ii,kk,jj)] * smem[smeIndex(i-1,k+1,j  )]
                    +d_a[index_a(17,ii,kk,jj)] * smem[smeIndex(i+1,k+1,j  )]
                    +d_a[index_a(18,ii,kk,jj)] * smem[smeIndex(i  ,k+1,j-1)]
                    +d_a[index_a(19,ii,kk,jj)] * smem[smeIndex(i  ,k+1,j+1)];
            }
        
            if(kk==kts){
	        k1=k-1;
                c[index3(ii,kk-1,jj)] =d_a[index_a(1 ,ii,kk-1,jj)] * smem[smeIndex(i  ,k1  ,j  )]
                    +d_a[index_a(2 ,ii,kk-1,jj)] * smem[smeIndex(i-1,k1  ,j  )]
                    +d_a[index_a(3 ,ii,kk-1,jj)] * smem[smeIndex(i+1,k1  ,j  )]
                    +d_a[index_a(4 ,ii,kk-1,jj)] * smem[smeIndex(i  ,k1  ,j-1)]
                    +d_a[index_a(5 ,ii,kk-1,jj)] * smem[smeIndex(i  ,k1  ,j+1)]
                    +d_a[index_a(6 ,ii,kk-1,jj)] * smem[smeIndex(i+1,k1  ,j+1)]
                    +d_a[index_a(7 ,ii,kk-1,jj)] * smem[smeIndex(i+1,k1  ,j-1)]
                    +d_a[index_a(8 ,ii,kk-1,jj)] * smem[smeIndex(i-1,k1  ,j-1)]
                    +d_a[index_a(9 ,ii,kk-1,jj)] * smem[smeIndex(i-1,k1  ,j+1)]
                    +d_a[index_a(15,ii,kk-1,jj)] * smem[smeIndex(i  ,k1+1,j  )]
                    +d_a[index_a(16,ii,kk-1,jj)] * smem[smeIndex(i-1,k1+1,j  )]
                    +d_a[index_a(17,ii,kk-1,jj)] * smem[smeIndex(i+1,k1+1,j  )]
                    +d_a[index_a(18,ii,kk-1,jj)] * smem[smeIndex(i  ,k1+1,j-1)]
                    +d_a[index_a(19,ii,kk-1,jj)] * smem[smeIndex(i  ,k1+1,j+1)];
	    }
            
            if(kk==kte){
                k1=k+1;
                c[index3(ii,kk+1,jj)] =d_a[index_a(1 , ii, kk+1, jj)] * smem[smeIndex( i, k1, j)]
                    +d_a[index_a(2 ,ii,kk+1,jj)] * smem[smeIndex(i-1,k1  ,j  )]
                    +d_a[index_a(3 ,ii,kk+1,jj)] * smem[smeIndex(i+1,k1  ,j  )]
                    +d_a[index_a(4 ,ii,kk+1,jj)] * smem[smeIndex(i  ,k1  ,j-1)]
                    +d_a[index_a(5 ,ii,kk+1,jj)] * smem[smeIndex(i  ,k1  ,j+1)]
                    +d_a[index_a(6 ,ii,kk+1,jj)] * smem[smeIndex(i+1,k1  ,j+1)]
                    +d_a[index_a(7 ,ii,kk+1,jj)] * smem[smeIndex(i+1,k1  ,j-1)]
                    +d_a[index_a(8 ,ii,kk+1,jj)] * smem[smeIndex(i-1,k1  ,j-1)]
                    +d_a[index_a(9 ,ii,kk+1,jj)] * smem[smeIndex(i-1,k1  ,j+1)]
                    +d_a[index_a(10,ii,kk+1,jj)] * smem[smeIndex(i  ,k1-1,j  )]
                    +d_a[index_a(11,ii,kk+1,jj)] * smem[smeIndex(i-1,k1-1,j  )]
                    +d_a[index_a(12,ii,kk+1,jj)] * smem[smeIndex(i+1,k1-1,j  )]
                    +d_a[index_a(13,ii,kk+1,jj)] * smem[smeIndex(i  ,k1-1,j-1)]
                    +d_a[index_a(14,ii,kk+1,jj)] * smem[smeIndex(i  ,k1-1,j+1)];
            }
        }
    }
     
}
/*
#define index3c(i,k,j) ( (i)-its+(k-kts)*(ite-its+1)+(j-jts)*(ite-its+1)*(kte-kts+1) + (jte-jts+1) * (ite-its+1))
#define index3s(i,k,j) ( (i)-its+(j-jts)*(ite-its+1))
#define index3e(i,k,j) ( (i)-its+(j-jts)*(ite-its+1) + (jte-jts+1)*(ite-its+1)*(kte-kts+2))
//
#define index_as(m,i,k,j) (m + (i-its) * 14 +   (j-jts)*14*(ite-its+1))
#define index_ae(m,i,k,j) (m + (i-its) * 14 +   (j-jts)*14*(ite-its+1) + 14* (jte-jts+1) * (ite-its+1) + 19 *  (jte-jts+1) * (ite-its+1) *  (kte-kts+1))
#define index_a1(m,i,k,j) (m + (i-its) * 19 + (k-kts)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)* (kte-kts+1) + 14 * (jte-jts+1) * (ite-its+1))



__global__ void construct_coo_cu(float *a, int *I,  int *ra,  int* rb, float *rc, int its, int ite, int jts, int jte, int kts, int kte, int jend){
    int i,k,j,m;
    int NX = ite - its + 3;
    int NY = kte - kts + 3;
    int NZ = jte - jts + 3;
    int NG=(NX-2)*NY*(NZ-2);

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
    int idx= blockDim.x * blockIdx.x + threadIdx.x;
    int kdx= blockDim.y * blockIdx.y + threadIdx.y;
    int jdx= blockDim.z * blockIdx.z + threadIdx.z;
    idx+=ibegin;
    kdx+=kts-1;
    jdx+=jbegin;
    int i,j,k;
    i=idx; k=kdx; j=jdx;
 

    if(jdx>=jbegin && jdx<=jend){
        if(kdx>=kts && kdx<=kte){
            if(idx>=ibegin && idx<=iend){
        for(m=0;m<19;m++){
            ra[index_a1(m,i,k,j)] = index3c(i,k,j);
            rb[index_a1(m,i,k,j)] = index3b(i + I[m*3+0], k + I[m*3+1], j + I[m*3+2]);
            rc[index_a1(m,i,k,j)] = a[index_a(m+1 ,i,k,j)];
        }
      }
    }
  }

       if(jdx>=jbegin && jdx<=jend){
        if(kdx==kts-1){
            if(idx>=ibegin && idx<=iend){
        for(m=0;m<9;m++){
            ra[index_as(m,i,k,j)] = index3s(i,k,j);
            rb[index_as(m,i,k,j)] = index3b(i + I[m*3+0], k + I[m*3+1], j + I[m*3+2]);
            rc[index_as(m,i,k,j)] = a[index_a(m+1 ,i,k,j)];
        }
        for(m=14;m<19;m++){
            ra[index_as(m-5,i,k,j)] = index3s(i,k,j);
            rb[index_as(m-5,i,k,j)] = index3b(i + I[m*3+0], k + I[m*3+1], j + I[m*3+2]);
            rc[index_as(m-5,i,k,j)] = a[index_a(m+1 ,i,k,j)];
        }
      }
  }
    if(jdx>=jbegin && jdx<=jend){
        if(kdx==kte+1){
            if(idx>=ibegin && idx<=iend){
        for(m=0;m<14;m++){
            ra[index_ae(m,i,k,j)] = index3e(i,k,j);
            rb[index_ae(m,i,k,j)] = index3b(i + I[m*3+0], k + I[m*3+1], j + I[m*3+2]);
            rc[index_ae(m,i,k,j)] = a[index_a(m+1 ,i,k,j)];
        }
      }
  }
}
*/
__global__ void cu_matrixpro_x(float *d_a,double *d_var,float *c, int its, int ite, int jts, int jte, int kts, int kte, int jend, int ims, int ime,  int jms, int jme){
    //dim3 threadsPerBlock(8, 4, 4);
    //dim3 numBlocks(22,91,90);
    
    int NX = ite - its + 3;
    int NY = kte - kts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
	//if(idx>=iend || jdx>=jend || kdx>=kte) return;
    //d_var[index(idx,jdx,kdx)]=h_var[k][j][i]=var(i,j,k);//
    //xDIM=8;yDIM=4;zDim=4;
    int idx= blockDim.x * blockIdx.x + threadIdx.x;
    int kdx= blockDim.y * blockIdx.y + threadIdx.y;
    int jdx= blockDim.z * blockIdx.z + threadIdx.z;
    idx+=ibegin;
    kdx+=kts-1;
    jdx+=jbegin;
    int i,j,k;
    i=idx; k=kdx; j=jdx;
    
    //  for(j=jbegin;j<=jend;j++){
    //    for(k=kts-1;k<kte;k++){
    //      for(i=ibegin;i<=iend;i++){
  
    if(jdx>=jbegin && jdx<=jend){
        if(kdx>=kts && kdx<=kte){
            if(idx>=ibegin && idx<=iend){
                 c[index3(idx,kdx,jdx)] =d_a[index_a(1 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j  )]
                    +d_a[index_a(2 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j  )]
                    +d_a[index_a(3 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j  )]
                    +d_a[index_a(4 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j-1)]
                    +d_a[index_a(5 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j+1)]
                    +d_a[index_a(6 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j+1)]
                    +d_a[index_a(7 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j-1)]
                    +d_a[index_a(8 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j-1)]
                    +d_a[index_a(9 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j+1)]
                    +d_a[index_a(10,idx,kdx,jdx)] * d_var[index_x(i  ,k-1,j  )]
                    +d_a[index_a(11,idx,kdx,jdx)] * d_var[index_x(i-1,k-1,j  )]
                    +d_a[index_a(12,idx,kdx,jdx)] * d_var[index_x(i+1,k-1,j  )]
                    +d_a[index_a(13,idx,kdx,jdx)] * d_var[index_x(i  ,k-1,j-1)]
                    +d_a[index_a(14,idx,kdx,jdx)] * d_var[index_x(i  ,k-1,j+1)]
                    +d_a[index_a(15,idx,kdx,jdx)] * d_var[index_x(i  ,k+1,j  )]
                    +d_a[index_a(16,idx,kdx,jdx)] * d_var[index_x(i-1,k+1,j  )]
                    +d_a[index_a(17,idx,kdx,jdx)] * d_var[index_x(i+1,k+1,j  )]
                    +d_a[index_a(18,idx,kdx,jdx)] * d_var[index_x(i  ,k+1,j-1)]
                    +d_a[index_a(19,idx,kdx,jdx)] * d_var[index_x(i  ,k+1,j+1)];
            }
        }
    }
    
    
       if(jdx>=jbegin && jdx<=jend){
        if(kdx==kts-1){
            if(idx>=ibegin && idx<=iend){
            c[index3(idx,kdx,jdx)] =d_a[index_a(1 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j  )]
                +d_a[index_a(2 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j  )]
                +d_a[index_a(3 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j  )]
                +d_a[index_a(4 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j-1)]
                +d_a[index_a(5 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j+1)]
                +d_a[index_a(6 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j+1)]
                +d_a[index_a(7 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j-1)]
                +d_a[index_a(8 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j-1)]
                +d_a[index_a(9 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j+1)]
                +d_a[index_a(15,idx,kdx,jdx)] * d_var[index_x(i  ,k+1,j  )]
                +d_a[index_a(16,idx,kdx,jdx)] * d_var[index_x(i-1,k+1,j  )]
                +d_a[index_a(17,idx,kdx,jdx)] * d_var[index_x(i+1,k+1,j  )]
                +d_a[index_a(18,idx,kdx,jdx)] * d_var[index_x(i  ,k+1,j-1)]
                +d_a[index_a(19,idx,kdx,jdx)] * d_var[index_x(i  ,k+1,j+1)];
		}
        }
    }

    if(jdx>=jbegin && jdx<=jend){
        if(kdx==kte+1){
            if(idx>=ibegin && idx<=iend){
            c[index3(idx,kdx,jdx)] =d_a[index_a(1 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j  )]
                +d_a[index_a(2 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j  )]
                +d_a[index_a(3 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j  )]
                +d_a[index_a(4 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j-1)]
                +d_a[index_a(5 ,idx,kdx,jdx)] * d_var[index_x(i  ,k  ,j+1)]
                +d_a[index_a(6 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j+1)]
                +d_a[index_a(7 ,idx,kdx,jdx)] * d_var[index_x(i+1,k  ,j-1)]
                +d_a[index_a(8 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j-1)]
                +d_a[index_a(9 ,idx,kdx,jdx)] * d_var[index_x(i-1,k  ,j+1)]
                +d_a[index_a(10,idx,kdx,jdx)] * d_var[index_x(i  ,k-1,j  )]
                +d_a[index_a(11,idx,kdx,jdx)] * d_var[index_x(i-1,k-1,j  )]
                +d_a[index_a(12,idx,kdx,jdx)] * d_var[index_x(i+1,k-1,j  )]
                +d_a[index_a(13,idx,kdx,jdx)] * d_var[index_x(i  ,k-1,j-1)]
                +d_a[index_a(14,idx,kdx,jdx)] * d_var[index_x(i  ,k-1,j+1)];
		}
        }
    }
}



__global__ void cu_matrixpro(float *d_a,float *d_var,float *c, int its, int ite, int jts, int jte, int kts, int kte, int jend){
    //dim3 threadsPerBlock(8, 4, 4);
    //dim3 numBlocks(22,91,90);
    
    int NX = ite - its + 3;
    int NY = kte - kts + 3;

    int ibegin = its;
    int iend   = ite;
    int jbegin = jts;
	//if(idx>=iend || jdx>=jend || kdx>=kte) return;
    //d_var[index(idx,jdx,kdx)]=h_var[k][j][i]=var(i,j,k);//
    //xDIM=8;yDIM=4;zDim=4;
    int idx= blockDim.x * blockIdx.x + threadIdx.x;
    int kdx= blockDim.y * blockIdx.y + threadIdx.y;
    int jdx= blockDim.z * blockIdx.z + threadIdx.z;
    idx+=ibegin;
    kdx+=kts-1;
    jdx+=jbegin;
    int i,j,k;
    i=idx; k=kdx; j=jdx;
    
    //  for(j=jbegin;j<=jend;j++){
    //    for(k=kts-1;k<kte;k++){
    //      for(i=ibegin;i<=iend;i++){
  
    if(jdx>=jbegin && jdx<=jend){
        if(kdx>=kts && kdx<=kte){
            if(idx>=ibegin && idx<=iend){
                 c[index3(idx,kdx,jdx)] =d_a[index_a(1 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j  )]
                    +d_a[index_a(2 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j  )]
                    +d_a[index_a(3 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j  )]
                    +d_a[index_a(4 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j-1)]
                    +d_a[index_a(5 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j+1)]
                    +d_a[index_a(6 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j+1)]
                    +d_a[index_a(7 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j-1)]
                    +d_a[index_a(8 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j-1)]
                    +d_a[index_a(9 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j+1)]
                    +d_a[index_a(10,idx,kdx,jdx)] * d_var[index3b(i  ,k-1,j  )]
                    +d_a[index_a(11,idx,kdx,jdx)] * d_var[index3b(i-1,k-1,j  )]
                    +d_a[index_a(12,idx,kdx,jdx)] * d_var[index3b(i+1,k-1,j  )]
                    +d_a[index_a(13,idx,kdx,jdx)] * d_var[index3b(i  ,k-1,j-1)]
                    +d_a[index_a(14,idx,kdx,jdx)] * d_var[index3b(i  ,k-1,j+1)]
                    +d_a[index_a(15,idx,kdx,jdx)] * d_var[index3b(i  ,k+1,j  )]
                    +d_a[index_a(16,idx,kdx,jdx)] * d_var[index3b(i-1,k+1,j  )]
                    +d_a[index_a(17,idx,kdx,jdx)] * d_var[index3b(i+1,k+1,j  )]
                    +d_a[index_a(18,idx,kdx,jdx)] * d_var[index3b(i  ,k+1,j-1)]
                    +d_a[index_a(19,idx,kdx,jdx)] * d_var[index3b(i  ,k+1,j+1)];
            }
        }
    }
    
    
       if(jdx>=jbegin && jdx<=jend){
        if(kdx==kts-1){
            if(idx>=ibegin && idx<=iend){
            c[index3(idx,kdx,jdx)] =d_a[index_a(1 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j  )]
                +d_a[index_a(2 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j  )]
                +d_a[index_a(3 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j  )]
                +d_a[index_a(4 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j-1)]
                +d_a[index_a(5 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j+1)]
                +d_a[index_a(6 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j+1)]
                +d_a[index_a(7 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j-1)]
                +d_a[index_a(8 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j-1)]
                +d_a[index_a(9 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j+1)]
                +d_a[index_a(15,idx,kdx,jdx)] * d_var[index3b(i  ,k+1,j  )]
                +d_a[index_a(16,idx,kdx,jdx)] * d_var[index3b(i-1,k+1,j  )]
                +d_a[index_a(17,idx,kdx,jdx)] * d_var[index3b(i+1,k+1,j  )]
                +d_a[index_a(18,idx,kdx,jdx)] * d_var[index3b(i  ,k+1,j-1)]
                +d_a[index_a(19,idx,kdx,jdx)] * d_var[index3b(i  ,k+1,j+1)];
		}
        }
    }

    if(jdx>=jbegin && jdx<=jend){
        if(kdx==kte+1){
            if(idx>=ibegin && idx<=iend){
            c[index3(idx,kdx,jdx)] =d_a[index_a(1 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j  )]
                +d_a[index_a(2 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j  )]
                +d_a[index_a(3 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j  )]
                +d_a[index_a(4 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j-1)]
                +d_a[index_a(5 ,idx,kdx,jdx)] * d_var[index3b(i  ,k  ,j+1)]
                +d_a[index_a(6 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j+1)]
                +d_a[index_a(7 ,idx,kdx,jdx)] * d_var[index3b(i+1,k  ,j-1)]
                +d_a[index_a(8 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j-1)]
                +d_a[index_a(9 ,idx,kdx,jdx)] * d_var[index3b(i-1,k  ,j+1)]
                +d_a[index_a(10,idx,kdx,jdx)] * d_var[index3b(i  ,k-1,j  )]
                +d_a[index_a(11,idx,kdx,jdx)] * d_var[index3b(i-1,k-1,j  )]
                +d_a[index_a(12,idx,kdx,jdx)] * d_var[index3b(i+1,k-1,j  )]
                +d_a[index_a(13,idx,kdx,jdx)] * d_var[index3b(i  ,k-1,j-1)]
                +d_a[index_a(14,idx,kdx,jdx)] * d_var[index3b(i  ,k-1,j+1)];
		}
        }
    }
    /*
    idx= blockDim.x * blockIdx.x + threadIdx.x;
    kdx= blockDim.y * blockIdx.y + threadIdx.y;
    jdx= blockDim.z * blockIdx.z + threadIdx.z;
    int ii=idx-2*(idx/10)+ibegin-1;
    int jj=jdx-2*(jdx/6)+jbegin-1;
    int kk=kdx-2*(kdx/6)+kts-1;
 
    __shared__ float smem[360];
    //if(blockIdx.x==0&&blockIdx.y==0&&blockIdx.z==0){

    i=threadIdx.x;
    k=threadIdx.y;
    j=threadIdx.z;
    //global memory to shared memory
    if(jj>=jbegin-1 && jj<=jend+1){ 
    if(kk>=kts-1 && kk<=kte+1){ 
	if(ii>=ibegin-1 && ii<=iend+1){
            smem[smeIndex(i,k,j)]=d_var[index3b(ii,kk,jj)];
	}
    }}
	__syncthreads();

    int jend1=jbegin+5;
    int iend1=ibegin+9;
    int kte1=kte+5;
    jend1=jend;iend1=iend;kte1=kte;

    if(jj>=jbegin && jj<=jend1 && j>0 && j< 5){
        if(kk>=kts && kk<=kte1 && k>0 && k< 5){
            if(ii>=ibegin && ii<=iend1 && i>0 && i< 9){
                 c[index3(ii,kk,jj)] =d_a[index_a(1 ,ii,kk,jj)] * smem[smeIndex(i  ,k  ,j  )]
                    +d_a[index_a(2 ,ii,kk,jj)] * smem[smeIndex(i-1,k  ,j  )]
                    +d_a[index_a(3 ,ii,kk,jj)] * smem[smeIndex(i+1,k  ,j  )]
                    +d_a[index_a(4 ,ii,kk,jj)] * smem[smeIndex(i  ,k  ,j-1)]
                    +d_a[index_a(5 ,ii,kk,jj)] * smem[smeIndex(i  ,k  ,j+1)]
                    +d_a[index_a(6 ,ii,kk,jj)] * smem[smeIndex(i+1,k  ,j+1)]
                    +d_a[index_a(7 ,ii,kk,jj)] * smem[smeIndex(i+1,k  ,j-1)]
                    +d_a[index_a(8 ,ii,kk,jj)] * smem[smeIndex(i-1,k  ,j-1)]
                    +d_a[index_a(9 ,ii,kk,jj)] * smem[smeIndex(i-1,k  ,j+1)]
                    +d_a[index_a(10,ii,kk,jj)] * smem[smeIndex(i  ,k-1,j  )]
                    +d_a[index_a(11,ii,kk,jj)] * smem[smeIndex(i-1,k-1,j  )]
                    +d_a[index_a(12,ii,kk,jj)] * smem[smeIndex(i+1,k-1,j  )]
                    +d_a[index_a(13,ii,kk,jj)] * smem[smeIndex(i  ,k-1,j-1)]
                    +d_a[index_a(14,ii,kk,jj)] * smem[smeIndex(i  ,k-1,j+1)]
                    +d_a[index_a(15,ii,kk,jj)] * smem[smeIndex(i  ,k+1,j  )]
                    +d_a[index_a(16,ii,kk,jj)] * smem[smeIndex(i-1,k+1,j  )]
                    +d_a[index_a(17,ii,kk,jj)] * smem[smeIndex(i+1,k+1,j  )]
                    +d_a[index_a(18,ii,kk,jj)] * smem[smeIndex(i  ,k+1,j-1)]
                    +d_a[index_a(19,ii,kk,jj)] * smem[smeIndex(i  ,k+1,j+1)];
            }
        }
    }
    */
}
#endif
