#include "hip/hip_runtime.h"
#ifndef SVRASR
#define SVRASR
#define index_yl(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk +(m)*(ni*nk*nj))
#define index_b(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk+(m-1)*ni*nk*nj)
__device__ void warp_shfl(int warplen,int i, int k, int j, int ni, int nj, int nk,int rr, float *yl_bench, float * b){
     float bl1,bl2;
    if(threadIdx.x<=warplen && i>=1 &&i<=ni){
        bl2=b[index_b(i,k,j,2+rr)];
        bl1=b[index_b(i,k,j,1)];
    }
    for(int l=1;l<=warplen;l++){
        float yl_left=__shfl_up(*yl_bench, 1,64);
        //if(i>=1&&i<=ni){
        if(threadIdx.x==l){
            *yl_bench+=-yl_left*bl2;
            if(rr==1) *yl_bench *=bl1;
        }
    }
}

__device__ void warp_shfl_32(int offset,int i, int k, int j, int ni, int nj, int nk,int rr, float *yl_bench, float * b){
     float bl1,bl2;
    if(threadIdx.x>=1+offset && threadIdx.x<offset+32 && i>=1 &&i<=ni){
        bl2=b[index_b(i,k,j,2+rr)];
        bl1=b[index_b(i,k,j,1)];
    }
    for(int l=1+offset;l<32+offset;l++){
        float yl_left=__shfl_up(*yl_bench, 1,64);
        //if(i>=1&&i<=ni){
        if(threadIdx.x==l){
            *yl_bench+=-yl_left*bl2;
            if(rr==1) *yl_bench *=bl1;
        }
    }

}

__device__ void warp_shfl_64(int localid,int warplen, int ll, int i, int k, int j, int ni, int nj, int nk, int rr, float *yl_bench, float * b, float yl_orig){
     float bl1,bl2,bl4;
    if( i>=1 &&i<=ni){
        bl2=b[index_b(i,k,j,2+rr)];
        bl1=b[index_b(i,k,j,1)];
        bl4=b[index_b(i,k,j,4+rr)];
    }
    if(threadIdx.x==0 && ll!=0) *yl_bench=yl_orig;
    if(threadIdx.x==0 && ll==0) *yl_bench=yl_orig- *yl_bench*bl4;
    for(int l=1;l<=warplen;l++){
        float yl_left=__shfl_up(*yl_bench, 1,64);
        //if(i>=1&&i<=ni){
        if(localid==l){
            *yl_bench= yl_orig- yl_left*bl2 - *yl_bench*bl4;
            if(rr==1) *yl_bench *=bl1;
        }
    }
}

#define k_index()\
        kn=k+1;\
        if(rr==1){\
            k=nk+1-k;\
            kn=nk+1-kn;\
        }

#define save_yl()\
    if(i>=1&&i<=ni){\
        if(threadIdx.x>0 || ll==0)yl[index_yl(i,k,j,m)]=yl_bench;\
        if(kk<nk){\
            if(threadIdx.x>0 || ll==0){\
                yl_bench=yl[index_yl(i,kn,j,m)]-yl_bench*b[index_b(i,kn,j,4+rr)];\
            }else{\
                yl_bench=yl[index_yl(i,kn,j,m)];\
            }\
            if(rr==1&&i==ni) yl_bench *=b[index_b(i,kn,j,1)];\
        }\
    }

__global__ void cu_svrasr( float* __restrict__ yl,int m, float* __restrict__ b, int ni, int nk, int nj){
//#include "defindex.hpp"
    
    int ii=  threadIdx.x+1;
    int j= blockDim.y * blockIdx.y + threadIdx.y +1;
    //int warpid=(ii-1) / 32;
    int i=ii,k,offset=1,l,ll,ii1=ii,kn; 
    float yl_bench,yl_bench_new,yl_left,yl_right;
        yl_bench=0;
        //int localid=threadIdx.x%32;
        //int k_offset=threadIdx.x/32;
        int warplen=blockDim.x;
        int warplen1=warplen-1;
        int cnt=(ni-2) / warplen1 ;
        int ni_fix=(cnt+1)* (warplen1) +1 ;
        //if(threadIdx.x>=32) return;
    if(j>=1&&j<=nj){
    for(int rr=0;rr<=1;rr++){
        for(int ll=0;ll<=cnt;ll++){
            __syncthreads();
            i=ii1+warplen1*ll;
            k=1;
            yl_bench=0;
            if(rr==1){
                i=ni+1-i;k=nk+1-k;
            }
            if( i>=1&&i<=ni) {
                yl_bench=yl[index_yl(i,k,j,m)];
            }
            float bl1,bl2,bl1n;
            for(int kk=1;kk<=nk;kk++){ 
                k=kk;
                kn=kk+1;
                if(rr==1){
                    k=nk+1-k;
                    kn=nk+1-kn;
                }
                yl_bench_new=0;yl_left=0;
                int l_fix=0;
                if(ll==cnt) l_fix=(ni_fix-ni);
                if(i>=1&&i<=ni){
                    bl2=b[index_b(i,k,j,2+rr)];
                    bl1=b[index_b(i,k,j,1)];
                }
                for(l=1;l<=warplen-1-l_fix;l++){
                    yl_left=__shfl_up(yl_bench, offset,warplen);
                    //if(i>=1&&i<=ni){
                    if(threadIdx.x==l){
                        yl_bench+=-yl_left*bl2;
                        if(rr==1) yl_bench *=bl1;
                    }
                } 
                __syncthreads();
                if(i>=1&&i<=ni){
                    if(threadIdx.x>0 || ll==0)yl[index_yl(i,k,j,m)]=yl_bench;
                    if(kk<nk){
                        if(threadIdx.x>0 || ll==0){
                            yl_bench=yl[index_yl(i,kn,j,m)]-yl_bench*b[index_b(i,kn,j,4+rr)];
                        }else{
                            yl_bench=yl[index_yl(i,kn,j,m)];
                        }
                        if(rr==1&&i==ni) yl_bench *=b[index_b(i,kn,j,1)];
                    }
                }
            }//k
        }//l if(i)
        if(rr==0){
            if(j==nj&&i==1){
            yl[index_yl(ni,nk,nj,m)]=yl[index_yl(ni,nk,nj,m)]*b[index_b(ni,nk,nj,1)];
            }
        }
        //__syncthreads();
    }//r
    }//j
    }

__global__ void cu_svrasr_merge( float* __restrict__ yl,int m, float* __restrict__ b, int ni, int nk, int nj){
//#include "defindex.hpp"
    
    int ii= blockDim.x * blockIdx.x + threadIdx.x+1;
    int j= blockDim.y * blockIdx.y + threadIdx.y +1;
    //int warpid=(ii-1) / 32;
    int i=ii,k,offset=1,l,ll,ii1=ii,kn; 
    float yl_bench,yl_bench_new,yl_left,yl_right;
        yl_bench=0;
        int warplen=64;
        int warplen1=63;
        int cnt=(ni-2) / warplen1 ;
        int ni_fix=(cnt+1)* (warplen1) +1 ;
        int l_fix;

        int laneid=0;
        int localid=threadIdx.x;
        int k_offset=-laneid;
        //if(threadIdx.x==63) return;
    if(j>=1&&j<=nj){
    for(int rr=0;rr<=1;rr++){
        for(ll=0;ll<cnt;ll++){
            __syncthreads();
            i=ii1+warplen1*ll;
            if(rr==1) i=ni+1-i;
            int kk=1;
            k=kk;k_index();
            yl_bench=0;
            if( i>=1&&i<=ni) {
                yl_bench=yl[index_yl(i,k,j,m)];
            }
            warp_shfl(warplen1, i, k, j, ni,nj,nk, rr, &yl_bench, b);
                if( threadIdx.x<=63){
                save_yl(); 
                }
            for(int kk1=2;kk1<=nk;kk1++){ 
                kk=kk1+k_offset;
                //kk=kk1;
                k=kk;k_index();
                float yl_orig=yl[index_yl(i,k,j,m)];
                warp_shfl_64(localid,warplen1, ll,i, k, j, ni,nj,nk, rr, &yl_bench, b, yl_orig);
                if(i>=1&&i<=ni){
                    if(threadIdx.x>0 || ll==0)yl[index_yl(i,k,j,m)]=yl_bench;
                }
            }//k
            }//l if(i)
            ll=cnt;
            __syncthreads();
            i=ii1+warplen1*ll;
            k=1;
            if(rr==1){
                i=ni+1-i;k=nk+1-k;
            }
            yl_bench=0;
            if( i>=1&&i<=ni) {
                yl_bench=yl[index_yl(i,k,j,m)];
            }
            for(int kk=1;kk<=nk;kk++){ 
                k=kk;
                k_index();
                l_fix=(ni_fix-ni);
                warp_shfl(warplen1 -l_fix, i, k, j, ni,nj,nk, rr, &yl_bench, b);
                __syncthreads();
                save_yl();
            }//k
            if(rr==0){
            if(j==nj&&i==1){
            yl[index_yl(ni,nk,nj,m)]=yl[index_yl(ni,nk,nj,m)]*b[index_b(ni,nk,nj,1)];
            }
        }
        //__syncthreads();
    }//r
    }//j
  //  }
    }
__global__ void pacu_svrasr( float* __restrict__ yl,int m, float* __restrict__ b, int ni, int nk, int nj){
//#include "defindex.hpp"
    
    int ii= blockDim.x * blockIdx.x + threadIdx.x+1;
    int j= blockDim.y * blockIdx.y + threadIdx.y +1;
    //int warpid=(ii-1) / 32;
    int i=ii,k,offset=1,l,ll,ii1=ii,kn; 
    float yl_bench,yl_bench_new,yl_left,yl_right;
        yl_bench=0;
        int warplen=64;
        int warplen1=62;
        int cnt=(ni-2) / warplen1 ;
        int ni_fix=(cnt+1)* (warplen1) +1 ;
        int l_fix;

        int laneid=threadIdx.x/32;
        int localid=threadIdx.x%32+laneid;
        int k_offset=-laneid;
        if(threadIdx.x==63) return;
    if(j>=1&&j<=nj){
    for(int rr=0;rr<=1;rr++){
        for(ll=0;ll<cnt;ll++){
            __syncthreads();
            i=ii1+warplen1*ll;
            if(rr==1) i=ni+1-i;
            int kk=1;
            k=kk;k_index();
            yl_bench=0;
            if( i>=1&&i<=ni) {
                yl_bench=yl[index_yl(i,k,j,m)];
            }
#define debug
#ifdef debug
            warp_shfl(warplen1, i, k, j, ni,nj,nk, rr, &yl_bench, b);
                if( threadIdx.x<=31){
                save_yl(); 
                }
            k=kk+1+k_offset;k_index();
            warp_shfl_32(0,i, k, j, ni,nj,nk, rr, &yl_bench, b);
                if(  threadIdx.x<=62 ){
                    yl[index_yl(i,k,j,m)]=yl_bench;
                }
            for(int kk1=3;kk1<=nk;kk1++){ 
                kk=kk1+k_offset;
                //kk=kk1;
                k=kk;k_index();
                float yl_orig=yl[index_yl(i,k,j,m)];
                warp_shfl_64(localid,31, ll,i, k, j, ni,nj,nk, rr, &yl_bench, b, yl_orig);
                if(i>=1&&i<=ni){
                    if(threadIdx.x>0 || ll==0)yl[index_yl(i,k,j,m)]=yl_bench;
                }
            }//k
                //kk=nk+k_offset;
                k=nk;
                k_index();
                if( threadIdx.x>=32 && threadIdx.x<=62){
                yl_bench=yl[index_yl(i,k,j,m)]-yl_bench*b[index_b(i,k,j,4+rr)];\
                }
                warp_shfl_32(31,i, k, j, ni,nj,nk, rr, &yl_bench, b);
                    //if(i>=1&&i<=ni) 
                if( threadIdx.x>=32 && threadIdx.x<=62){
                    yl[index_yl(i,k,j,m)]=yl_bench;
                }
#endif
        }//l if(i)
            ll=cnt;
            __syncthreads();
            i=ii1+warplen1*ll;
            k=1;
            if(rr==1){
                i=ni+1-i;k=nk+1-k;
            }
            yl_bench=0;
            if( i>=1&&i<=ni) {
                yl_bench=yl[index_yl(i,k,j,m)];
            }
            for(int kk=1;kk<=nk;kk++){ 
                k=kk;
                k_index();
                l_fix=(ni_fix-ni);
                warp_shfl(warplen1 -l_fix, i, k, j, ni,nj,nk, rr, &yl_bench, b);
                __syncthreads();
                save_yl();
            }//k
            if(rr==0){
            if(j==nj&&i==1){
            yl[index_yl(ni,nk,nj,m)]=yl[index_yl(ni,nk,nj,m)]*b[index_b(ni,nk,nj,1)];
            }
        }
        //__syncthreads();
    }//r
    }//j
  //  }
    }
#endif
