#ifdef __sw_slave__
#include <crts.h>
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h"
#include "athread_param.h"

//#define ADV_DEBUG

//ATHREAD global definition 
__thread_local int blockIndex;

__inline__ double min1(double a, double b, double c, double d, double e, double f) { 
    double i = a; 
    if (i > b) i = b; 
    if (i > c) i = c; 
    if (i > d) i = d; 
    if (i > e) i = e; 
    if (i > f) i = f; 

    return i;
}

__inline__ double min2(double a, double b) { 
    double i = a; 
    if (i > b) i = b; 

    return i;
}

__inline__ double max1(double a, double b, double c, double d, double e, double f) { 
    double i = a; 
    if (i < b) i = b; 
    if (i < c) i = c; 
    if (i < d) i = d; 
    if (i < e) i = e; 
    if (i < f) i = f; 

    return i;
}

__inline__ double max2(double a, double b) { 
    double i = a; 
    if (i < b) i = b; 

    return i;
}


 void  advection_tracer_1_0(advection_Para *Paras){

     int blockId_x,blockId_y,blockId_z;
     int threadId_x,threadId_y,threadId_z;
#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif
  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  advection_Para* Para_s = (advection_Para*)ldm_malloc(sizeof(advection_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection_Para));

  //int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_dyu,sizeof(double) * imt * jmt, MEM_TO_LDM);
  int imt = Para_s->imt;
  int jmt = Para_s->jmt;
  int km = Para_s->km;

  //double d_wkd[imt],d_wkdjm1[imt];
  double d_v_sface[imt],d_wkb[imt],d_u_wface[imt];
  double dd_wkd[imt*2];

     //cuda grid and block affine translation
     int blockSum   = Para_s->cuGrid[0] * Para_s->cuGrid[1] * Para_s->cuGrid[2];
     int threadSum  = Para_s->cuBlock[0] * Para_s->cuBlock[1] * Para_s->cuBlock[2];
     int quotient   = blockSum / 64;
     int remainder  = blockSum % 64;
     int blockStart = quotient * blockIndex;
     if(remainder > 0)
         blockStart += (remainder > blockIndex ? blockIndex : remainder);
     int blockEnd = blockStart + quotient;
     if(remainder > 0)
         blockEnd   += (remainder > blockIndex ? 1 : 0);

     int i,j,k;

     for(int block = blockStart; block < blockEnd; block++) {
         blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
         blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
         blockId_x = block % Para_s->cuGrid[0];
         j = blockId_y * Para_s->cuBlock[1];
         k = blockId_z * Para_s->cuBlock[2];

         CRTS_dma_get(d_wkb, &Para_s->d_wkb[k * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(dd_wkd, &Para_s->d_wkd[k *jmt *imt+(j-1) *imt], sizeof(double) * Para_s->cuBlock[0] * 2);

         //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
         for (i = 0; i < threadSum; i++) {
             threadId_x = i % Para_s->cuBlock[0];
             //threadId_z = threadId / (Para_s->cuBlock[0] * Para_s->cuBlock[1]);
             //threadId_y = (threadId - threadId_z * (Para_s->cuBlock[0] * Para_s->cuBlock[1])) / Para_s->cuBlock[0];

            // int i = blockId_x * Para_s->cuBlock[0] + threadId_x;

             if(1<=i&&i<(imt-1)&& k< km)
             {
                 d_v_sface[i]=(d_wkb[i]+d_wkb[i+1]) *Para_s->d_hts[j *imt+i] *0.25;
             }
             if(1<=j&&j<jmt-1&& k< km)
             {
                 //d_u_wface[i]=(d_wkdjm1[i]+d_wkd[i]) *Para_s->d_htw[j *imt+i] *0.25;
                 d_u_wface[i]=(dd_wkd[i]+dd_wkd[imt+i]) *Para_s->d_htw[j *imt+i] *0.25;
             }

         }
         CRTS_dma_put(&Para_s->d_v_sface[k * imt * jmt + j * imt], d_v_sface, sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_put(&Para_s->d_u_wface[k * imt * jmt + j * imt], d_u_wface, sizeof(double) * Para_s->cuBlock[0]);
     }

#ifdef ADV_DEBUG
     penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_0 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(advection_Para));
}


 void  advection_tracer_2(advection2_Para *Paras){

     int blockId_x,blockId_y,blockId_z;
     int threadId_x,threadId_y,threadId_z;

#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif

  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  advection2_Para* Para_s = (advection2_Para*)ldm_malloc(sizeof(advection2_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

 // int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_tarea_r,sizeof(double) * imt * jmt, MEM_TO_LDM);
//  ret = CRTS_memcpy_sldm(array1_2d,Para_s->d_hue,sizeof(double) * imt * jmt, MEM_TO_LDM);
 // ret = CRTS_memcpy_sldm(array2_2d,Para_s->d_hts,sizeof(double) * imt * jmt, MEM_TO_LDM);
  int imt = Para_s->imt;
  int jmt = Para_s->jmt;
  int km = Para_s->km;

  double d_dts = Para_s->d_dts;
  double *d_odzt = Para_s->d_odzt;

  double adv_xy1,adv_xy2,adv_xy3,adv_xy4;   
   
  double wt1,wt2;   
  
  double adv_zz,adv_za,adv_zb1,adv_zb2;    
  
  double adv_x0,adv_y0,adv_xx,adv_yy;    
  
  double adv_c1,adv_c2,adv_zc; 
  
  //double d_at[imt],d_atjm1[imt],d_atjp1[imt];
  //double d_v_sface[imt],d_v_sfacejm1[imt];
  //double d_vit[imt],d_vitjm1[imt],d_vitjp1[imt];
  double d_u_wface[imt],d_ws[imt],d_wskp1[imt],d_atmax[imt],d_atmin[imt],d_at00[imt],d_atkm1[imt],d_atkp1[imt];

  double dd_at[imt*3];
  double dd_v_sface[2*imt];
  double dd_vit[imt*3];

  wt1=-1.0e10;    
  wt2=+1.0e10;

     //cuda grid and block affine translation
     int blockSum   = Para_s->cuGrid[0] * Para_s->cuGrid[1] * Para_s->cuGrid[2];
     int threadSum  = Para_s->cuBlock[0] * Para_s->cuBlock[1] * Para_s->cuBlock[2];
     int quotient   = blockSum / 64;
     int remainder  = blockSum % 64;
     int blockStart = quotient * blockIndex;
     if(remainder > 0)
         blockStart += (remainder > blockIndex ? blockIndex : remainder);
     int blockEnd = blockStart + quotient;
     if(remainder > 0)
         blockEnd   += (remainder > blockIndex ? 1 : 0);

     int i,j,k;

     for(int block = blockStart; block < blockEnd; block++) {
         blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
         blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
         blockId_x = block % Para_s->cuGrid[0];
         j = blockId_y * Para_s->cuBlock[1];
         k = blockId_z * Para_s->cuBlock[2];

         CRTS_dma_get(d_u_wface, &Para_s->d_u_wface[k * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_ws, &Para_s->d_ws[k * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_wskp1, &Para_s->d_ws[(k+1) * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atkm1, &Para_s->d_at[(k-1) * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atkp1, &Para_s->d_at[(k+1) * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);

         CRTS_dma_get(dd_at, &Para_s->d_at[k * imt * jmt + (j-1) * imt], sizeof(double) * Para_s->cuBlock[0] * 3);
         CRTS_dma_get(dd_v_sface, &Para_s->d_v_sface[k * imt * jmt + (j-1) * imt], sizeof(double) * Para_s->cuBlock[0] * 2);
         CRTS_dma_get(dd_vit, &Para_s->d_vit[k * imt * jmt + (j-1) * imt], sizeof(double) * Para_s->cuBlock[0] * 3);

         for (i = 0; i < threadSum; i++) {
             if(1<=i&&i<(imt-1)&&1<=j&&j<(jmt-1)&& k< km){

                 adv_x0=(dd_at[imt+i+1]+dd_at[imt+i]) *d_u_wface[i+1] *Para_s->d_tarea_r[j *imt+i]-
                         (dd_at[imt+i]+dd_at[imt+i-1]) *d_u_wface[i] *Para_s->d_tarea_r[j *imt+i];
                 adv_y0=((dd_at[2*imt+i]+dd_at[imt+i]) *dd_v_sface[imt+i]-(dd_at[imt+i]+dd_at[i]) *dd_v_sface[i]) *
                         Para_s->d_tarea_r[j *imt+i];
                 adv_xy1=-d_dts *(dd_at[imt+i+1]-dd_at[imt+i]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         d_u_wface[i+1] *d_u_wface[i+1] /(Para_s->d_htw[j *imt+i+1] *
                         Para_s->d_hun[j *imt+i+1]);
                 adv_xy2=d_dts *(dd_at[imt+i]-dd_at[imt+i-1]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         d_u_wface[i] *d_u_wface[i] /(Para_s->d_htw[j *imt+i] *
                         Para_s->d_hun[j *imt+i]);
                 adv_xy3=-d_dts *(dd_at[2*imt+i]-dd_at[imt+i]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         dd_v_sface[imt+i] *dd_v_sface[imt+i] /(Para_s->d_hts[j *imt+i] *
                         Para_s->d_hue[j *imt+i]);
                 adv_xy4=d_dts *(dd_at[imt+i]-dd_at[i]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         dd_v_sface[i] *dd_v_sface[i] /(Para_s->d_hts[(j-1) *imt+i] *
                         Para_s->d_hue[(j-1) *imt+i]);
                 adv_c1=-dd_at[imt+i] *(d_u_wface[i+1]-d_u_wface[i]) *Para_s->d_tarea_r[j *imt+i] *2.0;
                 adv_c2=-dd_at[imt+i] *(dd_v_sface[imt+i]-dd_v_sface[i]) *Para_s->d_tarea_r[j *imt+i] *2.0;

                 if(k>0 && k < (km-1))
                 {
                     adv_za=0.5 *Para_s->d_odzp[k] *d_ws[i] *(dd_at[imt+i]+d_atkm1[i])-
                             0.5 *Para_s->d_odzp[k] *d_wskp1[i] *(dd_at[imt+i]+d_atkp1[i]);
                     adv_zb1=-0.5 *Para_s->d_odzp[k] *d_ws[i] *d_ws[i] *d_odzt[k] *
                             (d_atkm1[i]-dd_at[imt+i]) *d_dts;
                     adv_zb2=0.5 *Para_s->d_odzp[k] *d_wskp1[i] *d_wskp1[i] *d_odzt[k+1] *
                             (dd_at[imt+i]-d_atkp1[i]) *d_dts;
                     adv_zc=-Para_s->d_odzp[k] *dd_at[imt+i] *(d_ws[i]-d_wskp1[i]);

                     d_atmax[i]=max2(max1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt1,
                             dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt1,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt1,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt1,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt1,
                             d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt1),
                              d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                              (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt1);

                     d_atmin[i]=min2(min1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt2,
                             dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt2,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt2,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt2,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt2,
                             d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt2),
                             d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt2);

                 }
                 else if(k==0)
                 {
                     adv_za=-0.5 *Para_s->d_odzp[0] *d_wskp1[i] *(d_atkp1[i]+dd_at[imt+i]);
                     adv_zb1=0.0;
                     adv_zb2=0.5 *Para_s->d_odzp[0] *d_wskp1[i] *d_wskp1[i] *d_odzt[1] *
                             (dd_at[imt+i]-d_atkp1[i]) *d_dts;
                     adv_zc=Para_s->d_odzp[0] *dd_at[imt+i] *d_wskp1[i];

                     d_atmax[i]=max1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt1,
                             dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt1,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt1,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt1,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt1,
                             d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt1);

                     d_atmin[i]=min1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt2,
                             dd_at[i] *Para_s->d_vit[(j-1) *imt+i]+(1.0e0-dd_vit[i]) *wt2,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt2,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt2,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt2,
                             d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt2);

                 }
                 else
                 {
                     adv_za=0.5 *Para_s->d_odzp[km-1] *d_ws[i] *(dd_at[imt+i]+d_atkm1[i]);
                     adv_zb1=-0.5 *Para_s->d_odzp[km-1] *d_ws[i] *d_ws[i] *d_odzt[km-1] *
                             (d_atkm1[i]-dd_at[imt+i]) *d_dts;
                     adv_zb2=0.0;
                     adv_zc=-Para_s->d_odzp[km-1] *dd_at[imt+i] *d_ws[i];

                     d_atmax[i]=max1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt1,
                            dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt1,
                            dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt1,
                            dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt1,
                            dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt1,
                            d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                            (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt1);

                     d_atmin[i]=min1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt2,
                            dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt2,
                            dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt2,
                            dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt2,
                            dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt2,
                            d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                            (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt2);

                 }

                 adv_xx=-(adv_x0+adv_xy1+adv_xy2+adv_c1);
                 adv_yy=-(adv_y0+adv_xy3+adv_xy4+adv_c2);
                 adv_zz=-(adv_zb1+adv_zb2+adv_za+adv_zc);
                 d_at00[i]=dd_at[imt+i]+(adv_xx+adv_yy+adv_zz) *d_dts;

             }

         }

         CRTS_dma_put(&Para_s->d_atmax[k * imt * jmt + j * imt], d_atmax, sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_put(&Para_s->d_atmin[k * imt * jmt + j * imt], d_atmin, sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_put(&Para_s->d_at00[k * imt * jmt + j * imt], d_at00, sizeof(double) * Para_s->cuBlock[0]);

     }

#ifdef ADV_DEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_2 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(advection2_Para));
}

 void  advection_tracer_3(advection2_Para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int threadId_x,threadId_y,threadId_z;
  int blockId_x,blockId_y,blockId_z;
#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif
  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  advection2_Para* Para_s = (advection2_Para*)ldm_malloc(sizeof(advection2_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

 // int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_tarea_r,sizeof(double) * imt * jmt, MEM_TO_LDM);
  //ret = CRTS_memcpy_sldm(array1_2d,Para_s->d_hue,sizeof(double) * imt * jmt, MEM_TO_LDM);
  //ret = CRTS_memcpy_sldm(array2_2d,Para_s->d_hts,sizeof(double) * imt * jmt, MEM_TO_LDM);
  int imt = Para_s->imt;
  int jmt = Para_s->jmt;
  int km = Para_s->km;

  double dd_at00[3*imt];

  double d_at00pk1[imt];
  double d_at00mk1[imt];

  double dd_atmax[3*imt];
  double d_atmaxpk1[imt];
  double d_atmaxmk1[imt];

  double dd_atmin[3*imt];
  double d_atminpk1[imt];
  double d_atminmk1[imt];

  double dd_at[3*imt];
  double d_atpk1[imt];
  double d_atmk1[imt];
  double d_u_wface[imt];

  double dd_v_wface[2*imt];
  double d_ws[imt];
  double d_wspk1[imt];
  double d_tarea_r[imt];

  double dd_hue[2*imt];

  double dd_hts[2*imt];
  double d_dts = Para_s->d_dts;
  double* d_odzt = Para_s->d_odzt;
  //double* d_adv_tt = Para_s->d_adv_tt;
  double d_adv_tt[imt];
  double adv_xy1,adv_xy2,adv_xy3,adv_xy4;
  double adv_zz,adv_za,adv_zb1,adv_zb2;
  double adv_x0,adv_y0,adv_xx,adv_yy;
  double adv_c1,adv_c2,adv_zc;

     //cuda grid and block affine translation
     int blockSum   = Para_s->cuGrid[0] * Para_s->cuGrid[1] * Para_s->cuGrid[2];
     int threadSum  = Para_s->cuBlock[0] * Para_s->cuBlock[1] * Para_s->cuBlock[2];
     int quotient   = blockSum / 64;
     int remainder  = blockSum % 64;
     int blockStart = quotient * blockIndex;
     if(remainder > 0)
         blockStart += (remainder > blockIndex ? blockIndex : remainder);
     int blockEnd = blockStart + quotient;
     if(remainder > 0)
         blockEnd   += (remainder > blockIndex ? 1 : 0);

     int i,j,k;

     for(int block = blockStart; block < blockEnd; block++) {
         blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
         blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
         blockId_x = block % Para_s->cuGrid[0];
         j = blockId_y * Para_s->cuBlock[1];
         k = blockId_z * Para_s->cuBlock[2];

         CRTS_dma_get(d_at00pk1, &Para_s->d_at00[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_at00mk1, &Para_s->d_at00[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atmaxpk1, &Para_s->d_atmax[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atmaxmk1, &Para_s->d_atmax[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atminpk1, &Para_s->d_atmin[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atminmk1, &Para_s->d_atmin[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atpk1, &Para_s->d_at[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_atmk1, &Para_s->d_at[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_u_wface, &Para_s->d_u_wface[k * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_ws, &Para_s->d_ws[k * jmt * imt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_wspk1, &Para_s->d_ws[(k + 1) * jmt * imt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(d_tarea_r, &Para_s->d_tarea_r[j*imt], sizeof(double) * Para_s->cuBlock[0]);

         CRTS_dma_get(dd_at00, &Para_s->d_at00[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
         CRTS_dma_get(dd_atmax, &Para_s->d_atmax[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
         CRTS_dma_get(dd_atmin, &Para_s->d_atmin[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
         CRTS_dma_get(dd_at, &Para_s->d_at[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
         CRTS_dma_get(dd_v_wface, &Para_s->d_v_sface[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 2);
         CRTS_dma_get(dd_hue, &Para_s->d_hue[(j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 2);
         CRTS_dma_get(dd_hts, &Para_s->d_hts[(j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 2);

         for (i = 0; i < threadSum; i++) {
             d_adv_tt[i] = 0.0;
             if (1 <= i && i < (imt - 1) && 1 <= j && j < (jmt - 1) && k==0) {

                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i + 1] > dd_atmax[imt+i+1] || dd_at00[imt+i + 1] < dd_atmin[imt+i+1]) && i <= imt - 3)) {
                     adv_xy1 = -(dd_at[imt+i + 1] -dd_at[imt+i]) * fabs(d_u_wface[i + 1]) * d_tarea_r[i];
                 } else {
                     adv_xy1 = -d_dts * (dd_at[imt+i + 1] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i + 1] * d_u_wface[i + 1] /
                               (Para_s->d_htw[j * imt + i + 1] * Para_s->d_hun[j * imt + i + 1]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i - 1] > dd_atmax[imt+i - 1] || dd_at00[imt+i - 1] < dd_atmin[imt+i - 1]) &&i >= 2)) {
                     adv_xy2 = (dd_at[imt+i] - dd_at[imt+i - 1]) * fabs(d_u_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy2 = d_dts * (dd_at[imt+i] - dd_at[imt+i - 1]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i] * d_u_wface[i] /
                               (Para_s->d_htw[j * imt + i] * Para_s->d_hun[j * imt + i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[2*imt+i] > dd_atmax[2*imt+i] || dd_at00[2*imt+i] < dd_atmin[2*imt+i]) && j <= jmt - 3)) {
                     adv_xy3 = -(dd_at[2*imt+i] - dd_at[imt+i]) * fabs(dd_v_wface[imt+i]) * d_tarea_r[i];
                 } else {
                     adv_xy3 = -d_dts * (dd_at[2*imt+i] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[imt+i] * dd_v_wface[imt+i] /
                               (dd_hts[imt+i] * dd_hue[imt+i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[i] > dd_atmax[i] || dd_at00[i] < dd_atmin[i]) && j >= 2)) {
                     adv_xy4 = (dd_at[imt+i] - dd_at[i]) * fabs(dd_v_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy4 = d_dts * (dd_at[imt+i] - dd_at[i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[i] * dd_v_wface[i] /
                               (dd_hts[i] * dd_hue[i]);
                 }
                 adv_zb1 = 0.0;
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((d_at00pk1[i] > d_atmaxpk1[i] ||
                       d_at00pk1[i] < d_atminpk1[i]) && k <= km - 2)) {
                     adv_zb2 = 0.5 * fabs(d_wspk1[i]) * Para_s->d_odzp[k] *
                               (dd_at[imt+i] - d_atpk1[i]);
                 } else {
                     adv_zb2 = 0.5 * Para_s->d_odzp[0] * Para_s->d_ws[1 * jmt * imt + j * imt + i] *
                               Para_s->d_ws[1 * jmt * imt + j * imt + i] * d_odzt[1] *
                               (Para_s->d_at[j * imt + i] - Para_s->d_at[1 * jmt * imt + j * imt + i]) * d_dts;
                 }
                 adv_za = -0.5 * Para_s->d_odzp[0] * Para_s->d_ws[1 * jmt * imt + j * imt + i] *
                          (Para_s->d_at[1 * jmt * imt + j * imt + i] + Para_s->d_at[j * imt + i]);
                 adv_zc = Para_s->d_odzp[0] * Para_s->d_at[j * imt + i] * Para_s->d_ws[1 * jmt * imt + j * imt + i];
                 adv_c1 = -dd_at[imt+i] * (d_u_wface[i + 1] - d_u_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_c2 = -dd_at[imt+i] * (dd_v_wface[imt+i] - dd_v_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_x0 = (dd_at[imt+i + 1] +dd_at[imt+i]) * d_u_wface[i + 1] * d_tarea_r[i] -
                          (dd_at[imt+i] + dd_at[imt+i - 1]) * d_u_wface[i] * d_tarea_r[i];
                 adv_y0 = ((dd_at[2*imt+i] + dd_at[imt+i]) * dd_v_wface[imt+i] -
                           (dd_at[imt+i] + dd_at[i]) * dd_v_wface[i]) * d_tarea_r[i];
                 adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
                 adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
                 adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
                 //d_adv_tt[k * jmt * imt + j * imt + i] = adv_xx + adv_yy + adv_zz;
                 d_adv_tt[i] = adv_xx + adv_yy + adv_zz;


             }
             if (k >= 1 && k < km - 1) {
                 if (1 <= j && j < (jmt - 1) && 1 <= i && i < (imt - 1)) {
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[imt+i + 1] > dd_atmax[imt+i+1] || dd_at00[imt+i + 1] < dd_atmin[imt+i+1]) && i <= imt - 3) {
                         adv_xy1 = -(dd_at[imt+i + 1] -dd_at[imt+i]) *fabs(d_u_wface[i + 1]) * d_tarea_r[i];
                     } else {
                         adv_xy1 = -d_dts * (dd_at[imt+i + 1] - dd_at[imt+i]) * 2.0 *
                                   d_tarea_r[i] * d_u_wface[i + 1] * d_u_wface[i + 1] /
                                   (Para_s->d_htw[j * imt + i + 1] * Para_s->d_hun[j * imt + i + 1]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[imt+i - 1] > dd_atmax[imt+i - 1] || dd_at00[imt+i - 1] < dd_atmin[imt+i - 1]) && i >= 2) {
                         adv_xy2 = (dd_at[imt+i] - dd_at[imt+i - 1]) * fabs(d_u_wface[i]) * d_tarea_r[i];
                     } else {
                         adv_xy2 = d_dts * (dd_at[imt+i] - dd_at[imt+i - 1]) * 2.0 *
                                   d_tarea_r[i] * d_u_wface[i] * d_u_wface[i] /
                                   (Para_s->d_htw[j * imt + i] * Para_s->d_hun[j * imt + i]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[2*imt+i] > dd_atmax[2*imt+i] || dd_at00[2*imt+i] < dd_atmin[2*imt+i]) && j <= jmt - 3) {
                         adv_xy3 = -(dd_at[2*imt+i] - dd_at[imt+i]) *
                                   fabs(dd_v_wface[imt+i]) * d_tarea_r[i];
                     } else {
                         adv_xy3 = -d_dts * (dd_at[2*imt+i] - dd_at[imt+i]) * 2.0 *
                                   d_tarea_r[i] * dd_v_wface[imt+i] * dd_v_wface[imt+i] /
                                   (dd_hts[imt+i] * dd_hue[imt+i]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[i] > dd_atmax[i] || dd_at00[i] < dd_atmin[i]) && j >= 2) {
                         adv_xy4 = (dd_at[imt+i] - dd_at[i]) * fabs(dd_v_wface[i]) * d_tarea_r[i];
                     } else {
                         adv_xy4 = d_dts * (dd_at[imt+i] - dd_at[i]) * 2.0 *
                                   d_tarea_r[i] * dd_v_wface[i] * dd_v_wface[i] /
                                   (dd_hts[i] * dd_hue[i]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          d_at00mk1[i] > d_atmaxmk1[i] || d_at00mk1[i] < d_atminmk1[i]) && k >= 1) {
                         adv_zb1 = -0.5 * fabs(d_ws[i]) * Para_s->d_odzp[k] * (d_atmk1[i] - dd_at[imt+i]);
                     } else {
                         adv_zb1 = -0.5 * Para_s->d_odzp[k] * d_ws[i] *
                                   d_ws[i] * d_odzt[k] * (d_atmk1[i] - dd_at[imt+i]) * d_dts;
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          d_at00pk1[i] > d_atmaxpk1[i] || d_at00pk1[i] < d_atminpk1[i]) && k <= km - 2) {
                         adv_zb2 = 0.5 * fabs(d_wspk1[i]) * Para_s->d_odzp[k] * (dd_at[imt+i] - d_atpk1[i]);
                     } else {
                         adv_zb2 = 0.5 * Para_s->d_odzp[k] * d_wspk1[i] * d_wspk1[i] * d_odzt[k + 1] *
                                   (dd_at[imt+i] - d_atpk1[i]) * d_dts;
                     }
                     adv_c1 = -dd_at[imt+i] * (d_u_wface[i + 1] - d_u_wface[i]) * d_tarea_r[i] * 2.0;
                     adv_c2 = -dd_at[imt+i] * (dd_v_wface[imt+i] - dd_v_wface[i]) * d_tarea_r[i] * 2.0;
                     adv_za = 0.5 * Para_s->d_odzp[k] * d_ws[i] *
                              (dd_at[imt+i] + d_atmk1[i]) - 0.5 * Para_s->d_odzp[k] * d_wspk1[i] *
                                                       (dd_at[imt+i] + d_atpk1[i]);
                     adv_zc = -Para_s->d_odzp[k] *dd_at[imt+i] * (d_ws[i] - d_wspk1[i]);
                     adv_x0 = (dd_at[imt+i + 1] + dd_at[imt+i]) * d_u_wface[i + 1] * d_tarea_r[i] -
                              (dd_at[imt+i] + dd_at[imt+i - 1]) * d_u_wface[i] * d_tarea_r[i];
                     adv_y0 = ((dd_at[2*imt+i] + dd_at[imt+i]) * dd_v_wface[imt+i] -
                               (dd_at[imt+i] + dd_at[i]) * dd_v_wface[i]) * d_tarea_r[i];
                     adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
                     adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
                     adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
                     // d_adv_tt[k * jmt * imt + j * imt + i] = adv_xx + adv_yy + adv_zz;
                     d_adv_tt[i] = adv_xx + adv_yy + adv_zz;
                 }
             }
             if (k == km - 1) {
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i + 1] > dd_atmax[imt+i+1] || dd_at00[imt+i + 1] < dd_atmin[imt+i+1]) && i <= imt - 3)) {
                     adv_xy1 = -(dd_at[imt+i + 1] - dd_at[imt+i]) *
                               fabs(d_u_wface[i + 1]) * d_tarea_r[i];
                 } else {
                     adv_xy1 = -d_dts * (dd_at[imt+i + 1] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i + 1] *
                               d_u_wface[i + 1] / (Para_s->d_htw[j * imt + i + 1] * Para_s->d_hun[j * imt + i + 1]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i - 1] > dd_atmax[imt+i - 1] || dd_at00[imt+i - 1] < dd_atmin[imt+i - 1]) && i >= 2)) {
                     adv_xy2 = (dd_at[imt+i] - dd_at[imt+i - 1]) * fabs(d_u_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy2 = d_dts * (dd_at[imt+i] - dd_at[imt+i - 1]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i] *
                               d_u_wface[i] / (Para_s->d_htw[j * imt + i] * Para_s->d_hun[j * imt + i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[2*imt+i] > dd_atmax[2*imt+i] || dd_at00[2*imt+i] < dd_atmin[2*imt+i]) && j <= jmt - 3)) {
                     adv_xy3 = -(dd_at[2*imt+i] - dd_at[imt+i]) *
                               fabs(dd_v_wface[imt+i]) * d_tarea_r[i];
                 } else {
                     adv_xy3 = -d_dts * (dd_at[2*imt+i] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[imt+i] * dd_v_wface[imt+i] /
                               (dd_hts[imt+i] * dd_hue[imt+i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[i] > dd_atmax[i] || dd_at00[i] < dd_atmin[i]) && j >= 2)) {
                     adv_xy4 = (dd_at[imt+i] - dd_at[i]) *
                               fabs(dd_v_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy4 = d_dts * (dd_at[imt+i] - dd_at[i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[i] * dd_v_wface[i] /
                               (dd_hts[i] * dd_hue[i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     d_at00mk1[i] > d_atmaxmk1[i] || d_at00mk1[i] < d_atminmk1[i]) {
                     adv_zb1 = -0.5 * fabs(d_ws[i]) * Para_s->d_odzp[k] * (d_atmk1[i] - dd_at[imt+i]);
                 } else {
                     adv_zb1 = -0.5 * Para_s->d_odzp[k] * d_ws[i] * d_ws[i] * d_odzt[k] *
                               (d_atmk1[i] - dd_at[imt+i]) * d_dts;
                 }
                 adv_zb2 = 0.0;
                 adv_c1 = -dd_at[imt+i] * (d_u_wface[i + 1] - d_u_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_c2 = -dd_at[imt+i] * (dd_v_wface[imt+i] - dd_v_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_za = 0.5 * Para_s->d_odzp[km - 1] * Para_s->d_ws[(km - 1) * jmt * imt + j * imt + i] *
                          (Para_s->d_at[(km - 1) * jmt * imt + j * imt + i] +
                           Para_s->d_at[(km - 2) * jmt * imt + j * imt + i]);
                 adv_zc = -Para_s->d_odzp[km - 1] * Para_s->d_at[(km - 1) * jmt * imt + j * imt + i] *
                          Para_s->d_ws[(km - 1) * jmt * imt + j * imt + i];
                 adv_x0 = (dd_at[imt+i + 1] +dd_at[imt+i]) * d_u_wface[i + 1] * d_tarea_r[i] -
                          (dd_at[imt+i] + dd_at[imt+i - 1]) * d_u_wface[i] * d_tarea_r[i];
                 adv_y0 = ((dd_at[2*imt+i] + dd_at[imt+i]) * dd_v_wface[imt+i] - (dd_at[imt+i] + dd_at[i]) *
                                                                   dd_v_wface[i]) * d_tarea_r[i];
                 adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
                 adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
                 adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
                 // d_adv_tt[k * jmt * imt + j * imt + i] = adv_xx + adv_yy + adv_zz;
                     d_adv_tt[i] = adv_xx + adv_yy + adv_zz;
             }

         }

         CRTS_dma_put(&Para_s->d_adv_tt[k * imt * jmt + j * imt], d_adv_tt, sizeof(double) * Para_s->cuBlock[0]);
     }


#ifdef ADV_DEBUG 
     penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_3 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(advection2_Para));
}

#endif // __sw_slave__
