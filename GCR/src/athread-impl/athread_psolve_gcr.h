#ifndef SRC_ATHREAD_IMPL_ATHREAD_PSOLVE_GCR_H_
#define SRC_ATHREAD_IMPL_ATHREAD_PSOLVE_GCR_H_
extern void athread_kernel1_host (int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, 
                           int NX, int NY, int NZ, float* p, double* x0);
extern void athread_kernel2_host (int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, float* r, float* f0);
extern void athread_kernel3_host (int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, float* p, float* r);

extern void athread_kernel4_host (int iter_max, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int NX, int NY, int NZ, float* p);

extern void athread_kernel5_host (int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           double ac, float* p, double* x0);

extern void athread_kernel6_host (int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           double ac, float* r, float* ap);

extern void athread_kernel7_host (int m, double* b, double* c2, double* aps);

extern void athread_kernel8_host (int l, int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           float* p, double* b);

extern void athread_kernel9_host (int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           float* r, float* f0, float* ap);

extern void athread_kernel_reduce_d_host (int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float* r, double* d);

extern void athread_kernel_reduce_c11_host (int m, int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float* r, float *ap, double* c11);

extern void athread_kernel_reduce_c12_host (int m, int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float *ap, double* c12);

extern void athread_kernel_reduce_cl_host (int l, int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float *ar, float* ap, double *cl);

extern void athread_kernel_reduce_norm2 (int len, float* r, double *result);

extern void athread_matrixpro_host (float *a, float *b,float *c,
    int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern void athread_matrixpro_double_host (float *a, double *b,float *c,
    int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern void athread_svrasr_host (float* yl, float* b, int m, int ni, int nk, int nj);

#endif // SRC_ATHREAD_IMPL_ATHREAD_PSOLVE_GCR_H_