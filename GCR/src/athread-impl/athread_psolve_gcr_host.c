#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN (athread_kernel1_slave) (struct KernelParam1 *);
extern void SLAVE_FUN (athread_kernel2_slave) (struct KernelParam2 *);
extern void SLAVE_FUN (athread_kernel3_slave) (struct KernelParam3 *);
extern void SLAVE_FUN (athread_kernel4_slave) (struct KernelParam4 *);
extern void SLAVE_FUN (athread_kernel5_slave) (struct KernelParam5 *);
extern void SLAVE_FUN (athread_kernel6_slave) (struct KernelParam6 *);
extern void SLAVE_FUN (athread_kernel7_slave) (struct KernelParam7 *);
extern void SLAVE_FUN (athread_kernel8_slave) (struct KernelParam8 *);
extern void SLAVE_FUN (athread_kernel9_slave) (struct KernelParam9 *);

extern void SLAVE_FUN (athread_kernel_reduce_d_slave) (struct KernelReudceDParam *);
extern void SLAVE_FUN (athread_kernel_reduce_c11_slave) (struct KernelReudceC11Param *);
extern void SLAVE_FUN (athread_kernel_reduce_c12_slave) (struct KernelReudceC12Param *);
extern void SLAVE_FUN (athread_kernel_reduce_cl_slave) (struct KernelReudceClParam *);
extern void SLAVE_FUN (athread_kernel_reduce_norm2_slave) (struct KernelReudceNorm2Param *);

extern void SLAVE_FUN (athread_matrixpro_kts_slave) (struct MatrixProParam *);
extern void SLAVE_FUN (athread_matrixpro_kts_kte_slave) (struct MatrixProParam *);
extern void SLAVE_FUN (athread_matrixpro_kte_slave) (struct MatrixProParam *);

extern void SLAVE_FUN (athread_matrixpro_double_kts_slave) (struct MatrixProDoubleParam *);
extern void SLAVE_FUN (athread_matrixpro_double_kts_kte_slave) (struct MatrixProDoubleParam *);
extern void SLAVE_FUN (athread_matrixpro_double_kte_slave) (struct MatrixProDoubleParam *);

extern void SLAVE_FUN (athread_svrasr_slave) (struct SvrasrParam *);

void athread_kernel1_host (int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, 
                           int NX, int NY, int NZ, float* p, double* x0) {
  struct KernelParam1 para;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.p   = p;
  para.x0  = x0;
	athread_spawn (athread_kernel1_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel2_host (int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, float* r, float* f0) {
  struct KernelParam2 para;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.r   = r;
  para.f0  = f0;
	athread_spawn (athread_kernel2_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel3_host (int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, float* p, float* r) {
  struct KernelParam3 para;
  para.m   = m;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.p   = p;
  para.r   = r;
	athread_spawn (athread_kernel3_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel4_host (int iter_max, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int NX, int NY, int NZ, float* p) {
  struct KernelParam4 para;
  para.iter_max = iter_max;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.p   = p;
	athread_spawn (athread_kernel4_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel5_host (int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           double ac, float* p, double* x0) {
  struct KernelParam5 para;
  para.m   = m;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.ac  = ac;
  para.p   = p;
  para.x0  = x0;
	athread_spawn (athread_kernel5_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel6_host (int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           double ac, float* r, float* ap) {
  struct KernelParam6 para;
  para.m   = m;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.ac  = ac;
  para.r   = r;
  para.ap  = ap;
	athread_spawn (athread_kernel6_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel7_host (int m, double* b, double* c2, double* aps) {
  struct KernelParam7 para;
  para.m   = m;
  para.b   = b;
  para.c2  = c2;
  para.aps = aps;
	athread_spawn (athread_kernel7_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel8_host (int l, int m, int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           float* p, double* b) {
  struct KernelParam8 para;
  para.l   = l;
  para.m   = m;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.p   = p;
  para.b   = b;
	athread_spawn (athread_kernel8_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel9_host (int its, int kts, int jts, 
                           int ite, int kte, int jte,
                           int ims, int ime, int jms, int jend,
                           int NX, int NY, int NZ, 
                           float* r, float* f0, float* ap) {
  struct KernelParam9 para;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.r   = r;
  para.f0  = f0;
  para.ap  = ap;
	athread_spawn (athread_kernel9_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel_reduce_d_host (int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float* r, double* d) {
  struct KernelReudceDParam para;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.r   = r;
  para.d   = d;
	athread_spawn (athread_kernel_reduce_d_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel_reduce_c11_host (int m, int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float* r, float *ap, double* c11) {
  struct KernelReudceC11Param para;
  para.m   = m;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.r   = r;
  para.ap  = ap;
  para.c11 = c11;
	athread_spawn (athread_kernel_reduce_c11_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel_reduce_c12_host (int m, int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float *ap, double* c12) {
  struct KernelReudceC12Param para;
  para.m   = m;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.ap  = ap;
  para.c12 = c12;
	athread_spawn (athread_kernel_reduce_c12_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel_reduce_cl_host (int l, int its, int kts, int jts, 
                                   int ite, int kte, int jte,
                                   int ims, int ime, int jms, int jend,
                                   int NX, int NY, int NZ, float *ar, float* ap, double *cl) {
  struct KernelReudceClParam para;
  para.l   = l;
  para.its = its;
  para.kts = kts;
  para.jts = jts;
  para.ite = ite;
  para.kte = kte;
  para.jte = jte;
  para.ims = ims;
  para.ime = ime;
  para.jms = jms;
  para.jend = jend;
  para.NX  = NX;
  para.NY  = NY;
  para.NZ  = NZ;
  para.ar  = ar;
  para.ap  = ap;
  para.cl  = cl;
	athread_spawn (athread_kernel_reduce_cl_slave, &para);
	athread_join ();
  return ;
}

void athread_kernel_reduce_norm2 (int len, float* r, double *result) {
  struct KernelReudceNorm2Param para;
  para.len    = len;
  para.r      = r;
  para.result = result;
	athread_spawn (athread_kernel_reduce_norm2_slave, &para);
	athread_join ();
  return ;
}

void athread_matrixpro_host (float *a,float *b,float *c,
    int its, int ite, int jts, int jte, int kts, int kte, int jend) {

  struct MatrixProParam para;
  para.its  = its;
  para.ite  = ite;
  para.jts  = jts;
  para.jte  = jte;
  para.kts  = kts;
  para.kte  = kte;
  para.jend = jend;
  para.a = a;
  para.b = b;
  para.c = c;

	athread_spawn (athread_matrixpro_kts_slave, &para);
	athread_join ();

	athread_spawn (athread_matrixpro_kts_kte_slave, &para);
	athread_join ();

	athread_spawn (athread_matrixpro_kte_slave, &para);
	athread_join ();
  return ;
}

void athread_matrixpro_double_host (float *a, double *b, float *c,
    int its, int ite, int jts, int jte, int kts, int kte, int jend) {

  struct MatrixProDoubleParam para;
  para.its  = its;
  para.ite  = ite;
  para.jts  = jts;
  para.jte  = jte;
  para.kts  = kts;
  para.kte  = kte;
  para.jend = jend;
  para.a = a;
  para.b = b;
  para.c = c;

	athread_spawn (athread_matrixpro_double_kts_slave, &para);
	athread_join ();

	athread_spawn (athread_matrixpro_double_kts_kte_slave, &para);
	athread_join ();

	athread_spawn (athread_matrixpro_double_kte_slave, &para);
	athread_join ();
  return ;
}

void athread_svrasr_host (float* yl, float* b, int m, int ni, int nk, int nj) {

  struct SvrasrParam para;
  para.yl = yl;
  para.b  = b;
  para.m  = m;
  para.ni = ni;
  para.nk = nk;
  para.nj = nj;

	athread_spawn (athread_svrasr_slave, &para);
	athread_join ();
  return ;
}