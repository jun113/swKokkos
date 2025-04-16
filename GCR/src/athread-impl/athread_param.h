#ifndef SRC_ATHREAD_IMPL_ATHREAD_PARAM_H_
#define SRC_ATHREAD_IMPL_ATHREAD_PARAM_H_

struct KernelParam1 {
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int NX;
  int NY;
  int NZ;
  float* p;
  double* x0;
};

struct KernelParam2 {
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* r;
  float* f0;
};

struct KernelParam3 {
  int m;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* p;
  float* r;
};

struct KernelParam4 {
  int iter_max;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int NX;
  int NY;
  int NZ;
  float* p;
};

struct KernelParam5 {
  int m;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  double ac;
  float* p;
  double* x0;
};

struct KernelParam6 {
  int m;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  double ac;
  float* r;
  float* ap;
};

struct KernelParam7 {
  int m;
  double* b;
  double* c2;
  double* aps;
};

struct KernelParam8 {
  int l;
  int m;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* p;
  double* b;
};

struct KernelParam9 {
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* r;
  float* f0;
  float* ap;
};

struct KernelReudceDParam {
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* r;
  double* d;
};

struct KernelReudceC11Param {
  int m;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* r;
  float* ap;
  double* c11;
};

struct KernelReudceC12Param {
  int m;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* ap;
  double* c12;
};

struct KernelReudceClParam {
  int l;
  int its;
  int kts;
  int jts;
  int ite;
  int kte;
  int jte;
  int ims;
  int ime;
  int jms;
  int jend;
  int NX;
  int NY;
  int NZ;
  float* ar;
  float* ap;
  double* cl;
};

struct KernelReudceNorm2Param {
  int len;
  float*  r;
  double* result;
};

struct MatrixProParam {
  int its;
  int ite;
  int jts;
  int jte;
  int kts;
  int kte;
  int jend;
  float* a;
  float* b;
  float* c;
};

struct MatrixProDoubleParam {
  int its;
  int ite;
  int jts;
  int jte;
  int kts;
  int kte;
  int jend;
  float* a;
  double* b;
  float* c;
};

struct SvrasrParam {
  float* yl;
  float* b;
  int m;
  int ni;
  int nk;
  int nj;
};

#endif  // SRC_ATHREAD_IMPL_ATHREAD_PARAM_H_