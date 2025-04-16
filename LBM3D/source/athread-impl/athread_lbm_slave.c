#include "../include/params.h"
#include "../include/athread_lbm_params.h"

#include "slave.h"

#include <math.h>

#define INDEX_3D(l1,l2,i0,i1,i2) ((i0)*(l1)*(l2) + (i1)*(l2) + (i2))
#define INDEX_4D(l1,l2,l3,i0,i1,i2,i3) ((i0)*(l1)*(l2)*(l3) + (i1)*(l2)*(l3) + (i2)*(l3) + (i3))

__thread_local int buf_reduce_int __attribute__ ((aligned(64)));

void collide_stream_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;
	double* fA = athread_params->fA;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int D3 = athread_params->D3;

  double tau = athread_params->params->tau;
  double tau_inv = 1. / athread_params->params->tau;
  double wo = 1. / 3.;
  double ws = 1. / 18.;
  double wd = 1. / 36.;
  double omtau_inv = 1.0 - tau_inv;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = ((start+len) < (D1-1)) ? (start+len) : (D1-1);

  // if (athread_tid == 0) {
  //   printf ("D1 = %d len = %d\n", D1, len);
  // }
  // athread_ssync_array ();
  // printf ("mytid = %d, start = %d, end = %d\n", athread_tid, start, end);
	for (int i = start; i < end; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
			for (int k = 1; k < D3 - 1; ++k) {
        double f0  = fA[INDEX_4D(nx, nz, 19, i, j, k, 0)];
        double f1  = fA[INDEX_4D(nx, nz, 19, i, j, k, 1)];
        double f2  = fA[INDEX_4D(nx, nz, 19, i, j, k, 2)];
        double f3  = fA[INDEX_4D(nx, nz, 19, i, j, k, 3)];
        double f4  = fA[INDEX_4D(nx, nz, 19, i, j, k, 4)];
        double f5  = fA[INDEX_4D(nx, nz, 19, i, j, k, 5)];
        double f6  = fA[INDEX_4D(nx, nz, 19, i, j, k, 6)];
        double f7  = fA[INDEX_4D(nx, nz, 19, i, j, k, 7)];
        double f8  = fA[INDEX_4D(nx, nz, 19, i, j, k, 8)];
        double f9  = fA[INDEX_4D(nx, nz, 19, i, j, k, 9)];
        double f10 = fA[INDEX_4D(nx, nz, 19, i, j, k, 10)];
        double f11 = fA[INDEX_4D(nx, nz, 19, i, j, k, 11)];
        double f12 = fA[INDEX_4D(nx, nz, 19, i, j, k, 12)];
        double f13 = fA[INDEX_4D(nx, nz, 19, i, j, k, 13)];
        double f14 = fA[INDEX_4D(nx, nz, 19, i, j, k, 14)];
        double f15 = fA[INDEX_4D(nx, nz, 19, i, j, k, 15)];
        double f16 = fA[INDEX_4D(nx, nz, 19, i, j, k, 16)];
        double f17 = fA[INDEX_4D(nx, nz, 19, i, j, k, 17)];
        double f18 = fA[INDEX_4D(nx, nz, 19, i, j, k, 18)];
   
        // compute density
        double density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
        double density_inv = 1. / density;
   
        // compute velocities
        double u = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16)) * density_inv;
        double v = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18)) * density_inv;
        double w = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17)) * density_inv;
   
        // update distribubtions
        double tw0r = wo * density * tau_inv; // w[0]*rho
        double twsr = ws * density * tau_inv; // w[1-6]*rho
        double twdr = wd * density * tau_inv; // w[7-18]*rho
   
        double tu = 3.0 * u;
        double tv = 3.0 * v;
        double tw = 3.0 * w;
   
        double a = 1.0 - 1.5 * (u * u + v * v + w * w);
   
        double udot7 = tu + tv;
        double udot9 = tu + tw;
        double udot11 = tv + tw;
        double udot13 = tu - tv;
        double udot15 = tu - tw;
        double udot17 = tv - tw;
   
        double feq0  = tw0r * a;
        double feq1  = twsr * (a + tu * (1.0 + 0.5 * tu));
        double feq2  = twsr * (a - tu * (1.0 - 0.5 * tu));
        double feq3  = twsr * (a + tv * (1.0 + 0.5 * tv));
        double feq4  = twsr * (a - tv * (1.0 - 0.5 * tv));
        double feq5  = twsr * (a + tw * (1.0 + 0.5 * tw));
        double feq6  = twsr * (a - tw * (1.0 - 0.5 * tw));
        double feq7  = twdr * (a + udot7 * (1.0 + 0.5 * udot7));
        double feq8  = twdr * (a - udot7 * (1.0 - 0.5 * udot7));
        double feq9  = twdr * (a + udot9 * (1.0 + 0.5 * udot9));
        double feq10 = twdr * (a - udot9 * (1.0 - 0.5 * udot9));
        double feq11 = twdr * (a + udot11 * (1.0 + 0.5 * udot11));
        double feq12 = twdr * (a - udot11 * (1.0 - 0.5 * udot11));
        double feq13 = twdr * (a + udot13 * (1.0 + 0.5 * udot13));
        double feq14 = twdr * (a - udot13 * (1.0 - 0.5 * udot13));
        double feq15 = twdr * (a + udot15 * (1.0 + 0.5 * udot15));
        double feq16 = twdr * (a - udot15 * (1.0 - 0.5 * udot15));
        double feq17 = twdr * (a + udot17 * (1.0 + 0.5 * udot17));
        double feq18 = twdr * (a - udot17 * (1.0 - 0.5 * udot17));
   
        double fB0  = omtau_inv * f0 + feq0;
        double fB1  = omtau_inv * f1 + feq1;
        double fB2  = omtau_inv * f2 + feq2;
        double fB3  = omtau_inv * f3 + feq3;
        double fB4  = omtau_inv * f4 + feq4;
        double fB5  = omtau_inv * f5 + feq5;
        double fB6  = omtau_inv * f6 + feq6;
        double fB7  = omtau_inv * f7 + feq7;
        double fB8  = omtau_inv * f8 + feq8;
        double fB9  = omtau_inv * f9 + feq9;
        double fB10 = omtau_inv * f10 + feq10;
        double fB11 = omtau_inv * f11 + feq11;
        double fB12 = omtau_inv * f12 + feq12;
        double fB13 = omtau_inv * f13 + feq13;
        double fB14 = omtau_inv * f14 + feq14;
        double fB15 = omtau_inv * f15 + feq15;
        double fB16 = omtau_inv * f16 + feq16;
        double fB17 = omtau_inv * f17 + feq17;
        double fB18 = omtau_inv * f18 + feq18;
   
        // stream updated distributions
        fB[INDEX_4D(nx, nz, 19, i, j, k, 0)] = fB0;
        fB[INDEX_4D(nx, nz, 19, i, j + 1, k, 1)] = fB1;
        fB[INDEX_4D(nx, nz, 19, i, j - 1, k, 2)] = fB2;
        fB[INDEX_4D(nx, nz, 19, i + 1, j, k, 3)] = fB3;
        fB[INDEX_4D(nx, nz, 19, i - 1, j, k, 4)] = fB4;
        fB[INDEX_4D(nx, nz, 19, i, j, k + 1, 5)] = fB5;
        fB[INDEX_4D(nx, nz, 19, i, j, k - 1, 6)] = fB6;
        fB[INDEX_4D(nx, nz, 19, i + 1, j + 1, k, 7)] = fB7;
        fB[INDEX_4D(nx, nz, 19, i - 1, j - 1, k, 8)] = fB8;
        fB[INDEX_4D(nx, nz, 19, i, j + 1, k + 1, 9)] = fB9;
        fB[INDEX_4D(nx, nz, 19, i, j - 1, k - 1, 10)] = fB10;
        fB[INDEX_4D(nx, nz, 19, i + 1, j, k + 1, 11)] = fB11;
        fB[INDEX_4D(nx, nz, 19, i - 1, j, k - 1, 12)] = fB12;
        fB[INDEX_4D(nx, nz, 19, i - 1, j + 1, k, 13)] = fB13;
        fB[INDEX_4D(nx, nz, 19, i + 1, j - 1, k, 14)] = fB14;
        fB[INDEX_4D(nx, nz, 19, i, j + 1, k - 1, 15)] = fB15;
        fB[INDEX_4D(nx, nz, 19, i, j - 1, k + 1, 16)] = fB16;
        fB[INDEX_4D(nx, nz, 19, i + 1, j, k - 1, 17)] = fB17;
        fB[INDEX_4D(nx, nz, 19, i - 1, j, k + 1, 18)] = fB18;
			}
		}
	}
	return ;
}

void bb_left_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D3 = athread_params->D3;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D1-1) ? (start+len) : (D1-1);

	const int j = 1;
	for (int i = start; i < end; ++i) {
		for (int k = 1; k < D3 - 1; ++k) {
      // 2, 8, 10, 14, 16
      fB[INDEX_4D(nx, nz, 19, i, j, k, 1)] = fB[INDEX_4D(nx, nz, 19, i, j-1, k, 2)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 7)] = fB[INDEX_4D(nx, nz, 19, i-1, j-1, k, 8)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 9)] = fB[INDEX_4D(nx, nz, 19, i, j-1, k-1, 10)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 13)] = fB[INDEX_4D(nx, nz, 19, i+1, j-1, k, 14)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 15)] = fB[INDEX_4D(nx, nz, 19, i, j-1, k+1, 16)];
		}
	}
	return ;
}

void bb_right_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int D3 = athread_params->D3;
	const int j = D2 - 2;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D1-1) ? (start+len) : (D1-1);

	for (int i = start; i < end; ++i) {
		for (int k = 1; k < D3 - 1; ++k) {
    	// 1, 7, 9, 13, 15
      fB[INDEX_4D(nx, nz, 19, i, j, k, 2)] = fB[INDEX_4D(nx, nz, 19, i, j+1, k, 1)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 8)] = fB[INDEX_4D(nx, nz, 19, i+1, j+1, k, 7)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 10)] = fB[INDEX_4D(nx, nz, 19, i, j+1, k+1, 9)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 14)] = fB[INDEX_4D(nx, nz, 19, i-1, j+1, k, 13)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 16)] = fB[INDEX_4D(nx, nz, 19, i, j+1, k-1, 15)];
		}
	}
	return ;
}

void bb_front_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  // int ny = params->ny;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int D3 = athread_params->D3;
	int k = D3 - 2;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D1-1) ? (start+len) : (D1-1);

	for (int i = start; i < end; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
    	// 5, 9, 11, 16, 18
      fB[INDEX_4D(nx, nz, 19, i, j, k, 6)] = fB[INDEX_4D(nx, nz, 19, i, j, k+1, 5)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 10)] = fB[INDEX_4D(nx, nz, 19, i, j+1, k+1, 9)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 12)] = fB[INDEX_4D(nx, nz, 19, i+1, j, k+1, 11)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 15)] = fB[INDEX_4D(nx, nz, 19, i, j-1, k+1, 16)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 17)] = fB[INDEX_4D(nx, nz, 19, i-1, j, k+1, 18)];
		}
	}
	return ;
}

void bb_back_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  // int ny = params->ny;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int k = 1;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D1-1) ? (start+len) : (D1-1);

	for (int i = start; i < end; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
    	// 6, 10, 12, 15, 17
      fB[INDEX_4D(nx, nz, 19, i, j, k, 5)] = fB[INDEX_4D(nx, nz, 19, i, j, k-1, 6)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 9)] = fB[INDEX_4D(nx, nz, 19, i, j-1, k-1, 10)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 11)] = fB[INDEX_4D(nx, nz, 19, i-1, j, k-1, 12)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 16)] = fB[INDEX_4D(nx, nz, 19, i, j+1, k-1, 15)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 18)] = fB[INDEX_4D(nx, nz, 19, i+1, j, k-1, 17)];
		}
	}
	return ;
}

void bb_bottom_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;

	int D2 = athread_params->D2;
	int D3 = athread_params->D3;

	int i = 1;

	// Range: [1, D2-1)
	int len = D2 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D2-1) ? (start+len) : (D2-1);

	for (int j = start; j < end; ++j) {
		for (int k = 1; k < D3 - 1; ++k) {
    	// 4, 8 , 12, 13, 18
      fB[INDEX_4D(nx, nz, 19, i, j, k, 3)] = fB[INDEX_4D(nx, nz, 19, i-1, j, k, 4)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 7)] = fB[INDEX_4D(nx, nz, 19, i-1, j-1, k, 8)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 11)] = fB[INDEX_4D(nx, nz, 19, i-1, j, k-1, 12)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 14)] = fB[INDEX_4D(nx, nz, 19, i-1, j+1, k, 13)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 17)] = fB[INDEX_4D(nx, nz, 19, i-1, j, k+1, 18)];
		}
	}
	return ;
}

void bb_top_slave (struct AthreadParams* athread_params) {
	double* fB = athread_params->fB;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int D3 = athread_params->D3;

  double wd = 1. / 36;
  double u = params->lid_u;
  double w = params->lid_w;
	int i = D1 - 2;

  double rhs_7 = 6.0 * wd * u;
  double rhs_11 = 6.0 * wd * w;
  double rhs_14 = -6.0 * wd * u;
  double rhs_17 = -6.0 * wd * w;

	// Range: [1, D2-1)
	int len = D2 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D2-1) ? (start+len) : (D2-1);

	for (int j = start; j < end; ++j) {
		for (int k = 1; k < D3 - 1; ++k) {
    	// 3, 7, 11, 14, 17
      fB[INDEX_4D(nx, nz, 19, i, j, k, 4)] = fB[INDEX_4D(nx, nz, 19, i+1, j, k, 3)];
      fB[INDEX_4D(nx, nz, 19, i, j, k, 8)] = fB[INDEX_4D(nx, nz, 19, i+1, j+1, k, 7)] - rhs_7;
      fB[INDEX_4D(nx, nz, 19, i, j, k, 12)] = fB[INDEX_4D(nx, nz, 19, i+1, j, k+1, 11)] - rhs_11;
      fB[INDEX_4D(nx, nz, 19, i, j, k, 13)] = fB[INDEX_4D(nx, nz, 19, i+1, j-1, k, 14)] - rhs_14;
      fB[INDEX_4D(nx, nz, 19, i, j, k, 18)] = fB[INDEX_4D(nx, nz, 19, i+1, j, k-1, 17)] - rhs_17;
		}
	}
	return ;
}

void compute_macroscopic_slave (struct AthreadParams* athread_params) {

	double* f = athread_params->fB;

	double* u = athread_params->u;
	double* v = athread_params->v;
	double* w = athread_params->w;
	double* rho = athread_params->rho;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int D3 = athread_params->D3;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = (start+len) < (D1-1) ? (start+len) : (D1-1);

	for (int i = start; i < end; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
			for (int k = 1; k < D3 - 1; ++k) {
        double f0  = f[INDEX_4D(nx, nz, 19, i, j, k, 0)];
        double f1  = f[INDEX_4D(nx, nz, 19, i, j, k, 1)];
        double f2  = f[INDEX_4D(nx, nz, 19, i, j, k, 2)];
        double f3  = f[INDEX_4D(nx, nz, 19, i, j, k, 3)];
        double f4  = f[INDEX_4D(nx, nz, 19, i, j, k, 4)];
        double f5  = f[INDEX_4D(nx, nz, 19, i, j, k, 5)];
        double f6  = f[INDEX_4D(nx, nz, 19, i, j, k, 6)];
        double f7  = f[INDEX_4D(nx, nz, 19, i, j, k, 7)];
        double f8  = f[INDEX_4D(nx, nz, 19, i, j, k, 8)];
        double f9  = f[INDEX_4D(nx, nz, 19, i, j, k, 9)];
        double f10 = f[INDEX_4D(nx, nz, 19, i, j, k, 10)];
        double f11 = f[INDEX_4D(nx, nz, 19, i, j, k, 11)];
        double f12 = f[INDEX_4D(nx, nz, 19, i, j, k, 12)];
        double f13 = f[INDEX_4D(nx, nz, 19, i, j, k, 13)];
        double f14 = f[INDEX_4D(nx, nz, 19, i, j, k, 14)];
        double f15 = f[INDEX_4D(nx, nz, 19, i, j, k, 15)];
        double f16 = f[INDEX_4D(nx, nz, 19, i, j, k, 16)];
        double f17 = f[INDEX_4D(nx, nz, 19, i, j, k, 17)];
        double f18 = f[INDEX_4D(nx, nz, 19, i, j, k, 18)];
        // compute density
        double density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
        double density_inv = 1. / density;

        // compute velocities
        double utmp = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16)) * density_inv;
        double vtmp = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18)) * density_inv;
        double wtmp = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17)) * density_inv;

        // write to device
        rho[INDEX_3D(nx, nz, i, j, k)] = density;
        u  [INDEX_3D(nx, nz, i, j, k)] = utmp;
        v  [INDEX_3D(nx, nz, i, j, k)] = vtmp;
        w  [INDEX_3D(nx, nz, i, j, k)] = wtmp;
			}
		}
	}
	return ;
}

void check_if_steady_state_slave (struct AthreadParams* athread_params) {
	int converged = 1;

	double* fA = athread_params->fA;

	double* u = athread_params->u;
	double* v = athread_params->v;
	double* w = athread_params->w;
	double* rho = athread_params->rho;

	struct Params* params = athread_params->params;

  int nx = params->nx;
  int nz = params->nz;
	double tol = params->tol;

	int D1 = athread_params->D1;
	int D2 = athread_params->D2;
	int D3 = athread_params->D3;

	// Range: [1, D1-1)
	int len = D1 - 2;
	len = (len + 63) / 64;
	int start = athread_tid * len + 1;
	int end = ((start+len) < (D1-1)) ? (start+len) : (D1-1);
	for (int i = start; i < end; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
			for (int k = 1; k < D3 - 1; ++k) {
    		// load distributions A
        double f0  = fA[INDEX_4D(nx, nz, 19, i, j, k, 0)];
        double f1  = fA[INDEX_4D(nx, nz, 19, i, j, k, 1)];
        double f2  = fA[INDEX_4D(nx, nz, 19, i, j, k, 2)];
        double f3  = fA[INDEX_4D(nx, nz, 19, i, j, k, 3)];
        double f4  = fA[INDEX_4D(nx, nz, 19, i, j, k, 4)];
        double f5  = fA[INDEX_4D(nx, nz, 19, i, j, k, 5)];
        double f6  = fA[INDEX_4D(nx, nz, 19, i, j, k, 6)];
        double f7  = fA[INDEX_4D(nx, nz, 19, i, j, k, 7)];
        double f8  = fA[INDEX_4D(nx, nz, 19, i, j, k, 8)];
        double f9  = fA[INDEX_4D(nx, nz, 19, i, j, k, 9)];
        double f10 = fA[INDEX_4D(nx, nz, 19, i, j, k, 10)];
        double f11 = fA[INDEX_4D(nx, nz, 19, i, j, k, 11)];
        double f12 = fA[INDEX_4D(nx, nz, 19, i, j, k, 12)];
        double f13 = fA[INDEX_4D(nx, nz, 19, i, j, k, 13)];
        double f14 = fA[INDEX_4D(nx, nz, 19, i, j, k, 14)];
        double f15 = fA[INDEX_4D(nx, nz, 19, i, j, k, 15)];
        double f16 = fA[INDEX_4D(nx, nz, 19, i, j, k, 16)];
        double f17 = fA[INDEX_4D(nx, nz, 19, i, j, k, 17)];
        double f18 = fA[INDEX_4D(nx, nz, 19, i, j, k, 18)];

        // compute velocities A
        double umom_A = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16));
        double vmom_A = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18));
        double wmom_A = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17));

        // load macroscopic properties B
        double rho_B = rho[INDEX_3D(nx, nz, i, j, k)];
        double umom_B = u[INDEX_3D(nx, nz, i, j, k)] * rho_B;
        double vmom_B = v[INDEX_3D(nx, nz, i, j, k)] * rho_B;
        double wmom_B = w[INDEX_3D(nx, nz, i, j, k)] * rho_B;

    		// convergence criteria
        int converged_momu = fabs(umom_B - umom_A) < fmax(tol * fabs(umom_A), 1e-12);
        int converged_momv = fabs(vmom_B - vmom_A) < fmax(tol * fabs(vmom_A), 1e-12);
        int converged_momw = fabs(wmom_B - wmom_A) < fmax(tol * fabs(wmom_A), 1e-12);
        converged &= converged_momu & converged_momv & converged_momw;
			}
		}
	}
	// reduce converged
  athread_ssync_array ();
  athread_redurt (&converged, &converged, 1, athread_int, OP_and, &buf_reduce_int, 1);
  if (athread_tid == 0) {
		*(athread_params->converged) = converged;
  }
  athread_ssync_array ();
	return ;
}

#undef INDEX_3D
#undef INDEX_4D