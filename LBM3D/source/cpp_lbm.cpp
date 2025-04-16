#include "include/typedefs.h"
#include "include/params.h"

#include <cmath>

inline int index_3d (const int &l1, const int &l2,
		const int &i0, const int &i1, const int &i2) {
	return i0*(l1*l2) + i1*(l2) + i2;
}

inline int index_4d (const int &l1, const int &l2, const int &l3,
		const int &i0, const int &i1, const int &i2, const int &i3) {
	return i0*(l1*l2*l3) + i1*(l2*l3) + i2*l3 + i3;
}

static void cpp_collide_stream (double* const fB, double* const fA, 
		const Params &params, const int &D1, const int &D2, const int &D3);

static void cpp_bb_left (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3);

static void cpp_bb_right (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3);

static void cpp_bb_front (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3);

static void cpp_bb_back (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3);

static void cpp_bb_bottom (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3);

static void cpp_bb_top (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3);

static void cpp_compute_macroscopic (double* const f,
		double* const u, double* const v, double* const w, double* const rho,
				const Params &params, const int &D1, const int &D2, const int &D3);

static void cpp_check_if_steady_state (double* const fA, 
		double* const u, double* const v, double* const w, double* const rho,
				const Params &params, int &converged, const int &D1, const int &D2, const int &D3);

void cpp_lbm_update (double* const fB, double* const fA, 
		double* const u, double* const v, double* const w, double* const rho,
				const Params &params, const int &step, int &converged, 
						const int &D1, const int &D2, const int &D3) {
	// collide_stream
	cpp_collide_stream (fB, fA, params, D1, D2, D3);
	// bb_left
	// cpp_bb_left (fB, params, D1, D2, D3);
	// bb_right
	// cpp_bb_right (fB, params, D1, D2, D3);
	// bb_front
	// cpp_bb_front (fB, params, D1, D2, D3);
	// bb_back
	// cpp_bb_back (fB, params, D1, D2, D3);
	// bb_bottom
	// cpp_bb_bottom (fB, params, D1, D2, D3);
	// bb_top
	// cpp_bb_top (fB, params, D1, D2, D3);
  if ((step + 1) % params.output_rate == 0) {
		// compute_macroscopic
		// cpp_compute_macroscopic (fB, u, v, w, rho, params, D1, D2, D3);
		// check_if_steady_state
    converged = 1;
		// cpp_check_if_steady_state (fA, u, v, w, rho, params, converged, D1, D2, D3);
	}
	return ;
}

static void cpp_collide_stream (double* const fB, double* const fA, 
		const Params &params, const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;

  const Double tau = params.tau;
  const Double tau_inv = 1. / params.tau;
  const Double wo = 1. / 3.;
  const Double ws = 1. / 18.;
  const Double wd = 1. / 36.;
  const Double omtau_inv = 1.0 - tau_inv;

	for (int i = 1; i < D1 - 1; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
			for (int k = 1; k < D3 - 1; ++k) {
        // load distributions
        Double f0  = fA[index_4d(nx, nz, 19, i, j, k, 0)];
        Double f1  = fA[index_4d(nx, nz, 19, i, j, k, 1)];
        Double f2  = fA[index_4d(nx, nz, 19, i, j, k, 2)];
        Double f3  = fA[index_4d(nx, nz, 19, i, j, k, 3)];
        Double f4  = fA[index_4d(nx, nz, 19, i, j, k, 4)];
        Double f5  = fA[index_4d(nx, nz, 19, i, j, k, 5)];
        Double f6  = fA[index_4d(nx, nz, 19, i, j, k, 6)];
        Double f7  = fA[index_4d(nx, nz, 19, i, j, k, 7)];
        Double f8  = fA[index_4d(nx, nz, 19, i, j, k, 8)];
        Double f9  = fA[index_4d(nx, nz, 19, i, j, k, 9)];
        Double f10 = fA[index_4d(nx, nz, 19, i, j, k, 10)];
        Double f11 = fA[index_4d(nx, nz, 19, i, j, k, 11)];
        Double f12 = fA[index_4d(nx, nz, 19, i, j, k, 12)];
        Double f13 = fA[index_4d(nx, nz, 19, i, j, k, 13)];
        Double f14 = fA[index_4d(nx, nz, 19, i, j, k, 14)];
        Double f15 = fA[index_4d(nx, nz, 19, i, j, k, 15)];
        Double f16 = fA[index_4d(nx, nz, 19, i, j, k, 16)];
        Double f17 = fA[index_4d(nx, nz, 19, i, j, k, 17)];
        Double f18 = fA[index_4d(nx, nz, 19, i, j, k, 18)];
   
        // compute density
        Double density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
        Double density_inv = 1. / density;
   
        // compute velocities
        Double u = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16)) * density_inv;
        Double v = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18)) * density_inv;
        Double w = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17)) * density_inv;
   
        // update distribubtions
        Double tw0r = wo * density * tau_inv; // w[0]*rho
        Double twsr = ws * density * tau_inv; // w[1-6]*rho
        Double twdr = wd * density * tau_inv; // w[7-18]*rho
   
        Double tu = 3.0 * u;
        Double tv = 3.0 * v;
        Double tw = 3.0 * w;
   
        Double a = 1.0 - 1.5 * (u * u + v * v + w * w);
   
        Double udot7 = tu + tv;
        Double udot9 = tu + tw;
        Double udot11 = tv + tw;
        Double udot13 = tu - tv;
        Double udot15 = tu - tw;
        Double udot17 = tv - tw;
   
        Double feq0  = tw0r * a;
        Double feq1  = twsr * (a + tu * (1.0 + 0.5 * tu));
        Double feq2  = twsr * (a - tu * (1.0 - 0.5 * tu));
        Double feq3  = twsr * (a + tv * (1.0 + 0.5 * tv));
        Double feq4  = twsr * (a - tv * (1.0 - 0.5 * tv));
        Double feq5  = twsr * (a + tw * (1.0 + 0.5 * tw));
        Double feq6  = twsr * (a - tw * (1.0 - 0.5 * tw));
        Double feq7  = twdr * (a + udot7 * (1.0 + 0.5 * udot7));
        Double feq8  = twdr * (a - udot7 * (1.0 - 0.5 * udot7));
        Double feq9  = twdr * (a + udot9 * (1.0 + 0.5 * udot9));
        Double feq10 = twdr * (a - udot9 * (1.0 - 0.5 * udot9));
        Double feq11 = twdr * (a + udot11 * (1.0 + 0.5 * udot11));
        Double feq12 = twdr * (a - udot11 * (1.0 - 0.5 * udot11));
        Double feq13 = twdr * (a + udot13 * (1.0 + 0.5 * udot13));
        Double feq14 = twdr * (a - udot13 * (1.0 - 0.5 * udot13));
        Double feq15 = twdr * (a + udot15 * (1.0 + 0.5 * udot15));
        Double feq16 = twdr * (a - udot15 * (1.0 - 0.5 * udot15));
        Double feq17 = twdr * (a + udot17 * (1.0 + 0.5 * udot17));
        Double feq18 = twdr * (a - udot17 * (1.0 - 0.5 * udot17));
   
        Double fB0  = omtau_inv * f0 + feq0;
        Double fB1  = omtau_inv * f1 + feq1;
        Double fB2  = omtau_inv * f2 + feq2;
        Double fB3  = omtau_inv * f3 + feq3;
        Double fB4  = omtau_inv * f4 + feq4;
        Double fB5  = omtau_inv * f5 + feq5;
        Double fB6  = omtau_inv * f6 + feq6;
        Double fB7  = omtau_inv * f7 + feq7;
        Double fB8  = omtau_inv * f8 + feq8;
        Double fB9  = omtau_inv * f9 + feq9;
        Double fB10 = omtau_inv * f10 + feq10;
        Double fB11 = omtau_inv * f11 + feq11;
        Double fB12 = omtau_inv * f12 + feq12;
        Double fB13 = omtau_inv * f13 + feq13;
        Double fB14 = omtau_inv * f14 + feq14;
        Double fB15 = omtau_inv * f15 + feq15;
        Double fB16 = omtau_inv * f16 + feq16;
        Double fB17 = omtau_inv * f17 + feq17;
        Double fB18 = omtau_inv * f18 + feq18;
   
        // stream updated distributions
        fB[index_4d(nx, nz, 19, i, j, k, 0)] = fB0;
        fB[index_4d(nx, nz, 19, i, j + 1, k, 1)] = fB1;
        fB[index_4d(nx, nz, 19, i, j - 1, k, 2)] = fB2;
        fB[index_4d(nx, nz, 19, i + 1, j, k, 3)] = fB3;
        fB[index_4d(nx, nz, 19, i - 1, j, k, 4)] = fB4;
        fB[index_4d(nx, nz, 19, i, j, k + 1, 5)] = fB5;
        fB[index_4d(nx, nz, 19, i, j, k - 1, 6)] = fB6;
        fB[index_4d(nx, nz, 19, i + 1, j + 1, k, 7)] = fB7;
        fB[index_4d(nx, nz, 19, i - 1, j - 1, k, 8)] = fB8;
        fB[index_4d(nx, nz, 19, i, j + 1, k + 1, 9)] = fB9;
        fB[index_4d(nx, nz, 19, i, j - 1, k - 1, 10)] = fB10;
        fB[index_4d(nx, nz, 19, i + 1, j, k + 1, 11)] = fB11;
        fB[index_4d(nx, nz, 19, i - 1, j, k - 1, 12)] = fB12;
        fB[index_4d(nx, nz, 19, i - 1, j + 1, k, 13)] = fB13;
        fB[index_4d(nx, nz, 19, i + 1, j - 1, k, 14)] = fB14;
        fB[index_4d(nx, nz, 19, i, j + 1, k - 1, 15)] = fB15;
        fB[index_4d(nx, nz, 19, i, j - 1, k + 1, 16)] = fB16;
        fB[index_4d(nx, nz, 19, i + 1, j, k - 1, 17)] = fB17;
        fB[index_4d(nx, nz, 19, i - 1, j, k + 1, 18)] = fB18;
			}
		}
	}
	return ;
}

static void cpp_bb_left (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;
	const int j = 1;
	for (int i = 1; i < D1 - 1; ++i) {
		for (int k = 1; k < D3 - 1; ++k) {
      // 2, 8, 10, 14, 16
      fB[index_4d(nx, nz, 19, i, j, k, 1)] = fB[index_4d(nx, nz, 19, i, j-1, k, 2)];
      fB[index_4d(nx, nz, 19, i, j, k, 7)] = fB[index_4d(nx, nz, 19, i-1, j-1, k, 8)];
      fB[index_4d(nx, nz, 19, i, j, k, 9)] = fB[index_4d(nx, nz, 19, i, j-1, k-1, 10)];
      fB[index_4d(nx, nz, 19, i, j, k, 13)] = fB[index_4d(nx, nz, 19, i+1, j-1, k, 14)];
      fB[index_4d(nx, nz, 19, i, j, k, 15)] = fB[index_4d(nx, nz, 19, i, j-1, k+1, 16)];
		}
	}
	return ;
}

static void cpp_bb_right (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;
	const int j = D2 - 2;
	for (int i = 1; i < D1 - 1; ++i) {
		for (int k = 1; k < D3 - 1; ++k) {
    	// 1, 7, 9, 13, 15
      fB[index_4d(nx, nz, 19, i, j, k, 2)] = fB[index_4d(nx, nz, 19, i, j+1, k, 1)];
      fB[index_4d(nx, nz, 19, i, j, k, 8)] = fB[index_4d(nx, nz, 19, i+1, j+1, k, 7)];
      fB[index_4d(nx, nz, 19, i, j, k, 10)] = fB[index_4d(nx, nz, 19, i, j+1, k+1, 9)];
      fB[index_4d(nx, nz, 19, i, j, k, 14)] = fB[index_4d(nx, nz, 19, i-1, j+1, k, 13)];
      fB[index_4d(nx, nz, 19, i, j, k, 16)] = fB[index_4d(nx, nz, 19, i, j+1, k-1, 15)];
		}
	}
	return ;
}

static void cpp_bb_front (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;
	const int k = D3 - 2;
	for (int i = 1; i < D1 - 1; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
    	// 5, 9, 11, 16, 18
      fB[index_4d(nx, nz, 19, i, j, k, 6)] = fB[index_4d(nx, nz, 19, i, j, k+1, 5)];
      fB[index_4d(nx, nz, 19, i, j, k, 10)] = fB[index_4d(nx, nz, 19, i, j+1, k+1, 9)];
      fB[index_4d(nx, nz, 19, i, j, k, 12)] = fB[index_4d(nx, nz, 19, i+1, j, k+1, 11)];
      fB[index_4d(nx, nz, 19, i, j, k, 15)] = fB[index_4d(nx, nz, 19, i, j-1, k+1, 16)];
      fB[index_4d(nx, nz, 19, i, j, k, 17)] = fB[index_4d(nx, nz, 19, i-1, j, k+1, 18)];
		}
	}
	return ;
}

static void cpp_bb_back (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;
	const int k = 1;
	for (int i = 1; i < D1 - 1; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
    	// 6, 10, 12, 15, 17
      fB[index_4d(nx, nz, 19, i, j, k, 5)] = fB[index_4d(nx, nz, 19, i, j, k-1, 6)];
      fB[index_4d(nx, nz, 19, i, j, k, 9)] = fB[index_4d(nx, nz, 19, i, j-1, k-1, 10)];
      fB[index_4d(nx, nz, 19, i, j, k, 11)] = fB[index_4d(nx, nz, 19, i-1, j, k-1, 12)];
      fB[index_4d(nx, nz, 19, i, j, k, 16)] = fB[index_4d(nx, nz, 19, i, j+1, k-1, 15)];
      fB[index_4d(nx, nz, 19, i, j, k, 18)] = fB[index_4d(nx, nz, 19, i+1, j, k-1, 17)];
		}
	}
	return ;
}

static void cpp_bb_bottom (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;
	const int i = 1;
	for (int j = 1; j < D2 - 1; ++j) {
		for (int k = 1; k < D3 - 1; ++k) {
    	// 4, 8 , 12, 13, 18
      fB[index_4d(nx, nz, 19, i, j, k, 3)] = fB[index_4d(nx, nz, 19, i-1, j, k, 4)];
      fB[index_4d(nx, nz, 19, i, j, k, 7)] = fB[index_4d(nx, nz, 19, i-1, j-1, k, 8)];
      fB[index_4d(nx, nz, 19, i, j, k, 11)] = fB[index_4d(nx, nz, 19, i-1, j, k-1, 12)];
      fB[index_4d(nx, nz, 19, i, j, k, 14)] = fB[index_4d(nx, nz, 19, i-1, j+1, k, 13)];
      fB[index_4d(nx, nz, 19, i, j, k, 17)] = fB[index_4d(nx, nz, 19, i-1, j, k+1, 18)];
		}
	}
	return ;
}

static void cpp_bb_top (double* const fB, const Params &params, 
		const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;
  const Double wd = 1. / 36;
  const Double u = params.lid_u;
  const Double w = params.lid_w;
	const int i = D1 - 2;

  Double rhs_7 = 6.0 * wd * u;
  Double rhs_11 = 6.0 * wd * w;
  Double rhs_14 = -6.0 * wd * u;
  Double rhs_17 = -6.0 * wd * w;

	for (int j = 1; j < D2 - 1; ++j) {
		for (int k = 1; k < D3 - 1; ++k) {
    	// 3, 7, 11, 14, 17
      fB[index_4d(nx, nz, 19, i, j, k, 4)] = fB[index_4d(nx, nz, 19, i+1, j, k, 3)];
      fB[index_4d(nx, nz, 19, i, j, k, 8)] = fB[index_4d(nx, nz, 19, i+1, j+1, k, 7)] - rhs_7;
      fB[index_4d(nx, nz, 19, i, j, k, 12)] = fB[index_4d(nx, nz, 19, i+1, j, k+1, 11)] - rhs_11;
      fB[index_4d(nx, nz, 19, i, j, k, 13)] = fB[index_4d(nx, nz, 19, i+1, j-1, k, 14)] - rhs_14;
      fB[index_4d(nx, nz, 19, i, j, k, 18)] = fB[index_4d(nx, nz, 19, i+1, j, k-1, 17)] - rhs_17;
		}
	}
	return ;
}

static void cpp_compute_macroscopic (double* const f,
		double* const u, double* const v, double* const w, double* const rho,
				const Params &params, const int &D1, const int &D2, const int &D3) {
  const int nx = params.nx;
  const int nz = params.nz;

	for (int i = 1; i < D1 - 1; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
			for (int k = 1; k < D3 - 1; ++k) {
        Double f0  = f[index_4d(nx, nz, 19, i, j, k, 0)];
        Double f1  = f[index_4d(nx, nz, 19, i, j, k, 1)];
        Double f2  = f[index_4d(nx, nz, 19, i, j, k, 2)];
        Double f3  = f[index_4d(nx, nz, 19, i, j, k, 3)];
        Double f4  = f[index_4d(nx, nz, 19, i, j, k, 4)];
        Double f5  = f[index_4d(nx, nz, 19, i, j, k, 5)];
        Double f6  = f[index_4d(nx, nz, 19, i, j, k, 6)];
        Double f7  = f[index_4d(nx, nz, 19, i, j, k, 7)];
        Double f8  = f[index_4d(nx, nz, 19, i, j, k, 8)];
        Double f9  = f[index_4d(nx, nz, 19, i, j, k, 9)];
        Double f10 = f[index_4d(nx, nz, 19, i, j, k, 10)];
        Double f11 = f[index_4d(nx, nz, 19, i, j, k, 11)];
        Double f12 = f[index_4d(nx, nz, 19, i, j, k, 12)];
        Double f13 = f[index_4d(nx, nz, 19, i, j, k, 13)];
        Double f14 = f[index_4d(nx, nz, 19, i, j, k, 14)];
        Double f15 = f[index_4d(nx, nz, 19, i, j, k, 15)];
        Double f16 = f[index_4d(nx, nz, 19, i, j, k, 16)];
        Double f17 = f[index_4d(nx, nz, 19, i, j, k, 17)];
        Double f18 = f[index_4d(nx, nz, 19, i, j, k, 18)];
        // compute density
        Double density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
        Double density_inv = 1. / density;

        // compute velocities
        Double utmp = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16)) * density_inv;
        Double vtmp = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18)) * density_inv;
        Double wtmp = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17)) * density_inv;

        // write to device
        rho[index_3d(nx, nz, i, j, k)] = density;
        u  [index_3d(nx, nz, i, j, k)] = utmp;
        v  [index_3d(nx, nz, i, j, k)] = vtmp;
        w  [index_3d(nx, nz, i, j, k)] = wtmp;
			}
		}
	}

	return ;
}

static void cpp_check_if_steady_state (double* const fA, 
		double* const u, double* const v, double* const w, double* const rho,
				const Params &params, int &converged, const int &D1, const int &D2, const int &D3) {

  const int nx = params.nx;
  const int nz = params.nz;

	const double tol = params.tol;

	converged = 1;
	for (int i = 1; i < D1 - 1; ++i) {
		for (int j = 1; j < D2 - 1; ++j) {
			for (int k = 1; k < D3 - 1; ++k) {
    		// load distributions A
        Double f0  = fA[index_4d(nx, nz, 19, i, j, k, 0)];
        Double f1  = fA[index_4d(nx, nz, 19, i, j, k, 1)];
        Double f2  = fA[index_4d(nx, nz, 19, i, j, k, 2)];
        Double f3  = fA[index_4d(nx, nz, 19, i, j, k, 3)];
        Double f4  = fA[index_4d(nx, nz, 19, i, j, k, 4)];
        Double f5  = fA[index_4d(nx, nz, 19, i, j, k, 5)];
        Double f6  = fA[index_4d(nx, nz, 19, i, j, k, 6)];
        Double f7  = fA[index_4d(nx, nz, 19, i, j, k, 7)];
        Double f8  = fA[index_4d(nx, nz, 19, i, j, k, 8)];
        Double f9  = fA[index_4d(nx, nz, 19, i, j, k, 9)];
        Double f10 = fA[index_4d(nx, nz, 19, i, j, k, 10)];
        Double f11 = fA[index_4d(nx, nz, 19, i, j, k, 11)];
        Double f12 = fA[index_4d(nx, nz, 19, i, j, k, 12)];
        Double f13 = fA[index_4d(nx, nz, 19, i, j, k, 13)];
        Double f14 = fA[index_4d(nx, nz, 19, i, j, k, 14)];
        Double f15 = fA[index_4d(nx, nz, 19, i, j, k, 15)];
        Double f16 = fA[index_4d(nx, nz, 19, i, j, k, 16)];
        Double f17 = fA[index_4d(nx, nz, 19, i, j, k, 17)];
        Double f18 = fA[index_4d(nx, nz, 19, i, j, k, 18)];

        // compute velocities A
        Double umom_A = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16));
        Double vmom_A = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18));
        Double wmom_A = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17));

        // load macroscopic properties B
        Double rho_B = rho[index_3d(nx, nz, i, j, k)];
        Double umom_B = u[index_3d(nx, nz, i, j, k)] * rho_B;
        Double vmom_B = v[index_3d(nx, nz, i, j, k)] * rho_B;
        Double wmom_B = w[index_3d(nx, nz, i, j, k)] * rho_B;

    		// convergence criteria
        int converged_momu = fabs(umom_B - umom_A) < std::max(tol * fabs(umom_A), 1e-12);
        int converged_momv = fabs(vmom_B - vmom_A) < std::max(tol * fabs(vmom_A), 1e-12);
        int converged_momw = fabs(wmom_B - wmom_A) < std::max(tol * fabs(wmom_A), 1e-12);
        converged &= converged_momu & converged_momv & converged_momw;
			}
		}
	}
	return ;
}