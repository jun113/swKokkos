#ifndef LBM_SRC_INCLUDE_ATHREAD_LBM_PARAMS_H_
#define LBM_SRC_INCLUDE_ATHREAD_LBM_PARAMS_H_

#include "params.h"

struct AthreadParams {
	double* fB;
	double* fA;
	double* u;
	double* v;
	double* w;
	double* rho;
	struct Params* params;
	int* converged;
	int D1;
	int D2;
	int D3;
};

#endif // LBM_SRC_INCLUDE_ATHREAD_LBM_PARAMS_H_