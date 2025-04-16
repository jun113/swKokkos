#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN (athread_axpy_slave_double) (struct para_axpy_double*);

void athread_axpy_double (int N, double alpha, double* x, double* y) {

	struct para_axpy_double param_axpy_double;

	param_axpy_double.N     = N;
	param_axpy_double.alpha = alpha;
	param_axpy_double.x     = x;
	param_axpy_double.y     = y;

	athread_spawn (athread_axpy_slave_double, &param_axpy_double);
	athread_join ();
	return ;
}

extern void SLAVE_FUN (athread_axpy_slave_float) (struct para_axpy_float*);

void athread_axpy_float (int N, float alpha, float* x, float* y) {

	struct para_axpy_float param_axpy_float;

	param_axpy_float.N     = N;
	param_axpy_float.alpha = alpha;
	param_axpy_float.x     = x;
	param_axpy_float.y     = y;

	athread_spawn (athread_axpy_slave_float, &param_axpy_float);
	athread_join ();
	return ;
}