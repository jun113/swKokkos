#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN (athread_dot_product_2d_slave_double) (struct para_dot_product_double*);

void athread_dot_product_2d_double (int I0, int I1, double* a, double* b, double* result) {
	struct para_dot_product_double param;
	param.I0     = I0;
	param.I1     = I1;
	param.a      = a;
	param.b      = b;
	param.result = result;
	athread_spawn (athread_dot_product_2d_slave_double, &param);
	athread_join ();
	return ;
}

extern void SLAVE_FUN (athread_dot_product_2d_slave_float) (struct para_dot_product_float*);

void athread_dot_product_2d_float (int I0, int I1, float* a, float* b, float* result) {
	struct para_dot_product_float param;
	param.I0     = I0;
	param.I1     = I1;
	param.a      = a;
	param.b      = b;
	param.result = result;
	athread_spawn (athread_dot_product_2d_slave_float, &param);
	athread_join ();
	return ;
}