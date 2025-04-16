#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN (athread_3d_mult_slave_double) (struct para_3d_mult_double*);

void athread_3d_mult_double (int I0, int I1, int I2, double* a, double* b, double* c) {
	struct para_3d_mult_double param;
	param.I0 = I0;
	param.I1 = I1;
	param.I2 = I2;
	param.a  = a;
	param.b  = b;
	param.c  = c;
	athread_spawn (athread_3d_mult_slave_double, &param);
	athread_join ();
	return ;
}

extern void SLAVE_FUN (athread_3d_mult_slave_float) (struct para_3d_mult_float*);
void athread_3d_mult_float (int I0, int I1, int I2, float* a, float* b, float* c) {
	struct para_3d_mult_float param;
	param.I0 = I0;
	param.I1 = I1;
	param.I2 = I2;
	param.a  = a;
	param.b  = b;
	param.c  = c;
	athread_spawn (athread_3d_mult_slave_float, &param);
	athread_join ();
	return ;
}