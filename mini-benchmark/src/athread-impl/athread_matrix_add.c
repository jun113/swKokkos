#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN (athread_matrix_add_slave_double) (struct para_matrix_add_double*);

void athread_matrix_add_double (int I0, int I1, double* a, double* b, double* c) {
	struct para_matrix_add_double param_matrix_add_double;

	param_matrix_add_double.I0 = I0;
	param_matrix_add_double.I1 = I1;
	param_matrix_add_double.a  = a;
	param_matrix_add_double.b  = b;
	param_matrix_add_double.c  = c;

	athread_spawn (athread_matrix_add_slave_double, &param_matrix_add_double);
	athread_join ();

	return ;
}

extern void SLAVE_FUN (athread_matrix_add_slave_float) (struct para_matrix_add_float*);
void athread_matrix_add_float (int I0, int I1, float* a, float* b, float* c) {
	struct para_matrix_add_float param;

	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	param.c  = c;

	athread_spawn (athread_matrix_add_slave_double, &param);
	athread_join ();

	return ;
}