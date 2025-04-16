#include "athread_param.h"

#include "slave.h"

__thread_local double reduce_buf_d  __attribute__ ((aligned(64)));
__thread_local float  reduce_buf_f  __attribute__ ((aligned(64)));

void athread_dot_product_2d_slave_double (struct para_dot_product_double* para) {
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;
	double* result = para->result;

	int N = I0 * I1;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	double reduce = 0.0;
	for (int i = start; i < end; ++i) {
		reduce += a[i] * b[i];
	}
	athread_ssync_array (); 
	athread_redurt (&reduce, &reduce, 1, athread_double, OP_add, &reduce_buf_d, 1);
	athread_ssync_array (); 
	if (athread_tid == 0) {
		*result = reduce;
	}
	return ;
}

void athread_dot_product_2d_slave_float (struct para_dot_product_float* para) {
	int I0   = para->I0;
	int I1   = para->I1;
	float* a = para->a;
	float* b = para->b;
	float* result = para->result;

	int N = I0 * I1;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	float reduce = 0.0;
	for (int i = start; i < end; ++i) {
		reduce += a[i] * b[i];
	}
	athread_ssync_array ();
	athread_redurt (&reduce, &reduce, 1, athread_float, OP_add, &reduce_buf_f, 1);
	if (athread_tid == 0) {
		*result = reduce;
	}
	// athread_ssync_array (); 
	return ;
}