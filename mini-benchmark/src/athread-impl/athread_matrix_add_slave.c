#include "athread_param.h"

#include "slave.h"

void athread_matrix_add_slave_double (struct para_matrix_add_double* para) {
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;
	double* c = para->c;

	int N = I0 * I1;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	for (int i = start; i < end; ++i) {
		c[i] = a[i] + b[i];
	}
	return ;

}

void athread_matrix_add_slave_float (struct para_matrix_add_float* para) {
	int I0   = para->I0;
	int I1   = para->I1;
	float* a = para->a;
	float* b = para->b;
	float* c = para->c;

	int N = I0 * I1;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	for (int i = start; i < end; ++i) {
		c[i] = a[i] + b[i];
	}
	return ;

}