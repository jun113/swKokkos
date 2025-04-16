#include "athread_param.h"

#include "slave.h"

void athread_3d_mult_slave_double (struct para_3d_mult_double* para) {
	int I0    = para->I0;
	int I1    = para->I1;
	int I2    = para->I2;
	double* a = para->a;
	double* b = para->b;
	double* c = para->c;

	int N = I0 * I1 * I2;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	for (int i = start; i < end; ++i) {
		c[i] = a[i] * b[i];
	}
	return ;

}

void athread_3d_mult_slave_float (struct para_3d_mult_float* para) {
	int I0   = para->I0;
	int I1   = para->I1;
	int I2   = para->I2;
	float* a = para->a;
	float* b = para->b;
	float* c = para->c;

	int N = I0 * I1 * I2;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	for (int i = start; i < end; ++i) {
		c[i] = a[i] * b[i];
	}
	return ;

}