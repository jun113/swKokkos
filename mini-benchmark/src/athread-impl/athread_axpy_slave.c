#include "athread_param.h"

#include "slave.h"

void athread_axpy_slave_double (struct para_axpy_double* para) {
	int N = para->N;
	double alpha = para->alpha;
	double* x    = para->x;
	double* y    = para->y;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	for (int i = start; i < end; ++i) {
		y[i] = alpha * x[i] + y[i];
	}
	return;
}

void athread_axpy_slave_float (struct para_axpy_float* para) {
	int N = para->N;
	float alpha = para->alpha;
	float* x    = para->x;
	float* y    = para->y;

	int num = (N + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < N ? (start + num) : N;

	for (int i = start; i < end; ++i) {
		y[i] = alpha * x[i] + y[i];
	}
	return;
}