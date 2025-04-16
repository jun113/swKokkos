#ifndef AXPY_SRC_HEAD_ATHREAD_PARAM_H_
#define AXPY_SRC_HEAD_ATHREAD_PARAM_H_
struct para_axpy_double {
	int N;
	double alpha;
	double* x;
	double* y;
};
struct para_axpy_float {
	int N;
	float alpha;
	float* x;
	float* y;
};

struct para_matrix_add_double {
	int I0;
	int I1;
	double *a;
	double *b;
	double *c;
};

struct para_matrix_add_float {
	int I0;
	int I1;
	float *a;
	float *b;
	float *c;
};

struct para_3d_mult_double {
	int I0;
	int I1;
	int I2;
	double *a;
	double *b;
	double *c;
};
struct para_3d_mult_float {
	int I0;
	int I1;
	int I2;
	float *a;
	float *b;
	float *c;
};

struct para_dot_product_double {
	int I0;
	int I1;
	double *a;
	double *b;
	double *result;
};

struct para_dot_product_float {
	int I0;
	int I1;
	float *a;
	float *b;
	float *result;
};

struct para_5_point_stencil_double {
	int N;
	int I0;
	int I1;
	double *a;
	double *b;
};

struct para_5_point_stencil_float {
	int N;
	int I0;
	int I1;
	float *a;
	float *b;
};

#endif // AXPY_SRC_HEAD_ATHREAD_PARAM_H_