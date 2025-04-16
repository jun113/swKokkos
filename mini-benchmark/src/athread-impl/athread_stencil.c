#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN (athread_5_point_stencil_slave_double) (struct para_5_point_stencil_double*);
extern void SLAVE_FUN (athread_5_point_stencil_slave_float) (struct para_5_point_stencil_float*);

void athread_5_point_stencil_double (int N, int I0, int I1, double* a, double* b) {
	struct para_5_point_stencil_double param;
	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	athread_spawn (athread_5_point_stencil_slave_double, &param);
	athread_join ();
	return ;
}

void athread_5_point_stencil_float (int N, int I0, int I1, float* a, float* b) {
	struct para_5_point_stencil_float param;

	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;

	athread_spawn (athread_5_point_stencil_slave_float, &param);
	athread_join ();
	return ;
}

extern void SLAVE_FUN (athread_5_point_stencil_slave_double_1) (struct para_5_point_stencil_double*);
extern void SLAVE_FUN (athread_5_point_stencil_slave_double_2) (struct para_5_point_stencil_double*);
extern void SLAVE_FUN (athread_5_point_stencil_slave_double_3) (struct para_5_point_stencil_double*);
extern void SLAVE_FUN (athread_5_point_stencil_slave_double_4) (struct para_5_point_stencil_double*);
extern void SLAVE_FUN (athread_5_point_stencil_slave_double_5) (struct para_5_point_stencil_double*);

void athread_5_point_stencil_double_1 (int N, int I0, int I1, double* a, double* b) {
	struct para_5_point_stencil_double param;
	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	athread_spawn (athread_5_point_stencil_slave_double_1, &param);
	athread_join ();
	return ;
}
void athread_5_point_stencil_double_2 (int N, int I0, int I1, double* a, double* b) {
	struct para_5_point_stencil_double param;
	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	athread_spawn (athread_5_point_stencil_slave_double_2, &param);
	athread_join ();
	return ;
}
void athread_5_point_stencil_double_3 (int N, int I0, int I1, double* a, double* b) {
	struct para_5_point_stencil_double param;
	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	athread_spawn (athread_5_point_stencil_slave_double_3, &param);
	athread_join ();
	return ;
}
void athread_5_point_stencil_double_4 (int N, int I0, int I1, double* a, double* b) {
	struct para_5_point_stencil_double param;
	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	athread_spawn (athread_5_point_stencil_slave_double_4, &param);
	athread_join ();
	return ;
}
void athread_5_point_stencil_double_5 (int N, int I0, int I1, double* a, double* b) {
	struct para_5_point_stencil_double param;
	param.N  = N;
	param.I0 = I0;
	param.I1 = I1;
	param.a  = a;
	param.b  = b;
	athread_spawn (athread_5_point_stencil_slave_double_5, &param);
	athread_join ();
	return ;
}