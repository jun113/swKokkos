#include "athread_param.h"

#include "slave.h"

void athread_5_point_stencil_slave_double (struct para_5_point_stencil_double* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;

// {
// 	// Code 1
// 	int range = N * N;
// 	int num = (range + 63) / 64;
// 	int start = athread_tid * num;
// 	int end = (start + num) < range ? (start + num) : range;
// 	for (int idx = start; idx < end; ++idx) {
// 		int i0 = 1 + idx / N;
// 		int i1 = 1 + idx % N;
// 		b[i0*I1 + i1] += a[(i0  )*I1 + (i1  )]
//                    + a[(i0-1)*I1 + (i1  )]
//                    + a[(i0+1)*I1 + (i1  )]
//                    + a[(i0  )*I1 + (i1+1)]
//                    + a[(i0  )*I1 + (i1-1)];
// 	}
// 	// End code 1
// }
// {
// 	// Code 2
// 	int range = N;
// 	int num = (range + 63) / 64;
// 	int start = athread_tid * num;
// 	int end = (start + num) < range ? (start + num) : range;
// 	for (int i0 = 1; i0 < N+1; ++i0) {
// 		for (int idx = start; idx < end; ++idx) {
// 			int i1 = idx + 1;
// 			b[i0*I1 + i1] +=  a[(i0  )*I1 + (i1  )]
//                     	+ a[(i0-1)*I1 + (i1  )]
//                    		+ a[(i0+1)*I1 + (i1  )]
//                    		+ a[(i0  )*I1 + (i1+1)]
//                    		+ a[(i0  )*I1 + (i1-1)];
// 		}
// 	}
// 	// End code 2
// }
// {
// 	// Code 3
// 	int range = N;
// 	int num = (range + 63) / 64;
// 	int start = athread_tid * num;
// 	int end = (start + num) < range ? (start + num) : range;
// 	for (int idx = start; idx < end; ++idx) {
// 		int i0 = idx + 1;
// 		for (int i1 = 1; i1 < N+1; ++i1) {
// 			b[i0*I1 + i1] +=  a[(i0  )*I1 + (i1  )]
//                     	+ a[(i0-1)*I1 + (i1  )]
//                    		+ a[(i0+1)*I1 + (i1  )]
//                    		+ a[(i0  )*I1 + (i1+1)]
//                    		+ a[(i0  )*I1 + (i1-1)];
// 		}
// 	}
// 	// End code 3
// }
{
	// Code 4
	for (int i0 = athread_tid + 1; i0 < N+1; i0 += 64) {
		for (int i1 = 1; i1 < N+1; ++i1) {
			b[i0*I1 + i1] +=  a[(i0  )*I1 + (i1  )]
                    	+ a[(i0-1)*I1 + (i1  )]
                   		+ a[(i0+1)*I1 + (i1  )]
                   		+ a[(i0  )*I1 + (i1+1)]
                   		+ a[(i0  )*I1 + (i1-1)];
		}
	}
	// End code 4
}
	return ;
}

void athread_5_point_stencil_slave_float (struct para_5_point_stencil_float* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	float* a = para->a;
	float* b = para->b;
	int range = N * N;
	int num = (range + 63) / 64;

	int start = athread_tid * num;
	int end = (start + num) < range ? (start + num) : range;

	for (int idx = start; idx < end; ++idx) {
		int i0 = 1 + idx / I1;
		int i1 = 1 + idx % I1;
		// int i0 = 1 + (idx >> 12);
		// int i1 = 1 + (idx & (I1 - 1));

		b[i0*I1 + i1] += a[(i0  )*I1 + (i1  )]
                   + a[(i0-1)*I1 + (i1  )]
                   + a[(i0+1)*I1 + (i1  )]
                   + a[(i0  )*I1 + (i1+1)]
                   + a[(i0  )*I1 + (i1-1)];
	}
	return ;
}

inline int fast_div (int a, int b) {
  int result;
  for (result = 0; a >= b; result += 1) {
    a -= b;
  }
  return result;
}

void athread_5_point_stencil_slave_double_1 (struct para_5_point_stencil_double* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;

	// Code 1
	int range = N * N;
	int num = (range + 63) / 64;
	int start = athread_tid * num;
	int end = (start + num) < range ? (start + num) : range;
	for (int idx = start; idx < end; ++idx) {
		int i0 = idx / N + 1;
		int i1 = idx % N + 1;

		b[i0*I1 + i1] += a[(i0  )*I1 + (i1  )]
                   + a[(i0-1)*I1 + (i1  )]
                   + a[(i0+1)*I1 + (i1  )]
                   + a[(i0  )*I1 + (i1+1)]
                   + a[(i0  )*I1 + (i1-1)];
	}
	// End code 1
	return ;
}

void athread_5_point_stencil_slave_double_2 (struct para_5_point_stencil_double* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;

	// Code 2
	int range = N;
	int num = (range + 63) / 64;
	int start = athread_tid * num;
	int end = (start + num) < range ? (start + num) : range;
	for (int i0 = 1; i0 < N+1; ++i0) {
		for (int idx = start; idx < end; ++idx) {
			int i1 = idx + 1;
			b[i0*I1 + i1] +=  a[(i0  )*I1 + (i1  )]
                    	+ a[(i0-1)*I1 + (i1  )]
                   		+ a[(i0+1)*I1 + (i1  )]
                   		+ a[(i0  )*I1 + (i1+1)]
                   		+ a[(i0  )*I1 + (i1-1)];
		}
	}
	// End code 2
	return ;
}

void athread_5_point_stencil_slave_double_3 (struct para_5_point_stencil_double* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;

	// Code 3
	int range = N;
	int num = (range + 63) / 64;
	int start = athread_tid * num;
	int end = (start + num) < range ? (start + num) : range;
	for (int idx = start; idx < end; ++idx) {
		int i0 = idx + 1;
		for (int i1 = 1; i1 < N+1; ++i1) {
			b[i0*I1 + i1] +=  a[(i0  )*I1 + (i1  )]
                    	+ a[(i0-1)*I1 + (i1  )]
                   		+ a[(i0+1)*I1 + (i1  )]
                   		+ a[(i0  )*I1 + (i1+1)]
                   		+ a[(i0  )*I1 + (i1-1)];
		}
	}
	// End code 3
	return ;
}

void athread_5_point_stencil_slave_double_4 (struct para_5_point_stencil_double* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;
	// Code 4
	for (int i0 = athread_tid + 1; i0 < N+1; i0 += 64) {
		for (int i1 = 1; i1 < N+1; ++i1) {
			b[i0*I1 + i1] +=  a[(i0  )*I1 + (i1  )]
                    	+ a[(i0-1)*I1 + (i1  )]
                   		+ a[(i0+1)*I1 + (i1  )]
                   		+ a[(i0  )*I1 + (i1+1)]
                   		+ a[(i0  )*I1 + (i1-1)];
		}
	}
	// End code 4
	return ;
}

void athread_5_point_stencil_slave_double_5 (struct para_5_point_stencil_double* para) {
	int N     = para->N;
	int I0    = para->I0;
	int I1    = para->I1;
	double* a = para->a;
	double* b = para->b;

	// Code 5
	int range = N * N;
	int num = (range + 63) / 64;
	int start = athread_tid * num;
	int end = (start + num) < range ? (start + num) : range;
	for (int idx = start; idx < end/2; ++idx) {
    int idx1 = idx  * 2;
    int idx2 = idx1 + 1;
		int i0 = 1 + idx1 / N;
		int i1 = 1 + idx1 % N;
		b[i0*I1 + i1] += a[(i0  )*I1 + (i1  )]
                   + a[(i0-1)*I1 + (i1  )]
                   + a[(i0+1)*I1 + (i1  )]
                   + a[(i0  )*I1 + (i1+1)]
                   + a[(i0  )*I1 + (i1-1)];
		i0 = 1 + idx2 / N;
		i1 = 1 + idx2 % N;
		b[i0*I1 + i1] += a[(i0  )*I1 + (i1  )]
                   + a[(i0-1)*I1 + (i1  )]
                   + a[(i0+1)*I1 + (i1  )]
                   + a[(i0  )*I1 + (i1+1)]
                   + a[(i0  )*I1 + (i1-1)];
	}
	// End code 5
	return ;
}