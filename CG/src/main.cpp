#include <Kokkos_Core.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "./biio.h"
#include "pcg.h"
#include "utils.h"

#include <iomanip>
#include <fstream>

#define MAX_ITERS 10000

#define rtol 1.0e-8
#define tolerance 1.0e-4

int main(int argc, char *argv[]) {
	struct timeval tv;
	double total_time_cost = .0;

	// input the matrix

	INFO("Read matrix!\n");

	char *filename = argv[1];
	int m, n, nnzR, isSymmetric;
	int *RowPtr;
	int *ColIdx;
	double *Val;
	read_Dmatrix_32(&m, &n, &nnzR, &RowPtr, &ColIdx, &Val, &isSymmetric, filename);
	if (m != n) {
		printf("unequal\n");
		return 0;
	}

	INFO("%d,%d,%d\n\n", m, n, nnzR);

  std::ofstream file("cg_results.csv", std::ios::app);
  if (!file) {
    std::cerr << "err in "<< filename << std::endl;
    return 1;
  }
	double t_start (0.0), t_end (0.0);

	// construct the csrmatrix
	CsrMatrix csr_matrix;
	csr_matrix.rows = m;
	csr_matrix.row_off = RowPtr;
	csr_matrix.cols = ColIdx;
	csr_matrix.data = Val;
	csr_matrix.data_size = nnzR;

	// construct the rhs

	double *X = (double *)malloc(sizeof(double) * (n));
	double *rhs = (double *)malloc(sizeof(double) * (m));
	double *Y_golden = (double *)malloc(sizeof(double) * (m));
	memset(X, 0, sizeof(double) * (n));
	memset(rhs, 0, sizeof(double) * (n));
	memset(Y_golden, 0, sizeof(double) * (n));

	for (int i = 0; i < n; i++) {
		X[i] = 1;
	}

	for (int i = 0; i < n; i++)
		for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++)
			Y_golden[i] += Val[j] * X[ColIdx[j]];

	// The righthand side vector is set to be the product of the input matrix and an all 1.0 vector,
	// and the initial value of the solution vector is set to zero
	memcpy(rhs, Y_golden, sizeof(double) * n);
	memset(X, 0, sizeof(double) * (n));

	int iter_cpp (0), iter_ath (0), iter_kok (0);
	// solve stage

	INFO("swcg solve start!\n\n");

	gettimeofday(&tv, NULL);
	t_start = (double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6;
	cg_solve(csr_matrix, rhs, X, MAX_ITERS, rtol, &iter_cpp);
	gettimeofday(&tv, NULL);
	t_end = (double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6;
	double t_cpp = t_end - t_start;

	memset(X, 0, sizeof(double) * (n));
  athread_init();

	gettimeofday(&tv, NULL);
	t_start = (double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6;
	cg_solve_athread(csr_matrix, rhs, X, MAX_ITERS, rtol, &iter_ath);
	gettimeofday(&tv, NULL);
	t_end = (double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6;
	double t_ath = t_end - t_start;
  athread_halt();

	memset(X, 0, sizeof(double) * (n));
  Kokkos::initialize();

	gettimeofday(&tv, NULL);
	t_start = (double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6;
	kokkos_cg_solve(csr_matrix, rhs, X, MAX_ITERS, rtol, &iter_kok);
	gettimeofday(&tv, NULL);
	t_end = (double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6;
	double t_kok = t_end - t_start;

	// 计算Gflops
  double gflops_cpp = (2.0 * nnzR + 5.0 * 2.0 * n) * iter_cpp / (t_cpp * 1e9);
  double gflops_ath = (2.0 * nnzR + 5.0 * 2.0 * n) * iter_ath / (t_ath * 1e9);
  double gflops_kok = (2.0 * nnzR + 5.0 * 2.0 * n) * iter_kok / (t_kok * 1e9);

  // 将结果写入CSV文件
  // write_csv(filename, total_time_cost, nnzR, iter, gflops, "cg_kokkos_results.csv");

  file << filename << "," << nnzR << "," << t_cpp << "," << gflops_cpp << "," << iter_cpp
                                  << "," << t_ath << "," << gflops_ath << "," << iter_ath
 				  << "," << t_kok << "," << gflops_kok << "," << iter_kok << std::endl;

	file.close();

	// INFO("swcg solve time: %.4lfs\n\n", total_time_cost);

	// 输出结果
	// INFO("Iterations: %d\n", iter);
	// INFO("Solution vector:\n");
	// for (int i = 0; i < 10; i++)
	// {
	// 	printf("%f\n", X[i]);
	// }

	// 调用验证函数
	// verify_solution(X, m, tolerance);

  Kokkos::finalize();
	return 0;
}
