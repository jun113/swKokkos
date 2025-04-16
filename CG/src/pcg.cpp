#include "inc/functor.hpp"

#include "pcg.h"

#include <Kokkos_Core.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <crts.h>


// 示例
typedef struct
{
	double *p;
	double *z;
	double beta;
	int cells;
} Para;

// 矩阵向量乘法：Ap = A * p
typedef struct
{
	const CsrMatrix *csr_matrix; // 输入的CSR矩阵
	const double *vec;			 // 输入向量 p
	double *result;				 // 输出向量 Ap
	int rows;					 // 矩阵行数（计算任务的总大小）
} Para_spmv;

// 更新解，残差和搜索方向
typedef struct
{
	double *x;		 // 输出向量，更新后的结果（例如 x_{j+1}, r_{j+1}, 或 p_{j+1}）
	const double *y; // 输入向量1（例如 x_j, r_j, 或 r_{j+1}）
	const double *z; // 输入向量2（例如 p_j, μ, 或 p_j）
	double alpha;	 // 标量 alpha (aj 或 βj，正值或负值控制操作)
	int cells;		 // 向量长度
} Para_axpy;

// 点积或范数计算：result = (vec1, vec2) 或 ||vec1||^2
typedef struct
{
	const double *vec1; // 输入向量1
	const double *vec2; // 输入向量2（如果为 NULL，则计算 ||vec1||^2）
	double *result;		// 每个线程的部分点积结果存储数组
	int cells;			// 向量的总长度（计算任务的总大小）
} Para_reduce;

extern "C" void slave_example(Para *para);

extern "C" void slave_spmv(Para_spmv *para);	 // 矩阵向量乘法
extern "C" void slave_axpy(Para_axpy *para);	 // 通用 AXPY 算子
extern "C" void slave_reduce(Para_reduce *para); // 点积/范数计算

// 从核athread移植   朴素版

void cg_solve_athread(const CsrMatrix &csr_matrix, const double *rhs, double *psi, int max_iters, double tolerance, int *iter_out)
{

	// static int isInit = 0;
	// if (isInit == 0)
	// {
	// 	// 从核初始化
	// 	CRTS_init();
	// 	isInit = 1;
	// }

	int n = csr_matrix.rows;		 // 矩阵大小
	int nnz = csr_matrix.row_off[n]; // 非零元素个数

	// 分配临时向量
	double *r = (double *)malloc(n * sizeof(double));		  // 残差向量 r
	double *p = (double *)malloc(n * sizeof(double));		  // 搜索方向向量 p
	double *Ap = (double *)malloc(n * sizeof(double));		  // 矩阵向量乘积 Ap
	double *localsum = (double *)malloc(64 * sizeof(double)); // 用于从核返回部分和结果

	// 初始化残差 r = rhs - A * psi，psi 假设初始化为 0，则 r = rhs
	memcpy(r, rhs, n * sizeof(double));
	memcpy(p, r, n * sizeof(double)); // 初始搜索方向向量 p = r

	// 计算初始残差范数 ||r||^2
	double residual = 0.0, initial_residual = 0.0;

	Para_reduce reduce_para = {r, r, localsum, n}; // 用于计算 ||r||^2
	athread_spawn(slave_reduce, &reduce_para);	   // 从核计算 ||r||^2
	athread_join();
	for (int i = 0; i < 64; i++)
	{
		residual += localsum[i];
	}
	initial_residual = residual;

	// 如果初始残差已经满足容许误差，则直接返回
	if (sqrt(residual) / sqrt(initial_residual) < tolerance)
	{
		*iter_out = 0;
		free(r);
		free(p);
		free(Ap);
		free(localsum);
		return;
	}

	double alpha, beta, rdotr_new;

	// CG 迭代
	for (int iter = 0; iter < max_iters; iter++)
	{
		// Step 1: 计算 Ap = A * p
		Para_spmv spmv_para = {&csr_matrix, p, Ap, n};
		athread_spawn(slave_spmv, &spmv_para); // 从核执行矩阵向量乘法
		athread_join();

		// Step 2: 计算 alpha = (r, r) / (p, Ap)
		double pAp = 0.0;
		reduce_para.vec1 = p;
		reduce_para.vec2 = Ap;
		athread_spawn(slave_reduce, &reduce_para); // 从核计算点积 (p, Ap)
		athread_join();
		for (int i = 0; i < 64; i++)
		{
			pAp += localsum[i];
		}
		alpha = residual / pAp;

		// Step 3: 更新 psi 和 r
		// psi = psi + alpha * p
		Para_axpy axpy_para = {psi, psi, p, alpha, n};
		athread_spawn(slave_axpy, &axpy_para);
		athread_join();

		// r = r - alpha * Ap
		axpy_para.x = r;
		axpy_para.y = r;
		axpy_para.z = Ap;
		axpy_para.alpha = -alpha; // 注意 alpha 为负值
		athread_spawn(slave_axpy, &axpy_para);
		athread_join();

		// Step 4: 计算新的残差范数 ||r||^2
		rdotr_new = 0.0;
		reduce_para.vec1 = r;
		reduce_para.vec2 = r;
		athread_spawn(slave_reduce, &reduce_para); // 从核计算 ||r||^2
		athread_join();
		for (int i = 0; i < 64; i++)
		{
			rdotr_new += localsum[i];
		}

		// 检查相对残差是否满足收敛条件
		if (sqrt(rdotr_new) / sqrt(initial_residual) < tolerance)
		{
			*iter_out = iter + 1;
			free(r);
			free(p);
			free(Ap);
			free(localsum);
			return;
		}

		// Step 5: 计算 beta = ||r_new||^2 / ||r_old||^2
		beta = rdotr_new / residual;

		// Step 6: 更新搜索方向向量 p
		// p = r + beta * p
		axpy_para.x = p;
		axpy_para.y = r;
		axpy_para.z = p;
		axpy_para.alpha = beta;
		athread_spawn(slave_axpy, &axpy_para);
		athread_join();

		// 更新 residual 为新的 ||r_new||^2
		residual = rdotr_new;
	}

	// 达到最大迭代次数，返回最后的结果
	*iter_out = max_iters;

	// 释放内存
	free(r);
	free(p);
	free(Ap);
	free(localsum);
}

//  主核版
void cg_solve(const CsrMatrix &csr_matrix, const double *rhs, double *psi, int max_iters, double tolerance, int *iter_out)
{
	int n = csr_matrix.rows;		 // 矩阵大小（行数）
	int nnz = csr_matrix.row_off[n]; // 非零元素数量

	// 分配临时向量
	double *r = (double *)malloc(n * sizeof(double));  // 残差向量 r
	double *p = (double *)malloc(n * sizeof(double));  // 搜索方向向量 p
	double *Ap = (double *)malloc(n * sizeof(double)); // 矩阵向量乘积 Ap

	// Step 1: 初始化

	// 初始化残差向量 r = rhs
	memcpy(r, rhs, n * sizeof(double));
	memcpy(p, r, n * sizeof(double)); // 初始方向向量 p = r

	// for (int i = 0; i < 64; i++)
	// {

	// 	printf("r[%d]: %e\n", i, r[i]);
	// }

	double residual = 0.0;
	double initial_residual = 0.0;
	for (int i = 0; i < n; i++)
	{
		residual += r[i] * r[i];
	}
	initial_residual = residual;

	// 打印初始残差
	// printf("Initial residual: %e\n", sqrt(initial_residual));

	if (sqrt(residual) / sqrt(initial_residual) < tolerance)
	{
		*iter_out = 0;
		free(r);
		free(p);
		free(Ap);
		return;
	}

	double alpha, beta, rdotr_old, rdotr_new;

	// Step 2: 共轭梯度迭代
	for (int iter = 0; iter < max_iters; iter++)
	{
		// Ap = A * p
		memset(Ap, 0, n * sizeof(double));
		for (int i = 0; i < n; i++)
		{
			int row_start = csr_matrix.row_off[i];
			int row_end = csr_matrix.row_off[i + 1];
			for (int j = row_start; j < row_end; j++)
			{
				Ap[i] += csr_matrix.data[j] * p[csr_matrix.cols[j]];
			}
		}

		// 计算 alpha = (r, r) / (p, Ap)
		double pAp = 0.0;
		for (int i = 0; i < n; i++)
		{
			pAp += p[i] * Ap[i];
		}
		alpha = residual / pAp;

		// 更新 psi 和 r
		for (int i = 0; i < n; i++)
		{
			psi[i] += alpha * p[i]; // psi = psi + alpha * p
			r[i] -= alpha * Ap[i];	// r = r - alpha * Ap
		}

		// 检查残差是否收敛
		rdotr_new = 0.0;
		for (int i = 0; i < n; i++)
		{
			rdotr_new += r[i] * r[i];
		}
		if (sqrt(rdotr_new) / sqrt(initial_residual) < tolerance)
		{
			*iter_out = iter + 1;
			free(r);
			free(p);
			free(Ap);
			return;
		}

		// 计算 beta = (r_new, r_new) / (r_old, r_old)
		beta = rdotr_new / residual;

		// 更新 p = r + beta * p
		for (int i = 0; i < n; i++)
		{
			p[i] = r[i] + beta * p[i];
		}

		// 更新 residual 为新的 rdotr_new
		residual = rdotr_new;
	}

	// 如果达到最大迭代次数，返回最后的迭代结果
	*iter_out = max_iters;

	// 释放内存
	free(r);
	free(p);
	free(Ap);
}

void kokkos_cg_solve (const CsrMatrix &csr_matrix, 
    const double *rhs, double *psi, int max_iters, double tolerance, int *iter_out) {

	const int n   = csr_matrix.rows;		   // 矩阵大小（行数）
	const int nnz = csr_matrix.row_off[n]; // 非零元素数量

  ViewDouble1D v_matrx_data (csr_matrix.data, csr_matrix.data_size);
  ViewInt1D v_row_off       (csr_matrix.row_off, n + 1);
  ViewInt1D v_cols          (csr_matrix.cols, nnz);
  ViewDouble1D v_psi        (psi, n);

	// 分配临时向量
  Kokkos::View<double *> v_r  ("r",  n);
  Kokkos::View<double *> v_p  ("p",  n);
  Kokkos::View<double *> v_Ap ("Ap", n);

	// Step 1: 初始化

	// 初始化残差向量 r = rhs
	memcpy(v_r.data(), rhs, n * sizeof(double));
  // 初始方向向量 p = r
	memcpy(v_p.data(), rhs, n * sizeof(double)); 

	// for (int i = 0; i < 64; i++) {
	// 	printf ("r[%d]: %e\n", i, v_r(i));
	// }

	double residual = 0.0;
	double initial_residual = 0.0;

  Kokkos::parallel_reduce ("residual += v_r(i) * v_r(i)", n, MyReduce(v_r), residual);
	initial_residual = residual;

	// 打印初始残差
	// printf("Initial residual: %e\n", sqrt(initial_residual));

	if (sqrt(residual) / sqrt(initial_residual) < tolerance) {
		*iter_out = 0;
		return;
	}

  double alpha, beta, rdotr;

	// Step 2: 共轭梯度迭代
	for (int iter = 0; iter < max_iters; iter++) {

    Kokkos::parallel_for ("Ap = A * p", n,
        FunctorSpmv (v_row_off, v_cols, v_p, v_Ap, v_matrx_data));

		// 计算 alpha = (r, r) / (p, Ap)
		double pAp = 0.0;
    Kokkos::parallel_reduce ("pAp += v_p(i) * v_Ap(i)", n, MyReduce(v_p, v_Ap), pAp);
		alpha = residual / pAp;

		// 更新 psi 和 r
    Kokkos::parallel_for ("psi = alpha * p + psi", n,
        FunctorAXPY<double>(alpha, v_psi, v_p, v_psi));
    Kokkos::parallel_for ("r = - alpha * Ap + r", n,
        FunctorAXPY<double>(-alpha, v_r, v_Ap, v_r));

		// 检查残差是否收敛
		rdotr = 0.0;
    Kokkos::parallel_reduce ("rdotr += v_r(i) * v_r(i)", n, MyReduce(v_r), rdotr);

		if (sqrt(rdotr) / sqrt(initial_residual) < tolerance) {
			*iter_out = iter + 1;
			return;
		}

		// 计算 beta = (r_new, r_new) / (r_old, r_old)
		beta = rdotr / residual;

		// 更新 p = r + beta * p
    Kokkos::parallel_for ("p = beta * p + r", n,
        FunctorAXPY<double>(beta, v_p, v_p, v_r));

		// 更新 residual 为新的 rdotr_new
		residual = rdotr;
	}
  // 如果达到最大迭代次数，返回最后的迭代结果
  *iter_out = max_iters;
  return ;
}

// End
// ==============================
/*
void kokkos_cg_solve (const CsrMatrix &csr_matrix, 
    const double *rhs, double *psi, int max_iters, double tolerance, int *iter_out) {

	int n   = csr_matrix.rows;		   // 矩阵大小（行数）
	int nnz = csr_matrix.row_off[n]; // 非零元素数量

  Kokkos::View<double *> v_matrx_data ("CsrMatrix data", csr_matrix.data_size);
  Kokkos::View<int *> v_row_off ("csr matrix row off", n + 1);
  Kokkos::View<int *> v_cols ("csr matrix cols", nnz);

  Kokkos::View<double *> v_psi ("psi", n);

  Kokkos::View<double *>::HostMirror h_v_matrx_data = Kokkos::create_mirror_view(v_matrx_data);
  Kokkos::View<int *>::HostMirror    h_v_row_off    = Kokkos::create_mirror_view(v_row_off);
  Kokkos::View<int *>::HostMirror    h_v_cols       = Kokkos::create_mirror_view(v_cols);
  Kokkos::View<double *>::HostMirror h_v_psi        = Kokkos::create_mirror_view(v_psi);

	memcpy(h_v_matrx_data.data(), csr_matrix.data,    csr_matrix.data_size * sizeof(double));
	memcpy(h_v_row_off.data(),    csr_matrix.row_off, (n+1) * sizeof(int));
	memcpy(h_v_cols.data(),       csr_matrix.cols,    nnz   * sizeof(int));
	memcpy(h_v_psi.data(),        psi,                n * sizeof(double));

  Kokkos::deep_copy (v_matrx_data, h_v_matrx_data);
  Kokkos::deep_copy (v_row_off, h_v_row_off);
  Kokkos::deep_copy (v_cols, h_v_cols);
  Kokkos::deep_copy (v_psi, h_v_psi);

	// 分配临时向量
  Kokkos::View<double *> v_r  ("r",  n);
  Kokkos::View<double *> v_p  ("p",  n);
  Kokkos::View<double *> v_Ap ("Ap", n);

  Kokkos::View<double *>::HostMirror h_v_r  = Kokkos::create_mirror_view(v_r);
  Kokkos::View<double *>::HostMirror h_v_p  = Kokkos::create_mirror_view(v_p);

	// Step 1: 初始化

	// 初始化残差向量 r = rhs
	memcpy(h_v_r.data(), rhs, n * sizeof(double));
  // 初始方向向量 p = r
	memcpy(h_v_p.data(), rhs, n * sizeof(double)); 

  Kokkos::deep_copy (v_r, h_v_r);
  Kokkos::deep_copy (v_p, h_v_p);

	for (int i = 0; i < 64; i++) {
		printf ("r[%d]: %e\n", i, h_v_r(i));
	}

	double residual = 0.0;
	double initial_residual = 0.0;

  Kokkos::parallel_reduce ("residual += v_r(i) * v_r(i)", n,
      MyReduce<double>(v_r), residual);
	initial_residual = residual;

	// 打印初始残差
	printf("Initial residual: %e\n", sqrt(initial_residual));

	if (sqrt(residual) / sqrt(initial_residual) < tolerance) {
		*iter_out = 0;
		return;
	}
  // double alpha, beta, rdotr_old, rdotr_new;
  double alpha, beta, rdotr;

	// Step 2: 共轭梯度迭代
	for (int iter = 0; iter < max_iters; iter++) {

    Kokkos::parallel_for ("Ap = A * p", n,
        FunctorSpmv (v_row_off, v_cols, v_p, v_Ap, v_matrx_data));

		// 计算 alpha = (r, r) / (p, Ap)
		double pAp = 0.0;
    Kokkos::parallel_reduce ("pAp += v_p(i) * v_Ap(i)", n,
        MyReduce<double>(v_p, v_Ap), pAp);
		alpha = residual / pAp;

		// 更新 psi 和 r
    Kokkos::parallel_for ("psi = alpha * p + psi", n,
        FunctorAXPY<double>(alpha, v_psi, v_p, v_psi));
    Kokkos::parallel_for ("r = - alpha * Ap + r", n,
        FunctorAXPY<double>(-alpha, v_r, v_Ap, v_r));

		// 检查残差是否收敛
		rdotr = 0.0;
    Kokkos::parallel_reduce ("rdotr += v_r(i) * v_r(i)", n,
        MyReduce<double>(v_r), rdotr);

		if (sqrt(rdotr) / sqrt(initial_residual) < tolerance) {
			*iter_out = iter + 1;
      Kokkos::deep_copy (h_v_psi, v_psi);
	    memcpy(psi, h_v_psi.data(), n * sizeof(double));
			return;
		}

		// 计算 beta = (r_new, r_new) / (r_old, r_old)
		beta = rdotr / residual;

		// 更新 p = r + beta * p
    Kokkos::parallel_for ("p = beta * p + r", n,
        FunctorAXPY<double>(beta, v_p, v_p, v_r));

		// 更新 residual 为新的 rdotr_new
		residual = rdotr;
	}
  // 如果达到最大迭代次数，返回最后的迭代结果
  *iter_out = max_iters;
  Kokkos::deep_copy (h_v_psi, v_psi);
	memcpy(psi, h_v_psi.data(), n * sizeof(double));
  
  return ;
}
*/