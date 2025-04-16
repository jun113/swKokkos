#include <math.h>
#include "pcg_def.h"
#include <string.h>
#include <crts.h>



typedef struct
{
	int rows;
	int *row_off;
	int *cols;
	double *data;
	int data_size;
} CsrMatrix;



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

// 矩阵向量乘法：Ap = A * p
void slave_spmv(Para_spmv *para)
{
	// 拉取参数到本地
	Para_spmv slavePara;
	CRTS_dma_get(&slavePara, para, sizeof(Para_spmv));

	const CsrMatrix *csr_matrix = slavePara.csr_matrix;
	const double *vec = slavePara.vec;
	double *result = slavePara.result;
	int rows = csr_matrix->rows;

	// 当前线程的任务划分
	int len = rows / 64;
	int rest = rows % 64;
	int start, end;
	if (CRTS_tid < rest)
	{
		len++;
		start = CRTS_tid * len;
	}
	else
	{
		start = CRTS_tid * len + rest;
	}
	end = start + len;

	// 矩阵向量乘法计算
	for (int i = start; i < end; i++)
	{
		result[i] = 0.0;
		int row_start = csr_matrix->row_off[i];
		int row_end = csr_matrix->row_off[i + 1];
		for (int j = row_start; j < row_end; j++)
		{
			result[i] += csr_matrix->data[j] * vec[csr_matrix->cols[j]];
		}
	}
}

// 更新解，残差和搜索方向
void slave_axpy(Para_axpy *para)
{
	// 拉取参数到本地
	Para_axpy slavePara;
	CRTS_dma_get(&slavePara, para, sizeof(Para_axpy));

	double *x = slavePara.x;
	const double *y = slavePara.y;
	const double *z = slavePara.z;
	double alpha = slavePara.alpha;
	int cells = slavePara.cells;

	// 当前线程任务划分
	int len = cells / 64;  // 每个线程分配的元素数
	int rest = cells % 64; // 剩余的元素数
	int start, end;
	if (CRTS_tid < rest)
	{
		len++;
		start = CRTS_tid * len;
	}
	else
	{
		start = CRTS_tid * len + rest;
	}
	end = start + len;

	// 执行通用 AXPY 计算
	for (int i = start; i < end; i++)
	{
		x[i] = y[i] + alpha * z[i]; // 通用加法逻辑
	}
}

// 点积或范数计算：result = (vec1, vec2) 或 ||vec1||^2
void slave_reduce(Para_reduce *para)
{
	// 拉取参数到本地
	Para_reduce slavePara;
	CRTS_dma_get(&slavePara, para, sizeof(Para_reduce));

	const double *vec1 = slavePara.vec1;
	const double *vec2 = slavePara.vec2;
	double *result = slavePara.result;
	int cells = slavePara.cells;

	// 当前线程的任务划分
	int len = cells / 64;
	int rest = cells % 64;
	int start, end;
	if (CRTS_tid < rest)
	{
		len++;
		start = CRTS_tid * len;
	}
	else
	{
		start = CRTS_tid * len + rest;
	}
	end = start + len;

	// 当前线程的部分和
	double local_sum = 0.0;
	for (int i = start; i < end; i++)
	{
		local_sum += vec1[i] * vec2[i]; // 点积
	}

	// 将结果写入每个线程对应的 result 位置
	result[CRTS_tid] = local_sum;
}
