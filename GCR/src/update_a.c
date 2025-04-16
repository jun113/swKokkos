#include <omp.h>
//#include <iostream>
#include<string.h> 
//#include "inc/gcr.h"

#define index_a1(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*(kte-kts+3))
//#include "inc/gcr.h"
#define index_a(m, i,k,j) ((i-its)  + (k-kts+1) * (ite-its+1) + (j-jts) * (ite-its+1) * (kte-kts+3) +(m-1)* (jte-jts+1) * (ite-its+1) * (kte-kts+3))
 extern  void update_a(float * __restrict__ a_helm, float * __restrict__ a_helm_new,   int its, int ite, int jts, int jte, int kts, int kte) {
    // 定义常量和变量
    const int nx = ite - its + 1, nz = jte - jts + 1, ny = kte - kts + 3;
    const int size = 19*nx * ny * nz;
    int i, j, k, t;

    // 分配内存
//float *a_helm_new = (float*) malloc(size * sizeof(float));
//a_helm_new = (float*) malloc(size * sizeof(float));

    // 使用OpenMP并行更新数据
    #pragma omp parallel for private(i, j, k, t)   num_threads(4)  schedule(runtime) //schedule(runtime)  
    for (j = jts; j <= jte; j++) {
        for (k = kts-1; k <= kte+1; k++) {
            for (t= 1; t <= 19; t++) {
                for (i = its; i <= ite; i++) {
                    a_helm_new[index_a(t,i,k,j)] = a_helm[index_a1(t,i,k,j)];
                }
            }
        }
    }
/*
    #pragma omp parallel for private(i, j, k, t)   num_threads(4)  schedule(dynamic) //schedule(runtime)  
    for (j = jts; j <= jte; j++) {
        for (k = kts-1; k <= kte+1; k++) {
            for (t= 1; t <= 19; t++) {
                for (i = its; i <= ite; i++) {
                    a_helm[index_a(t,i,k,j)] = a_helm_new[index_a(t,i,k,j)];
                }
            }
        }
    }
*/

    // 将更新后的数据拷贝回原数组
//    memcpy(a_helm, a_helm_new, size * sizeof(float));
    // 释放内存
//    free(a_helm_new);
}


