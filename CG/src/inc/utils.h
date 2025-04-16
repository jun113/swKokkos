#ifndef _UTILS_
#define _UTILS_


#include <math.h>
#include <stdio.h>

// 验证求解结果是否正确
void verify_solution(const double *psi, int n, double tolerance) {
    for (int i = 0; i < n; i++) {
        double diff = fabs(psi[i] - 1.0); // 计算 psi[i] 和 1.0 的绝对差值
        if (diff > tolerance) {
            // 如果任意位置的差值大于容许值，报错并返回
            printf("Verification failed at position %d: psi[%d] = %f, expected = 1.0, diff = %f\n", i, i, psi[i], diff);
            return;
        }
    }
    // 如果所有位置都在容许范围内，输出验证通过
    printf("Verification passed: all values within tolerance %.10e\n", tolerance);
}

// 写入CSV的函数
void write_csv(const char *full_path, double time_cost, int nnz, int iter, double gflops, const char *csv_filename)
{
    // 提取文件名部分
    const char *matrix_name = strrchr(full_path, '/');
    if (matrix_name)
        matrix_name++; // 跳过 '/'
    else
        matrix_name = full_path; // 如果没有路径，则直接使用原始字符串

    // 去掉文件扩展名（.mtx）
    char clean_name[256];
    strncpy(clean_name, matrix_name, sizeof(clean_name) - 1);
    clean_name[sizeof(clean_name) - 1] = '\0'; // 确保字符串以 NULL 结尾
    char *dot = strrchr(clean_name, '.');
    if (dot && strcmp(dot, ".mtx") == 0)
    {
        *dot = '\0'; // 删除 .mtx 扩展名
    }

    // 打开文件以追加模式写入
    FILE *file = fopen(csv_filename, "a");
    if (file == NULL)
    {
        printf("Error: Unable to open the CSV file.\n");
        return;
    }

    // 定义缓冲区并格式化数据到字符串
    char buffer[256]; // 假设一行不会超过 256 个字符
    sprintf(buffer, "%s,%.6f,%d,%d,%.6f\n", clean_name, time_cost, nnz, iter, gflops);

    // 将缓冲区内容写入文件
    fwrite(buffer, strlen(buffer), 1, file);

    // 关闭文件
    fclose(file);
}

#endif
