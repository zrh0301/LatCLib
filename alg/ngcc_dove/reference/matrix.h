#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <stddef.h>
#include "gf256.h"

// 矩阵按行优先存储
typedef struct {
    size_t rows;
    size_t cols;
    uint8_t *data;
} gf256_matrix_t;

gf256_matrix_t* matrix_create(size_t rows, size_t cols);
void matrix_free(gf256_matrix_t *mat);

void matrix_add(gf256_matrix_t *C, const gf256_matrix_t *A, const gf256_matrix_t *B);
void matrix_mul(gf256_matrix_t *C, const gf256_matrix_t *A, const gf256_matrix_t *B);
void matrix_transpose(gf256_matrix_t *dst, const gf256_matrix_t *src);

// 对应算法中的 Upper 函数：下三角元素合并到上三角，下三角置0
void matrix_upper(gf256_matrix_t *mat);

// 矩阵求逆（高斯-约当消元），返回0成功，-1不可逆
int matrix_inv(gf256_matrix_t *dst, const gf256_matrix_t *src);

// 行向量 × 矩阵：res = v * A
void row_vec_mul_matrix(uint8_t *res, const uint8_t *v, const gf256_matrix_t *A);

// 矩阵 × 列向量：res = A * v
void matrix_mul_col_vec(uint8_t *res, const gf256_matrix_t *A, const uint8_t *v);

// 二次型：res = v^T * A * v
uint8_t quadratic_form(const uint8_t *v, const gf256_matrix_t *A);

// 用字节数组填充上三角矩阵（上三角元素按行优先排列）
void matrix_fill_upper(gf256_matrix_t *mat, const uint8_t *data);

// 导出上三角矩阵为字节数组
void matrix_dump_upper(const gf256_matrix_t *mat, uint8_t *data);

// 按行优先填充/导出完整矩阵
void matrix_fill(gf256_matrix_t *mat, const uint8_t *data);
void matrix_dump(const gf256_matrix_t *mat, uint8_t *data);


int matrix_solve(uint8_t *x, const gf256_matrix_t *A, const uint8_t *b);

#endif