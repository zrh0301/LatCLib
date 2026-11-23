/// @file matrix.c
/// @brief DOVE 中使用的 GF(256) 矩阵基本运算。
///
/// 实现策略：
///   - 公开 API（matrix.h）保持不变，调用方代码无需修改。
///   - 内部实现直接基于 DOVE/2/ref 目录下 UOV 参考实现中的基本运算：
///       * 单元素原语：gf256_mul / gf256_inv（来自 uov/gf16.h）
///       * 向量原语  ：gf256v_add / gf256v_madd / gf256v_mul_scalar /
//                    gf256v_set_zero / gf256v_conditional_add
///       * 矩阵原语  ：gf256mat_prod_ref / gf256mat_gaussian_elim_ref /
//                    gf256mat_back_substitute_ref
///   - 与原实现使用相同的不可约多项式（0x11b），因此数学结果一致。
///   - 矩阵仍以行优先（row-major）方式存储；
///     在调用 UOV 列优先（column-major）原语时，通过参数映射转换语义。
///

#include "matrix.h"
#include "uov/blas_matrix_ref.h"   // UOV 的 gf256mat_*_ref 原语
#include <stdlib.h>
#include <string.h>

// ============================================================================
// 矩阵对象的创建 / 释放 / 数据导入导出
// ============================================================================

gf256_matrix_t* matrix_create(size_t rows, size_t cols) {
    gf256_matrix_t *mat = malloc(sizeof(gf256_matrix_t));
    if (!mat) return NULL;
    mat->rows = rows;
    mat->cols = cols;
    mat->data = calloc(rows * cols, sizeof(uint8_t));
    if (!mat->data) {
        free(mat);
        return NULL;
    }
    return mat;
}

void matrix_free(gf256_matrix_t *mat) {
    if (mat) {
        free(mat->data);
        free(mat);
    }
}

void matrix_fill(gf256_matrix_t *mat, const uint8_t *data) {
    memcpy(mat->data, data, mat->rows * mat->cols);
}

void matrix_dump(const gf256_matrix_t *mat, uint8_t *data) {
    memcpy(data, mat->data, mat->rows * mat->cols);
}

// 上三角矩阵（按行优先展开为 1 维数组）
void matrix_fill_upper(gf256_matrix_t *mat, const uint8_t *data) {
    size_t n = mat->rows;
    size_t idx = 0;
    gf256v_set_zero(mat->data, (unsigned)(n * n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            mat->data[i * n + j] = data[idx++];
        }
    }
}

void matrix_dump_upper(const gf256_matrix_t *mat, uint8_t *data) {
    size_t n = mat->rows;
    size_t idx = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            data[idx++] = mat->data[i * n + j];
        }
    }
}

// ============================================================================
// 矩阵加法
// ============================================================================

void matrix_add(gf256_matrix_t *C, const gf256_matrix_t *A, const gf256_matrix_t *B) {
    // C = A + B = A XOR B  （C 已有内容将被覆盖）
    unsigned total = (unsigned)(A->rows * A->cols);
    memcpy(C->data, A->data, total);          // C = A
    gf256v_add(C->data, B->data, total);      // C ^= B  => C = A XOR B
}

// ============================================================================
// 矩阵转置
// ============================================================================

void matrix_transpose(gf256_matrix_t *dst, const gf256_matrix_t *src) {
    for (size_t i = 0; i < src->rows; i++) {
        for (size_t j = 0; j < src->cols; j++) {
            dst->data[j * src->rows + i] = src->data[i * src->cols + j];
        }
    }
}

// ============================================================================
// Upper 变换：下三角元素合并到上三角，下三角置 0
// ============================================================================

void matrix_upper(gf256_matrix_t *mat) {
    size_t n = mat->rows;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            mat->data[i * n + j] ^= mat->data[j * n + i];
            mat->data[j * n + i] = 0;
        }
    }
}

// ============================================================================
// 行向量 × 矩阵：res = v * A
//
// 数学：res[j] = sum_i v[i] * A[i][j]
//
// 利用恒等式 v * A = A^T * v（作为列向量看待，元素相同），
// 将 A 的行优先存储视作 A^T 的列优先存储，直接调用 UOV 的 gf256mat_prod_ref：
//   gf256mat_prod_ref(c, matA, n_A_vec_byte, n_A_width, b) 计算 c = matA * b，
//   其中 matA 按列优先存储、每列 n_A_vec_byte 字节、共 n_A_width 列。
//   此处令 matA = A->data，n_A_vec_byte = A->cols（每"列"含 cols 字节 = 结果长度），
//   n_A_width = A->rows（列数 = v 长度），则得到 c = A^T * v = v * A。
// ============================================================================

void row_vec_mul_matrix(uint8_t *res, const uint8_t *v, const gf256_matrix_t *A) {
    gf256mat_prod_ref(res, A->data, (unsigned)A->cols, (unsigned)A->rows, v);
}

// ============================================================================
// 矩阵 × 列向量：res = A * v
//
// 数学：res[i] = sum_j A[i][j] * v[j]
//
// 由于 A 为行优先存储，UOV 的 gf256mat_prod_ref 需要列优先存储；
// 为避免每次都做转置拷贝，这里直接按行做内积：
//   res[i] = XOR_j (A[i][j] * v[j])
// 其中乘法和累加均使用 UOV 的基本原语（gf256_mul）。
// ============================================================================

void matrix_mul_col_vec(uint8_t *res, const gf256_matrix_t *A, const uint8_t *v) {
    size_t rows = A->rows, cols = A->cols;
    gf256v_set_zero(res, (unsigned)rows);
    for (size_t i = 0; i < rows; i++) {
        uint8_t acc = 0;
        for (size_t j = 0; j < cols; j++) {
            acc ^= gf256_mul(A->data[i * cols + j], v[j]);
        }
        res[i] = acc;
    }
}

// ============================================================================
// 矩阵乘法：C = A * B
//
// 实现：按行计算，对 A 的每一行 i，调用 row_vec_mul_matrix 计算 C[i] = A[i] * B。
// A 的第 i 行在 A->data 中连续存放（行优先），可直接作为 v 传入。
// 内部通过 UOV 的 gf256mat_prod_ref 完成行向量乘矩阵。
// ============================================================================

void matrix_mul(gf256_matrix_t *C, const gf256_matrix_t *A, const gf256_matrix_t *B) {
    size_t m = A->rows, k = A->cols, n = B->cols;
    for (size_t i = 0; i < m; i++) {
        row_vec_mul_matrix(&C->data[i * n], &A->data[i * k], B);
    }
}

// ============================================================================
// 二次型：res = v^T * A * v
// 步骤：tmp = A * v（使用 matrix_mul_col_vec），再 res = v^T * tmp（使用 gf256_dot）。
// ============================================================================

uint8_t quadratic_form(const uint8_t *v, const gf256_matrix_t *A) {
    size_t n = A->rows;
    uint8_t tmp[256];
    matrix_mul_col_vec(tmp, A, v);
    return gf256_dot(v, tmp, n);
}

// ============================================================================
// 线性方程组求解：求解 A * x = b
// 使用 UOV 的 gf256mat_gaussian_elim_ref + gf256mat_back_substitute_ref。
// 注：
//   - UOV 原语采用列优先存储且会修改矩阵；因此在内部缓冲中转一份。
//   - 原语要求 len <= MAX_H = 216；
//     DOVE 中的 DOVE_O 在 128/256/512 参数集下分别为 44/96/216，未超限。
// ============================================================================

int matrix_solve(uint8_t *x, const gf256_matrix_t *A, const uint8_t *b) {
    size_t n = A->rows;
    if (A->cols != n) return -1;        // 仅支持方阵

    // 拷贝 b 与 A 到工作缓冲（按列优先布局 A：mat[n * n]）
    uint8_t *sqmat = (uint8_t *)malloc(n * n);
    uint8_t *constant = (uint8_t *)malloc(n);
    if (!sqmat || !constant) {
        free(sqmat); free(constant);
        return -1;
    }
    // 列优先转置：A[i][j] -> sqmat[j*n + i]
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            sqmat[j * n + i] = A->data[i * n + j];
        }
        constant[i] = b[i];
    }

    // 调用 UOV 高斯消元
    unsigned succ = gf256mat_gaussian_elim_ref(sqmat, constant, (unsigned)n);
    if (!succ) {
        free(sqmat); free(constant);
        return -1;
    }
    // 回代
    gf256mat_back_substitute_ref(constant, sqmat, (unsigned)n);

    memcpy(x, constant, n);
    free(sqmat);
    free(constant);
    return 0;
}

// ============================================================================
// 矩阵求逆：使用高斯-约当消元法 [A | I] -> [I | A^-1]
// 内部使用 UOV 的单元素原语（gf256_mul / gf256_inv）。
// 注：DOVE 当前未调用本函数，但保留以保持 API 兼容。
// ============================================================================

int matrix_inv(gf256_matrix_t *dst, const gf256_matrix_t *src) {
    size_t n = src->rows;
    if (src->cols != n) return -1;

    // 构造 [A | I]（列优先：aug 行为 2n 列，每列 2n 字节）
    uint8_t *aug = (uint8_t *)calloc(n * 2 * n, sizeof(uint8_t));
    if (!aug) return -1;

    for (size_t i = 0; i < n; i++) {
        memcpy(&aug[i * 2 * n],     &src->data[i * n], n);
        aug[i * 2 * n + n + i] = 1;
    }

    // 高斯-约当消元
    for (size_t col = 0; col < n; col++) {
        // 找主元
        size_t pivot = col;
        for (; pivot < n; pivot++) {
            if (aug[pivot * 2 * n + col] != 0) break;
        }
        if (pivot == n) {
            free(aug);
            return -1;
        }

        // 交换行
        if (pivot != col) {
            for (size_t j = col; j < 2 * n; j++) {
                uint8_t tmp = aug[col * 2 * n + j];
                aug[col * 2 * n + j] = aug[pivot * 2 * n + j];
                aug[pivot * 2 * n + j] = tmp;
            }
        }

        // 主元行归一化
        uint8_t inv_p = gf256_inv(aug[col * 2 * n + col]);
        for (size_t j = col; j < 2 * n; j++) {
            aug[col * 2 * n + j] = gf256_mul(aug[col * 2 * n + j], inv_p);
        }

        // 消去其余行
        for (size_t row = 0; row < n; row++) {
            if (row == col) continue;
            uint8_t factor = aug[row * 2 * n + col];
            if (factor == 0) continue;
            for (size_t j = col; j < 2 * n; j++) {
                aug[row * 2 * n + j] ^= gf256_mul(factor, aug[col * 2 * n + j]);
            }
        }
    }

    // 提取逆矩阵
    for (size_t i = 0; i < n; i++) {
        memcpy(&dst->data[i * n], &aug[i * 2 * n + n], n);
    }
    free(aug);
    return 0;
}
