#include "gf16_matrix.h"

#include <stdint.h>

static int valid_mul(const gf16_t *out, const gf16_t *a, const gf16_t *b,
                    size_t rows, size_t inner, size_t cols) {
    return out != NULL && a != NULL && b != NULL && rows != 0 && inner != 0 && cols != 0 &&
           rows <= SIZE_MAX / cols && rows <= SIZE_MAX / inner && inner <= SIZE_MAX / cols;
}

int gf16_matrix_mul(gf16_t *out, const gf16_t *a, const gf16_t *b,
                    size_t rows, size_t inner, size_t cols) {
    if (!valid_mul(out, a, b, rows, inner, cols)) {
        return -1;
    }

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            gf16_t sum = 0;
            for (size_t k = 0; k < inner; ++k) {
                sum ^= gf16_mul(a[i * inner + k], b[k * cols + j]);
            }
            out[i * cols + j] = sum;
        }
    }
    return 0;
}

int gf16_matrix_mul_add(gf16_t *out, const gf16_t *a, const gf16_t *b,
                        size_t rows, size_t inner, size_t cols) {
    if (!valid_mul(out, a, b, rows, inner, cols)) {
        return -1;
    }

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            gf16_t sum = 0;
            for (size_t k = 0; k < inner; ++k) {
                sum ^= gf16_mul(a[i * inner + k], b[k * cols + j]);
            }
            out[i * cols + j] ^= sum;
        }
    }
    return 0;
}

static gf16_t determinant_terms(const gf16_t *matrix, size_t n, size_t row,
                                unsigned int used, gf16_t product) {
    if (row == n) {
        return product;
    }

    gf16_t det = 0;
    for (size_t col = 0; col < n; ++col) {
        unsigned int bit = 1u << col;
        if ((used & bit) == 0) {
            det ^= determinant_terms(matrix, n, row + 1, used | bit,
                                     gf16_mul(product, matrix[row * n + col]));
        }
    }
    return det;
}

int gf16_matrix_determinant(gf16_t *out, const gf16_t *matrix, size_t n) {
    if (out == NULL || matrix == NULL || n == 0 || n > 5) {
        return -1;
    }
    *out = determinant_terms(matrix, n, 0, 0, 1);
    return 0;
}
