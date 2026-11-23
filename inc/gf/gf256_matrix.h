#ifndef PQC_GF256_MATRIX_H
#define PQC_GF256_MATRIX_H
#include "gf/gf256.h"
/* Row-major; nonzero dimensions, output does not overlap inputs.
 * Return -1 for NULL/zero/size overflow, output unchanged on failure. */
int gf256_matrix_mul(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t rows,size_t inner,size_t cols);
int gf256_matrix_mul_add(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t rows,size_t inner,size_t cols);

/* Row-major square matrix. Return -1 for invalid dimension, singular matrix
 * or allocation failure. Output is unchanged on failure. Exact in-place
 * solves (out == rhs) and inversions (out == matrix) are supported. */
int gf256_matrix_solve(uint8_t *out, const uint8_t *matrix, const uint8_t *rhs,
                       size_t n);
int gf256_matrix_inv(uint8_t *out, const uint8_t *matrix, size_t n);

/* Row-major matrix-vector multiplication with disjoint output. */
void gf256_matrix_mul_vec(uint8_t *out, const uint8_t *matrix, const uint8_t *v,
                          size_t rows, size_t cols);
/* v^T A v for a row-major n-by-n matrix. */
uint8_t gf256_matrix_quadratic(const uint8_t *v, const uint8_t *matrix, size_t n);
#endif
