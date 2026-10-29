#ifndef PQC_GF16_MATRIX_H
#define PQC_GF16_MATRIX_H

#include <stddef.h>
#include "gf16.h"

/* Row-major matrices. Output storage must not overlap either input. Return
 * 0 on success or -1 for invalid pointers/dimensions; failure leaves output
 * unchanged. */
int gf16_matrix_mul(gf16_t *out, const gf16_t *a, const gf16_t *b,
                    size_t rows, size_t inner, size_t cols);
/* Accumulate out ^= a*b over GF(16). */
int gf16_matrix_mul_add(gf16_t *out, const gf16_t *a, const gf16_t *b,
                        size_t rows, size_t inner, size_t cols);

/* Determinant by Leibniz expansion. Supports dimensions 1 through 5, with
 * matrix-dependent memory access independent of the matrix values. */
int gf16_matrix_determinant(gf16_t *out, const gf16_t *matrix, size_t n);

#endif
