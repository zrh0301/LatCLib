#ifndef PQC_GF256_MATRIX_H
#define PQC_GF256_MATRIX_H
#include "gf256.h"
/* Row-major; nonzero dimensions, output does not overlap inputs.
 * Return -1 for NULL/zero/size overflow, output unchanged on failure. */
int gf256_matrix_mul(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t rows,size_t inner,size_t cols);
int gf256_matrix_mul_add(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t rows,size_t inner,size_t cols);
#endif
