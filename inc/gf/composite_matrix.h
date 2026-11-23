#ifndef PQC_GF_COMPOSITE_MATRIX_H
#define PQC_GF_COMPOSITE_MATRIX_H

#include "gf/composite.h"

typedef enum {
    GF_COMPOSITE_256,
    GF_COMPOSITE_16
} gf_composite_field;

/* Invert a nonzero-diagonal triangular matrix stored by row with its
 * diagonal. The input and output each contain n*(n+1)/2 unpacked entries. */
void gf_composite_triangular_inv(uint8_t *out, const uint8_t *in, size_t n,
                                 int upper, gf_composite_field field);

/* out = L*U, or U*L when upper_first is set. L and U use unpacked triangular
 * storage. Output is column-major, packed in the GF(16) case. */
void gf_composite_triangular_product(uint8_t *out, const uint8_t *lower,
                                      const uint8_t *upper, size_t n, int upper_first,
                                      gf_composite_field field);

/* Solve A*x=b with VDOO's pivot rule. If a diagonal entry is zero, combine all
 * later rows into that row. If the combined pivot is also zero, return -1.
 * Input matrix is row-major and all entries (including x and b) are unpacked.
 * workspace holds n*(n+1) bytes. On failure x is unchanged. */
int gf_composite_matrix_solve(uint8_t *x, const uint8_t *a, const uint8_t *b,
                              size_t n, uint8_t *workspace, gf_composite_field field);

/* Multiply a column-major matrix by a vector, with independent strides for
 * columns. GF(16) inputs and output are nibble-packed. */
void gf_composite_col_matvec(uint8_t *out, const uint8_t *a, const uint8_t *x,
                             size_t rows, size_t cols, size_t stride,
                             gf_composite_field field);

/* Batched entries are consecutive vectors of batch_bytes bytes. Matrices of
 * batched entries are row-major. Matrices of scalars are column-major, using
 * scalar_stride bytes per column and packed nibbles for GF(16).
 * tri_madd computes out += upper*scalars, with upper packed by row.
 * transpose_madd computes out += scalars^T * batch. */
void gf_composite_batch_tri_madd(uint8_t *out, const uint8_t *upper,
                                 const uint8_t *scalars, size_t n, size_t scalar_stride,
                                 size_t cols, size_t batch_bytes, gf_composite_field field);
void gf_composite_batch_transpose_madd(uint8_t *out, const uint8_t *scalars,
                                       size_t inner, size_t scalar_stride, size_t rows,
                                       const uint8_t *batch, size_t cols, size_t batch_bytes,
                                       gf_composite_field field);

/* Evaluate a batch of upper-triangular quadratic forms. x is public. Each
 * upper-triangular entry holds batch_bytes packed coefficients. Scratch
 * holds 16*batch_bytes or 256*batch_bytes bytes, depending on field. */
void gf_composite_batch_quad(uint8_t *out, const uint8_t *upper, const uint8_t *x,
                             size_t n, size_t batch_bytes, uint8_t *scratch,
                             gf_composite_field field);

#endif
