#ifndef PQC_GF31_MATRIX_H
#define PQC_GF31_MATRIX_H

#include <stddef.h>
#include <stdint.h>

#include "gf/gf31.h"

/* GF(31) vector operators. One field element per byte. */
void gf31_vector_add(gf31_t *c, const gf31_t *d, int len);
gf31_t gf31_vector_dot(const gf31_t *a, const gf31_t *b, size_t len);
void gf31_vector_mul_scalar(const gf31_t *v, gf31_t scalar, size_t len, gf31_t *res);
/* Lazy helpers: the sums stay unreduced (<= 60) and are folded by
 * gf31_vector_reduce or gf31_vector_add_mod afterwards. */
void gf31_vector_add_no_mod(const gf31_t *a, const gf31_t *b, size_t len, gf31_t *res);
void gf31_vector_reduce(const gf31_t *v, size_t len, gf31_t *res);
void gf31_vector_add_mod(const gf31_t *a, const gf31_t *b, size_t len, gf31_t *res);
void gf31_vector_sub_mod(const gf31_t *a, const gf31_t *b, size_t len, gf31_t *res);

/* GF(31^2) = GF(31)[x]/(x^2 - 3), the degree-2 tower used by TSUOV.
 *
 * Planar layout: a vector of n pairs is stored as the n low coefficients
 * followed by the n high coefficients. Matrices are column-major with the
 * same split and with each column padded to a multiple of 32 bytes, which is
 * the layout produced by the reference expansion routines. */
void gf31_t2_matrix_mul(const uint8_t *a, const uint8_t *b, uint8_t *c,
                        int a_rows, int a_cols, int b_cols, int c_rows);
/* Row vector times a matrix. */
void gf31_t2_matrix_vec_mul(const uint8_t *v, const uint8_t *m, uint8_t *c,
                            int v_cols, int m_cols, int c_cols);
/* Coordinate map used by the TSUOV E operator on the planar layout:
 * (a0, a1) -> (a0 + a1, 8 * a1). */
void gf31_t2_phi_planar(uint8_t *data, uint8_t *dst, size_t n_pairs);

#endif
