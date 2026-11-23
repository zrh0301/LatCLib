#ifndef PQC_GF_BYTE_MATRIX_H
#define PQC_GF_BYTE_MATRIX_H

#include <stddef.h>
#include <stdint.h>

/* Each row-major entry occupies one byte. Exact aliasing is allowed for add.
 * The other operations require disjoint input and output except fold_upper,
 * which updates its square matrix in place. */
void gf_byte_add(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t len);
void gf_byte_transpose(uint8_t *out, const uint8_t *in, size_t rows, size_t cols);
void gf_byte_fold_upper(uint8_t *matrix, size_t n);
void gf_byte_pack_upper(uint8_t *out, const uint8_t *in, size_t n);
void gf_byte_unpack_upper(uint8_t *out, const uint8_t *in, size_t n);

/* Accumulate a square row-major batch of matrices into their packed upper
 * triangles, folding each strict lower entry into its transposed entry. */
void gf_byte_batch_fold_upper_add(uint8_t *out, const uint8_t *in, size_t n,
                                  size_t batch_bytes);

#endif
