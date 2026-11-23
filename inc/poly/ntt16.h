#ifndef PQC_NTT16_H
#define PQC_NTT16_H

#include <stddef.h>
#include <stdint.h>

/* 128-point negacyclic NTT over q = 3329 for int16_t coefficient vectors.
 * The forward transform accepts |a_i| <= q and Barrett-reduces the output
 * into (-q, q). The inverse includes the 1/128 scaling and accepts the
 * same range, so invntt(ntt(a)) == a (mod q). */
void pqc_ntt128(int16_t *a);
void pqc_invntt128(int16_t *a);

/* Apply the transform to blocks consecutive independent 128-coefficient
 * vectors. */
void pqc_ntt128_blocks(int16_t *a, size_t blocks);
void pqc_invntt128_blocks(int16_t *a, size_t blocks);

/* Product in the split NTT domain, modulo x^(128*blocks) + 1. Each
 * lane is a polynomial in x reduced by x^blocks minus the lane's root.
 * Inputs are reduced NTT coefficients in [-3329, 3329]. Output has one
 * inverse Montgomery factor and is lazily reduced. Call fp16_poly_tomont
 * before taking the inverse transform. Exact aliasing is allowed.
 * Accepts 4, 8 or 16 blocks. Returns -1 for other sizes with output unchanged. */
int pqc_ntt128_blocks_mont_mul(int16_t *out, const int16_t *a, const int16_t *b,
                                size_t blocks);

#endif
