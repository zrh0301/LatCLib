#ifndef PQC_POLY16_H
#define PQC_POLY16_H

#include <stddef.h>
#include <stdint.h>

#include "poly/fp16.h"

/* Coefficient-vector kernels over the 16-bit prime field. They are the int16_t
 * counterpart of the int32_t poly operators and keep the same lazy-reduction
 * contract: no kernel fully reduces unless it says so, and every range is
 * documented here so callers can chain operations safely. */

/* r = a + b coefficientwise. Exact in-place output is allowed. Inputs must
 * satisfy |a_i + b_i| <= INT16_MAX, which holds for the lazy ranges produced by
 * fp16_barrett_reduce and one level of addition. */
void fp16_poly_add(int16_t *r, const int16_t *a, const int16_t *b, size_t n);

/* r = a - b coefficientwise, same range contract as fp16_poly_add. */
void fp16_poly_sub(int16_t *r, const int16_t *a, const int16_t *b, size_t n);

/* Barrett reduction of every coefficient into (-q, q). Inputs must satisfy
 * |a_i| <= barrett*q, the bound of fp16_barrett_reduce. */
void fp16_poly_barrett(int16_t *r, size_t n, const fp16_ctx *ctx);

/* Add q where needed: fp16_poly_caddq maps (-q, q) into [0, q), and
 * fp16_poly_caddq2 maps (-2q, q) into [0, q). */
void fp16_poly_caddq(int16_t *r, size_t n, const fp16_ctx *ctx);
void fp16_poly_caddq2(int16_t *r, size_t n, const fp16_ctx *ctx);

/* Multiply every coefficient by R = 2^16 modulo q, i.e. move it into the
 * Montgomery domain. Each step is fp16_montgomery_reduce(a_i * mont_r2), so the
 * input may be any int16_t value. */
void fp16_poly_tomont(int16_t *r, size_t n, const fp16_ctx *ctx);

/* Centered binomial sampler over a byte buffer, as specified by the ICCS NGCC
 * reference samplers. Bits are read little-endian within each byte. eta is
 * split into a power-of-two high part hi (the largest power of two <= eta) and
 * lo = eta - hi: the first n*hi/4 bytes encode a canonical CBD_hi, and the
 * following lo bytes-per-coefficient block encodes a canonical CBD_lo, which is
 * added coefficientwise. For eta in {1, 2, 4} this is the plain contiguous CBD;
 * eta 3, 5 and 6 use the split layout of the reference implementations.
 * Buffer size must be eta*n/4 bytes, and eta must be at most 15. Output
 * coefficients lie in [-eta, eta]. */
void fp16_poly_cbd(int16_t *r, size_t n, const uint8_t *buf, unsigned eta);

/* Sample up to n coefficients by reading pairs of 12-bit little-endian
 * candidates from complete 3-byte groups. Returns the number accepted. */
size_t fp16_poly_rej_uniform12(int16_t *out, size_t n, const uint8_t *in,
                               size_t len, const fp16_ctx *ctx);

/* Little-endian coefficient bit strings, using 1 <= bits <= 15 and
 * ceil(n*bits/8) bytes. Encode retains the low bits of each coefficient.
 * Compress rounds and maps canonical coefficients in [0, q) to bits bits.
 * Decompress applies the corresponding rounded scaling back to modulus q. */
void fp16_poly_encode(uint8_t *out, const int16_t *in, size_t n, unsigned bits);
void fp16_poly_decode(int16_t *out, const uint8_t *in, size_t n, unsigned bits);
void fp16_poly_compress(uint8_t *out, const int16_t *in, size_t n,
                         unsigned bits, const fp16_ctx *ctx);
void fp16_poly_decompress(int16_t *out, const uint8_t *in, size_t n,
                           unsigned bits, const fp16_ctx *ctx);

/* Map little-endian message bits to 0 or ceil(q/2). Each chunk of
 * chunk_bytes message bytes is represented by consecutive copies of its
 * 8*chunk_bytes coefficients. copies is 2 or 4. msg_bytes must be a multiple
 * of nonzero chunk_bytes. The coefficient buffer has 8*copies*msg_bytes entries.
 * Decode accepts coefficients in [-q, q) and decides each bit by comparing
 * the sum of distances to floor(q/2) with floor(copies*q/4). */
void fp16_poly_encode_repeated_bits(int16_t *out, const uint8_t *msg,
                                    size_t msg_bytes, size_t chunk_bytes,
                                    unsigned copies, const fp16_ctx *ctx);
void fp16_poly_decode_repeated_bits(uint8_t *out, const int16_t *coeffs,
                                    size_t msg_bytes, size_t chunk_bytes,
                                    unsigned copies, const fp16_ctx *ctx);

#endif
