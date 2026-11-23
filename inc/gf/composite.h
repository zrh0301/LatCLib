#ifndef PQC_GF_COMPOSITE_H
#define PQC_GF_COMPOSITE_H

#include <stddef.h>
#include <stdint.h>

/* VDOO tower coordinates. GF(4) = GF(2)[x]/(x²+x+1),
 * GF(16) = GF(4)[y]/(y²+y+x), GF(256) = GF(16)[z]/(z²+z+xy).
 * Each extension stores the constant coefficient in the low bits.
 * GF(16) scalar inputs are nibbles. The multiplicative inverse of 0 is 0. */
uint8_t gf16_composite_mul(uint8_t a, uint8_t b);
uint8_t gf16_composite_inv(uint8_t a);
uint8_t gf256_composite_mul(uint8_t a, uint8_t b);
uint8_t gf256_composite_inv(uint8_t a);

/* Accumulate acc ^= scalar * in. len counts bytes. GF(16) has two nibbles
 * per byte, low nibble first. Exact in-place operation is supported. */
void gf16_composite_packed_madd(uint8_t *acc, const uint8_t *in, uint8_t scalar,
                                size_t len);
void gf256_composite_madd(uint8_t *acc, const uint8_t *in, uint8_t scalar,
                          size_t len);

#endif
