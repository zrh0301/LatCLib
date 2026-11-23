#ifndef PQC_GF16_H
#define PQC_GF16_H

#include <stdint.h>

/* GF(16) in the polynomial basis modulo x^4 + x + 1. Inputs and outputs
 * are nibble values in [0, 15]. Operations have fixed control flow. */
typedef uint8_t gf16_t;

gf16_t gf16_add(gf16_t a, gf16_t b);
gf16_t gf16_sub(gf16_t a, gf16_t b);
gf16_t gf16_mul(gf16_t a, gf16_t b);
/* Returns 0 for zero, otherwise the multiplicative inverse. */
gf16_t gf16_inv(gf16_t a);

#endif
