#ifndef PQC_GF31_H
#define PQC_GF31_H

#include <stdint.h>

/* GF(31), the prime field used by the TSUOV multivariate signature.
 * Elements are stored as bytes in [0, 30]. */
typedef uint8_t gf31_t;

/* Reduce a non-negative integer. The fold uses 2^5 = 1 (mod 31) and is exact
 * for 0 <= z <= 1022, which covers every sum, difference and product of two
 * field elements as well as any input byte. */
gf31_t gf31_reduce(int z);
gf31_t gf31_add(gf31_t a, gf31_t b);
gf31_t gf31_sub(gf31_t a, gf31_t b);
gf31_t gf31_mul(gf31_t a, gf31_t b);
/* Returns 0 for zero, otherwise the multiplicative inverse. */
gf31_t gf31_inv(gf31_t a);

extern const gf31_t pqc_gf31_inv_table[31];

#endif
