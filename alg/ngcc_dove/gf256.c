/* DOVE GF(256) vector layer backed by the library kernels.
 *
 * DOVE's reference gf256.c implemented these helpers on top of the UOV
 * 32-bit word primitives. The field is the same (0x11b), so the whole layer
 * is expressed with the public gf/gf256.h API instead. Bit-exact results are
 * covered by the KAT test. */

#include "gf256.h"
#include <string.h>

uint8_t gf256_add(uint8_t a, uint8_t b)
{
    return (uint8_t)(a ^ b);
}

uint8_t gf256_dot(const uint8_t *a, const uint8_t *b, size_t len)
{
    uint8_t sum = 0;
    for (size_t i = 0; i < len; i++)
        sum ^= gf256_mul(a[i], b[i]);
    return sum;
}

void gf256v_add(uint8_t *b, const uint8_t *a, unsigned num_byte)
{
    for (unsigned i = 0; i < num_byte; i++)
        b[i] ^= a[i];
}

void gf256v_madd(uint8_t *accu_c, const uint8_t *a, uint8_t gf256_b,
                 unsigned num_byte)
{
    gf256_vector_madd(accu_c, a, gf256_b, num_byte);
}

void gf256v_mul_scalar(uint8_t *a, uint8_t b, unsigned num_byte)
{
    for (unsigned i = 0; i < num_byte; i++)
        a[i] = gf256_mul(a[i], b);
}

void gf256v_set_zero(uint8_t *a, unsigned num_byte)
{
    memset(a, 0, num_byte);
}

void gf256v_conditional_add(uint8_t *accu_b, uint8_t condition,
                            const uint8_t *a, unsigned num_byte)
{
    uint8_t mask = (uint8_t)(0u - (condition & 1u));
    for (unsigned i = 0; i < num_byte; i++)
        accu_b[i] ^= (uint8_t)(a[i] & mask);
}
