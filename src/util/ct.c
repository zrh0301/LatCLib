#include "util/ct.h"

int pqc_ct_equal(const uint8_t *a, const uint8_t *b, size_t n)
{
    uint32_t diff = 0;

    for (size_t i = 0; i < n; ++i) {
        diff |= (uint32_t)(a[i] ^ b[i]);
    }
    /* Map 0 -> 1 and any nonzero value -> 0 without a branch. */
    return (int)(1u ^ ((diff | (0u - diff)) >> 31));
}

void pqc_ct_select(uint8_t *out, const uint8_t *x, size_t n, uint8_t b)
{
    uint8_t mask = (uint8_t)(0u - (uint32_t)(b & 1u));

    for (size_t i = 0; i < n; ++i) {
        out[i] ^= (uint8_t)(mask & (uint8_t)(x[i] ^ out[i]));
    }
}
