#include "poly/fp16.h"

#include <stddef.h>

/* Newton iteration doubles the number of correct low bits, so four steps
 * lift the inverse of an odd q from one bit to sixteen bits. */
static uint16_t inverse_mod_2_16(int16_t q)
{
    uint32_t x = 1;
    uint32_t qv = (uint32_t)(uint16_t)q;

    for (int i = 0; i < 4; ++i) {
        x *= 2u - qv * x;
    }
    return (uint16_t)x;
}

static int is_prime(int16_t q)
{
    int32_t n = q;

    if (n < 3) {
        return 0;
    }
    for (int32_t d = 2; d * d <= n; ++d) {
        if (n % d == 0) {
            return 0;
        }
    }
    return 1;
}

int fp16_init(fp16_ctx *ctx, int16_t q)
{
    fp16_ctx out;

    if (ctx == NULL || !is_prime(q)) {
        return -1;
    }

    out.q = q;
    out.qinv = inverse_mod_2_16(q);
    out.barrett = (uint16_t)(32768u / (uint32_t)q);
    /* R^2 mod q, computed as (R mod q)^2 mod q so nothing exceeds 32 bits. */
    {
        uint32_t r = 65536u % (uint32_t)q;
        out.mont_r2 = (uint16_t)((r * r) % (uint32_t)q);
    }
    *ctx = out;
    return 0;
}
