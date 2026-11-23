#include "poly_sample.h"
#include <stdlib.h>

int poly_rand_with_rng(lssi_t *a, int coefs, pqc_rng_func rng, void *rng_ctx)
{
    uint8_t buf[MAX_POLY_COEFS * 4];
    if (a == NULL || coefs < 0 || coefs > MAX_POLY_COEFS ||
        pqc_randombytes(buf, (size_t)coefs * 4, rng, rng_ctx) != 0) {
        return -1;
    }
    for (int i = 0; i < coefs; ++i) {
        uint32_t x = ((uint32_t)buf[4*i] << 24) |
                     ((uint32_t)buf[4*i+1] << 16) |
                     ((uint32_t)buf[4*i+2] << 8) | buf[4*i+3];
        a[i] = (int32_t)((int64_t)x -
                        (int64_t)(x >> 31) * INT64_C(4294967296));
    }
    return 0;
}

void poly_rand(lssi_t *a, int coefs)
{
    if (poly_rand_with_rng(a, coefs, NULL, NULL) != 0) abort();
}

unsigned int sample_CBD(int32_t *a, unsigned int len, const uint8_t *buf,
                        unsigned int buflen, lp_ctx *ctx)
{
    unsigned int count = 0;
    (void)ctx;
    for (unsigned int pos = 0; pos < buflen && count < len; ++pos) {
        uint32_t t = buf[pos];
        uint32_t d = (t & 0x55u) + ((t >> 1) & 0x55u);
        for (unsigned int i = 0; i < 2 && count < len; ++i) {
            a[count++] = (int32_t)((d >> (4*i)) & 3u) -
                         (int32_t)((d >> (4*i+2)) & 3u);
        }
    }
    return count;
}

int poly_uniform_CBD_with_rng(poly *a, lp_ctx *ctx,
                              pqc_rng_func rng, void *rng_ctx)
{
    uint8_t buf[MAX_POLY_COEFS / 2];
    if (a == NULL || ctx == NULL || ctx->ring.n != MAX_POLY_COEFS ||
        pqc_randombytes(buf, sizeof(buf), rng, rng_ctx) != 0) {
        return -1;
    }
    (void)sample_CBD(*a, MAX_POLY_COEFS, buf, sizeof(buf), ctx);
    return 0;
}

void poly_uniform_CBD(poly *a, lp_ctx *ctx)
{
    if (poly_uniform_CBD_with_rng(a, ctx, NULL, NULL) != 0) abort();
}
