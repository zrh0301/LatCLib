#include "fp32.h"

#include <stddef.h>

int fp32_init(fp32_ctx *ctx, int32_t q)
{
    if (ctx == NULL || q < 3 || (q & 1) == 0) {
        return -1;
    }
    for (int32_t d = 3; d <= q / d; d += 2) {
        if (q % d == 0) {
            return -1;
        }
    }

    uint32_t qinv = 1;
    for (unsigned int i = 0; i < 5; ++i) {
        qinv = (uint32_t)((uint64_t)qinv *
                         (UINT64_C(2) - (uint64_t)(uint32_t)q * qinv));
    }
    int32_t r = (int32_t)(UINT64_C(4294967296) % (uint32_t)q);
    fp32_ctx result = {
        q, qinv, r, (int32_t)((uint64_t)r * (uint32_t)r % (uint32_t)q)
    };
    *ctx = result;
    return 0;
}

int32_t fp32_reduce_raw(int64_t a, int32_t q)
{
    int64_t r = a % q;
    return (int32_t)(r + (int64_t)(r < 0) * q);
}

int32_t fp32_center_raw(int64_t a, int32_t q)
{
    int32_t r = fp32_reduce_raw(a, q);
    return r - (int32_t)(r > q / 2) * q;
}

int32_t fp32_reduce(int64_t a, const fp32_ctx *ctx)
{
    return fp32_reduce_raw(a, ctx->q);
}

int32_t fp32_center(int64_t a, const fp32_ctx *ctx)
{
    return fp32_center_raw(a, ctx->q);
}

int32_t fp32_add(int32_t a, int32_t b, const fp32_ctx *ctx)
{
    return fp32_reduce((int64_t)a + b, ctx);
}

int32_t fp32_sub(int32_t a, int32_t b, const fp32_ctx *ctx)
{
    return fp32_reduce((int64_t)a - b, ctx);
}

int32_t fp32_mul(int32_t a, int32_t b, const fp32_ctx *ctx)
{
    return fp32_reduce((int64_t)a * b, ctx);
}

int fp32_inv(int32_t *out, int32_t a, const fp32_ctx *ctx)
{
    int32_t base = fp32_reduce(a, ctx);
    if (out == NULL || base == 0) {
        return -1;
    }
    int32_t result = 1;
    uint32_t exponent = (uint32_t)ctx->q - 2;
    while (exponent != 0) {
        if ((exponent & 1) != 0) {
            result = fp32_mul(result, base, ctx);
        }
        base = fp32_mul(base, base, ctx);
        exponent >>= 1;
    }
    *out = result;
    return 0;
}

int32_t fp32_to_mont(int32_t a, const fp32_ctx *ctx)
{
    return fp32_reduce(fp32_mont_reduce((int64_t)a * ctx->mont_r2, ctx), ctx);
}

int32_t fp32_from_mont(int32_t a, const fp32_ctx *ctx)
{
    return fp32_reduce(fp32_mont_reduce(a, ctx), ctx);
}
