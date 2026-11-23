#ifndef FP16_H
#define FP16_H

#include <stdint.h>

/* Prime field with an odd prime modulus 3 <= q <= INT16_MAX, used by
 * algorithms whose coefficient vectors are stored as int16_t. This is the
 * 16-bit counterpart of fp32.h: Montgomery reduction uses R = 2^16, and
 * arithmetic follows the same lazy-reduction conventions.
 *
 * All kernels are inlined because they sit inside NTT butterflies. They keep
 * the exact result of the Kyber-style reference reductions: the value of
 * `>>` on negative operands is architecture-defined, and these kernels are
 * only used on two's-complement targets where it is an arithmetic shift. */
typedef struct {
    int16_t q;        /* odd prime, 0 < q <= INT16_MAX */
    uint16_t qinv;    /* q^{-1} modulo 2^16 */
    uint16_t barrett; /* floor(2^15 / q) */
    uint16_t mont_r2; /* R^2 modulo q, R = 2^16 */
} fp16_ctx;

/* Validate q and derive qinv/barrett. Return 0 on success or -1 for an even
 * or out-of-range modulus; failure leaves ctx unchanged. */
int fp16_init(fp16_ctx *ctx, int16_t q);

/* Montgomery reduction of a/R modulo q, R = 2^16.
 * Precondition: -2^15*q <= a <= 2^15*q. Returns (-q, q). */
static inline int16_t fp16_montgomery_reduce(int32_t a, const fp16_ctx *ctx)
{
    int32_t t = (int16_t)((int32_t)(int16_t)a * (int32_t)ctx->qinv);
    return (int16_t)((a - t * (int32_t)ctx->q) >> 16);
}

/* Barrett reduction of a into (-q, q).
 * Precondition: |a| <= barrett*q for the configured modulus. */
static inline int16_t fp16_barrett_reduce(int16_t a, const fp16_ctx *ctx)
{
    int32_t u = ((int32_t)a * (int32_t)ctx->barrett + (1 << 14)) >> 15;
    u *= (int32_t)ctx->q;
    return (int16_t)(a - u);
}

/* Conditional add of q: [-q, q) -> [0, q). */
static inline int16_t fp16_caddq(int16_t x, const fp16_ctx *ctx)
{
    return (int16_t)(x + ((x >> 15) & ctx->q));
}

/* Conditional add of q twice: [-2q, q) -> [0, q). */
static inline int16_t fp16_caddq2(int16_t x, const fp16_ctx *ctx)
{
    int16_t r = fp16_caddq(x, ctx);
    return fp16_caddq(r, ctx);
}

#endif /* FP16_H */
