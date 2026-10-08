#ifndef FP32_H
#define FP32_H

#include <stdint.h>

/* Prime field with 3 <= q <= INT32_MAX. All constants are immutable after init. */
typedef struct {
    int32_t q;
    uint32_t qinv; /* POSITIVE q^{-1} modulo R = 2^32. */
    int32_t mont_r; /* R modulo q, in [0, q). */
    int32_t mont_r2; /* R^2 modulo q, in [0, q). */
} fp32_ctx;

/* Validate a public odd prime q. Return 0 on success, -1 on invalid input.
 * Failure leaves ctx unchanged. Initialization is variable-time. */
int fp32_init(fp32_ctx *ctx, int32_t q);

/* These operations require a successfully initialized context. They accept
 * arbitrary signed inputs of the declared width and return canonical [0, q)
 * values, except center, which returns [-(q-1)/2, (q-1)/2].
 * This portable reference backend uses integer remainder and does not promise
 * constant-time execution on every target. */
int32_t fp32_reduce(int64_t a, const fp32_ctx *ctx);
int32_t fp32_center(int64_t a, const fp32_ctx *ctx);
/* Raw normalization helpers for callers that already have a public modulus.
 * Require q > 0; primality and Montgomery constants are not needed. center_raw
 * returns [-floor((q-1)/2), floor(q/2)], choosing +q/2 for an even-modulus tie. */
int32_t fp32_reduce_raw(int64_t a, int32_t q);
int32_t fp32_center_raw(int64_t a, int32_t q);
int32_t fp32_add(int32_t a, int32_t b, const fp32_ctx *ctx);
int32_t fp32_sub(int32_t a, int32_t b, const fp32_ctx *ctx);
int32_t fp32_mul(int32_t a, int32_t b, const fp32_ctx *ctx);

/* Return 0 and write the canonical inverse, or -1 if out is NULL or a == 0
 * modulo q. Failure leaves out unchanged. The exponent depends only on q. */
int fp32_inv(int32_t *out, int32_t a, const fp32_ctx *ctx);

/* Canonical conversions: to_mont(a) = a*R mod q; from_mont(a) = a/R mod q. */
int32_t fp32_to_mont(int32_t a, const fp32_ctx *ctx);
int32_t fp32_from_mont(int32_t a, const fp32_ctx *ctx);

/* Fast lazy Montgomery reduction. q must be positive, odd and <= INT32_MAX,
 * and qinv must be the POSITIVE inverse of q modulo R = 2^32.
 * Precondition: -2^31*q <= a <= 2^31*q.
 * Returns a/R modulo q in [-q, q], without canonicalization.
 * The unsigned low word is mapped explicitly into [-2^31, 2^31-1]. The
 * subtraction fits int64_t under the input bound, and is divisible by R;
 * division therefore needs neither a signed right shift nor rounding rules.
 * This routine contains no data-dependent branches or table accesses. */
static inline int32_t fp32_mont_reduce_raw(int64_t a, int32_t q,
                                         uint32_t qinv)
{
    uint32_t low = (uint32_t)((uint64_t)(uint32_t)a * qinv);
    int64_t t = (int64_t)low - (int64_t)(low >> 31) * INT64_C(4294967296);
    return (int32_t)((a - t * q) / INT64_C(4294967296));
}

static inline int32_t fp32_mont_reduce(int64_t a, const fp32_ctx *ctx)
{
    return fp32_mont_reduce_raw(a, ctx->q, ctx->qinv);
}

#endif /* FP32_H */
