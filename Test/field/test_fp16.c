#include "poly/fp16.h"

#include <stddef.h>
#include <stdio.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "fp16 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

/* Independent int64_t oracle; a and q*q products stay far below 2^63. */
static int32_t reference_mod(int64_t a, int32_t q)
{
    int64_t r = a % q;
    if (r < 0) {
        r += q;
    }
    return (int32_t)r;
}

/* Extended Euclid, deliberately different from the field's Newton iteration. */
static int32_t reference_inverse(int32_t a, int32_t q)
{
    int64_t r = q, next_r = reference_mod(a, q);
    int64_t t = 0, next_t = 1;
    while (next_r != 0) {
        int64_t quotient = r / next_r;
        int64_t remainder = r - quotient * next_r;
        int64_t coefficient = t - quotient * next_t;
        r = next_r;
        next_r = remainder;
        t = next_t;
        next_t = coefficient;
    }
    return reference_mod(t, q);
}

static uint64_t next_random(void)
{
    static uint64_t state = UINT64_C(0x3136703136703136);
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

/* The reference reductions are only valid on two's-complement targets; the
 * kernels document the same contract. */
static int16_t reference_montgomery(int32_t a, const fp16_ctx *ctx)
{
    int64_t t = (int64_t)(int16_t)a * (int64_t)ctx->qinv;
    t &= 0xffff;
    if (t >= 32768) {
        t -= 65536;
    }
    return (int16_t)(((int64_t)a - t * ctx->q) / 65536);
}

static int16_t reference_barrett(int16_t a, const fp16_ctx *ctx)
{
    int64_t u = ((int64_t)a * ctx->barrett + (1 << 14)) >> 15;
    return (int16_t)((int64_t)a - u * ctx->q);
}

static int check_montgomery(int64_t a, const fp16_ctx *ctx)
{
    int16_t result = fp16_montgomery_reduce((int32_t)a, ctx);
    int32_t r_inverse = reference_inverse(reference_mod(65536, ctx->q), ctx->q);

    CHECK(result == reference_montgomery((int32_t)a, ctx));
    CHECK(result >= -ctx->q && result <= ctx->q);
    CHECK(reference_mod(result, ctx->q) ==
          reference_mod(a * r_inverse, ctx->q));
    return 0;
}

static int check_field(int16_t q)
{
    fp16_ctx ctx;
    CHECK(fp16_init(&ctx, q) == 0);
    CHECK(ctx.q == q);
    CHECK((uint32_t)((uint32_t)(uint16_t)q * ctx.qinv) % 65536u == 1u);
    CHECK(ctx.barrett == 32768u / (uint32_t)q);

    /* Montgomery contract: result = a*R^{-1} mod q for the full documented
     * input range, including both boundaries and negative inputs. */
    const int64_t bound = (int64_t)32768 * q;
    const int64_t mont_inputs[] = {
        -bound, -bound + 1, bound - 1, bound,
        -(int64_t)q - 1, -q, -1, 0, 1, q - 1, q, (int64_t)q + 1,
        INT32_MIN, INT32_MAX
    };
    for (size_t i = 0; i < sizeof(mont_inputs) / sizeof(mont_inputs[0]); ++i) {
        if (mont_inputs[i] >= -bound && mont_inputs[i] <= bound) {
            CHECK(check_montgomery(mont_inputs[i], &ctx) == 0);
        }
    }
    for (unsigned int i = 0; i < 4000; ++i) {
        int64_t wide = (int64_t)(next_random() % (2 * (uint64_t)bound + 1)) - bound;
        CHECK(check_montgomery(wide, &ctx) == 0);
    }

    /* Barrett reduces into (-q, q) and preserves the residue. The documented
     * preconditions admit |a| <= barrett*q. */
    const int32_t barrett_limit = (int32_t)ctx.barrett * q;
    const int16_t barrett_inputs[] = {
        (int16_t)-barrett_limit, (int16_t)(-q), (int16_t)(-1), 0,
        (int16_t)(q / 2), (int16_t)(q - 1), (int16_t)barrett_limit
    };
    for (size_t i = 0; i < sizeof(barrett_inputs) / sizeof(barrett_inputs[0]); ++i) {
        int16_t a = barrett_inputs[i];
        int16_t result = fp16_barrett_reduce(a, &ctx);
        CHECK(result == reference_barrett(a, &ctx));
        CHECK(result >= -q && result <= q);
        CHECK(reference_mod(result, q) == reference_mod(a, q));
    }

    /* Conditional addition normalises the lazy ranges used by the NTT. */
    int32_t start = -2 * (int32_t)q;
    if (start < -32768) {
        start = -32768;
    }
    for (int32_t a = start; a < q; ++a) {
        int16_t result = fp16_caddq((int16_t)a, &ctx);
        if (a >= -q) {
            CHECK(result >= 0 && result < q);
            CHECK(reference_mod(result, q) == reference_mod(a, q));
        }
        int16_t twice = fp16_caddq2((int16_t)a, &ctx);
        CHECK(twice >= 0 && twice < q);
        CHECK(reference_mod(twice, q) == reference_mod(a, q));
    }
    return 0;
}

int main(void)
{
    const int16_t invalid[] = {
        -32768, -3329, -1, 0, 1, 2, 4, 9, 15, 341, 561, 1105, 32767
    };
    fp16_ctx ctx = {17, 1, 1, 1};
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        CHECK(fp16_init(&ctx, invalid[i]) == -1);
        CHECK(ctx.q == 17 && ctx.qinv == 1 && ctx.barrett == 1 && ctx.mont_r2 == 1);
    }
    CHECK(fp16_init(NULL, 3329) == -1);

    const int16_t primes[] = {3, 5, 17, 641, 3329, 32749};
    for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); ++i) {
        CHECK(check_field(primes[i]) == 0);
    }

    CHECK(fp16_init(&ctx, 3329) == 0);
    CHECK(ctx.qinv == 62209 && ctx.barrett == 9);
    CHECK(fp16_init(&ctx, 641) == 0);
    CHECK(ctx.qinv == 15745 && ctx.barrett == 51);
    puts("fp16: prime field, Montgomery/Barrett reductions and caddq passed");
    return 0;
}
