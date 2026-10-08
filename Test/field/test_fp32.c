#include "fp32.h"

#include <limits.h>
#include <stdio.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "fp32 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

/* Independent int64_t oracle; products of two int32_t residues fit. */
static int32_t reference_mod(int64_t a, int32_t q)
{
    int64_t r = a % q;
    if (r < 0) {
        r += q;
    }
    return (int32_t)r;
}

/* Extended Euclid, deliberately different from the field's Fermat inversion. */
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
    static uint64_t state = UINT64_C(0x7265666572656e63);
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static int32_t random_i32(void)
{
    uint32_t bits = (uint32_t)next_random();
    return (int32_t)((int64_t)bits -
                     (int64_t)(bits >> 31) * INT64_C(4294967296));
}

static int check_pair(int32_t a, int32_t b, const fp32_ctx *ctx)
{
    CHECK(fp32_add(a, b, ctx) == reference_mod((int64_t)a + b, ctx->q));
    CHECK(fp32_sub(a, b, ctx) == reference_mod((int64_t)a - b, ctx->q));
    CHECK(fp32_mul(a, b, ctx) == reference_mod((int64_t)a * b, ctx->q));

    int32_t am = fp32_to_mont(a, ctx);
    int32_t bm = fp32_to_mont(b, ctx);
    CHECK(am == reference_mod((int64_t)a * ctx->mont_r, ctx->q));
    CHECK(fp32_from_mont(am, ctx) == reference_mod(a, ctx->q));
    int32_t product = fp32_mont_reduce((int64_t)am * bm, ctx);
    CHECK(fp32_from_mont(product, ctx) == reference_mod((int64_t)a * b, ctx->q));
    return 0;
}

static int check_montgomery(int64_t a, int32_t r_inverse, const fp32_ctx *ctx)
{
    int32_t result = fp32_mont_reduce(a, ctx);
    int32_t expected = reference_mod((int64_t)reference_mod(a, ctx->q) *
                                    r_inverse, ctx->q);
    CHECK(result >= -ctx->q && result <= ctx->q);
    CHECK(reference_mod(result, ctx->q) == expected);
    CHECK(fp32_mont_reduce_raw(a, ctx->q, ctx->qinv) == result);
    return 0;
}

static int check_field(int32_t q)
{
    fp32_ctx ctx;
    CHECK(fp32_init(&ctx, q) == 0);
    CHECK(ctx.q == q);
    CHECK((uint32_t)((uint64_t)(uint32_t)q * ctx.qinv) == 1);
    CHECK(ctx.mont_r == reference_mod(INT64_C(4294967296), q));
    CHECK(ctx.mont_r2 == reference_mod((int64_t)ctx.mont_r * ctx.mont_r, q));

    const int64_t reduction_inputs[] = {
        INT64_MIN, INT64_MIN + 1, INT64_MAX, INT64_MAX - 1,
        INT32_MIN, INT32_MAX, -(int64_t)q * q, (int64_t)q * q,
        -(int64_t)q - 1, -q, -1, 0, 1, q - 1, q, (int64_t)q + 1
    };
    for (size_t i = 0; i < sizeof(reduction_inputs) / sizeof(reduction_inputs[0]); ++i) {
        int64_t a = reduction_inputs[i];
        CHECK(fp32_reduce(a, &ctx) == reference_mod(a, q));
        CHECK(fp32_reduce_raw(a, q) == reference_mod(a, q));
        int32_t centered = fp32_center(a, &ctx);
        CHECK(fp32_center_raw(a, q) == centered);
        CHECK(centered >= -(q / 2) && centered <= q / 2);
        CHECK(reference_mod(centered, q) == reference_mod(a, q));
    }

    const int32_t inputs[] = {
        INT32_MIN, INT32_MIN + 1, INT32_MAX, -q, -(q / 2), -1,
        0, 1, q / 2, q - 1, q
    };
    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); ++i) {
        for (size_t j = 0; j < sizeof(inputs) / sizeof(inputs[0]); ++j) {
            CHECK(check_pair(inputs[i], inputs[j], &ctx) == 0);
        }
        int32_t inverse = 123;
        if (reference_mod(inputs[i], q) == 0) {
            CHECK(fp32_inv(&inverse, inputs[i], &ctx) == -1);
            CHECK(inverse == 123);
        } else {
            CHECK(fp32_inv(&inverse, inputs[i], &ctx) == 0);
            CHECK(inverse == reference_inverse(inputs[i], q));
        }
    }
    CHECK(fp32_inv(NULL, 1, &ctx) == -1);

    int32_t r_inverse = reference_inverse(ctx.mont_r, q);
    int64_t bound = INT64_C(2147483648) * q;
    const int64_t mont_inputs[] = {
        -bound, -bound + 1, bound - 1, bound,
        INT32_MIN, INT32_MAX, -q, -1, 0, 1, q
    };
    for (size_t i = 0; i < sizeof(mont_inputs) / sizeof(mont_inputs[0]); ++i) {
        CHECK(check_montgomery(mont_inputs[i], r_inverse, &ctx) == 0);
    }
    /* The positive boundary reduces lazily to q, not to canonical zero. */
    CHECK(fp32_mont_reduce(bound, &ctx) == q);

    for (unsigned int i = 0; i < 4000; ++i) {
        int32_t a = random_i32(), b = random_i32();
        CHECK(check_pair(a, b, &ctx) == 0);
        int64_t wide = (int64_t)(next_random() % (UINT64_C(2) * (uint64_t)bound + 1)) - bound;
        CHECK(check_montgomery(wide, r_inverse, &ctx) == 0);
        if (i < 128 && reference_mod(a, q) != 0) {
            int32_t inverse;
            CHECK(fp32_inv(&inverse, a, &ctx) == 0);
            CHECK(inverse == reference_inverse(a, q));
        }
    }
    return 0;
}

int main(void)
{
    const int32_t invalid[] = {
        INT32_MIN, -1, 0, 1, 2, 4, 9, 15, 341, 561, 1105, 2147483645
    };
    fp32_ctx ctx = {17, UINT32_C(4042322161), 1, 1};
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        CHECK(fp32_init(&ctx, invalid[i]) == -1);
        CHECK(ctx.q == 17 && ctx.qinv == UINT32_C(4042322161));
        CHECK(ctx.mont_r == 1 && ctx.mont_r2 == 1);
    }
    CHECK(fp32_init(NULL, 3329) == -1);

    const int32_t raw_moduli[] = {1, 2, 4, 256};
    const int64_t raw_inputs[] = {INT64_MIN, INT64_MAX, -256, -1, 0, 1, 128, 256};
    for (size_t i = 0; i < sizeof(raw_moduli) / sizeof(raw_moduli[0]); ++i) {
        int32_t q = raw_moduli[i];
        for (size_t j = 0; j < sizeof(raw_inputs) / sizeof(raw_inputs[0]); ++j) {
            int32_t centered = fp32_center_raw(raw_inputs[j], q);
            CHECK(fp32_reduce_raw(raw_inputs[j], q) == reference_mod(raw_inputs[j], q));
            CHECK(centered >= -((q - 1) / 2) && centered <= q / 2);
            CHECK(reference_mod(centered, q) == reference_mod(raw_inputs[j], q));
        }
    }

    const int32_t primes[] = {3, 5, 17, 3329, 65537, 8380417, INT32_MAX};
    for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); ++i) {
        CHECK(check_field(primes[i]) == 0);
    }
    CHECK(fp32_init(&ctx, 3329) == 0);
    CHECK(ctx.qinv == UINT32_C(1806234369) && ctx.mont_r == 1353);
    CHECK(fp32_init(&ctx, 8380417) == 0);
    CHECK(ctx.qinv == UINT32_C(58728449) && ctx.mont_r == 4193792);
    puts("fp32: field arithmetic, inversion and Montgomery contracts passed");
    return 0;
}
