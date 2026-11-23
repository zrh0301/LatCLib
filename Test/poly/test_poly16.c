#include "poly/poly16.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "poly16 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

#define Q 3329

static uint64_t next_random(void)
{
    static uint64_t state = UINT64_C(0x706f6c793136706f);
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static int32_t mod_q(int64_t a)
{
    int64_t r = a % Q;
    return (int32_t)(r < 0 ? r + Q : r);
}

/* Independent bit-at-a-time CBD over one contiguous segment. */
static int16_t cbd_reference(const uint8_t *buf, unsigned eta, size_t i)
{
    const size_t bit = i * 2u * (size_t)eta;
    int16_t v = 0;

    for (unsigned k = 0; k < eta; ++k) {
        const size_t p = bit + k;
        v += (int16_t)((buf[p / 8] >> (p % 8)) & 1u);
    }
    for (unsigned k = 0; k < eta; ++k) {
        const size_t p = bit + eta + k;
        v -= (int16_t)((buf[p / 8] >> (p % 8)) & 1u);
    }
    return v;
}

/* Reference split sampler: high power-of-two segment first, then the rest. */
static int16_t cbd_split_reference(const uint8_t *buf, size_t n, unsigned eta, size_t i)
{
    unsigned hi = 1;
    while ((hi << 1) <= eta) {
        hi <<= 1;
    }
    int16_t v = cbd_reference(buf, hi, i);
    if (eta != hi) {
        v += cbd_reference(buf + (n * hi) / 4u, eta - hi, i);
    }
    return v;
}

int main(void)
{
    fp16_ctx ctx;
    int16_t a[512], b[512], r[512];
    uint8_t buf[6 * 512 / 4 + 8];

    CHECK(fp16_init(&ctx, Q) == 0);
    CHECK(ctx.mont_r2 == 1353);

    for (int i = 0; i < 512; ++i) {
        a[i] = (int16_t)(next_random() % Q);
        b[i] = (int16_t)(next_random() % Q);
    }

    /* add/sub keep exact integer values and are inverses of each other. */
    fp16_poly_add(r, a, b, 512);
    for (int i = 0; i < 512; ++i) {
        CHECK(r[i] == (int16_t)(a[i] + b[i]));
    }
    fp16_poly_sub(r, r, b, 512);
    for (int i = 0; i < 512; ++i) {
        CHECK(r[i] == a[i]);
    }

    /* Barrett reduction lands in (-q, q) and preserves the residue. */
    memcpy(r, a, sizeof(r));
    fp16_poly_barrett(r, 512, &ctx);
    for (int i = 0; i < 512; ++i) {
        CHECK(r[i] > -Q && r[i] < Q);
        CHECK(mod_q(r[i]) == mod_q(a[i]));
    }

    /* caddq/caddq2 normalise the lazy ranges into [0, q). */
    memcpy(r, a, sizeof(r));
    fp16_poly_caddq(r, 512, &ctx);
    for (int i = 0; i < 512; ++i) {
        CHECK(r[i] >= 0 && r[i] < Q);
        CHECK(mod_q(r[i]) == mod_q(a[i]));
    }
    for (int i = 0; i < 512; ++i) {
        r[i] = (int16_t)(a[i] - Q);
    }
    fp16_poly_caddq2(r, 512, &ctx);
    for (int i = 0; i < 512; ++i) {
        CHECK(r[i] >= 0 && r[i] < Q);
        CHECK(mod_q(r[i]) == mod_q(a[i] - Q));
    }

    /* to-montgomery multiplies by R = 2^16 modulo q. */
    memcpy(r, a, sizeof(r));
    fp16_poly_tomont(r, 512, &ctx);
    for (int i = 0; i < 512; ++i) {
        CHECK(mod_q(r[i]) == mod_q((int64_t)a[i] * 65536));
    }

    /* Sampler: eta 1..6 against the bit-by-bit reference, using two calls so
     * the buffer offsets are exercised on independent data too. */
    for (unsigned eta = 1; eta <= 6; ++eta) {
        const size_t bytes = (size_t)eta * 512 / 4;
        for (size_t i = 0; i < bytes; ++i) {
            buf[i] = (uint8_t)next_random();
        }
        fp16_poly_cbd(r, 512, buf, eta);
        for (size_t i = 0; i < 512; ++i) {
            CHECK(r[i] == cbd_split_reference(buf, 512, eta, i));
        }
    }

    /* eta 0 is the degenerate all-zero sample, and must not read the buffer. */
    fp16_poly_cbd(r, 512, NULL, 0);
    for (int i = 0; i < 512; ++i) {
        CHECK(r[i] == 0);
    }

    puts("poly16: arithmetic, to-montgomery and split CBD passed");
    return 0;
}
