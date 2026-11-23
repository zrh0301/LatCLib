#include "poly/poly16.h"

/* Canonical CBD_eta over one contiguous bit segment: coefficient i consumes the
 * 2*eta bits starting at bit 2*eta*i, the first eta bits contribute positively
 * and the next eta negatively. Bits are little-endian inside each byte, so a
 * 32-bit rolling accumulator covers a whole field plus one top-up byte. */
/* The rolling accumulator tops up at most one byte beyond the widest field
 * supported here; eta <= 15 keeps the field within 16 bits. */
#define CBD_MAX_FIELD_BITS (2u * 8u)
_Static_assert(CBD_MAX_FIELD_BITS + 8u <= 32u, "CBD accumulator overflow");

static void cbd_segment(int16_t *r, size_t n, const uint8_t *buf, unsigned eta,
                        int accumulate)
{
    const unsigned width = 2u * eta;
    const uint32_t mask = (1u << eta) - 1u;
    const uint32_t field = (1u << width) - 1u;
    uint32_t acc = 0;
    unsigned have = 0;
    size_t pos = 0;

    for (size_t i = 0; i < n; ++i) {
        while (have < width) {
            acc |= (uint32_t)buf[pos++] << have;
            have += 8;
        }
        const uint32_t v = acc & field;
        acc >>= width;
        have -= width;

        const int16_t s = (int16_t)(__builtin_popcount(v & mask) -
                                    __builtin_popcount((v >> eta) & mask));
        r[i] = accumulate ? (int16_t)(r[i] + s) : s;
    }
}

void fp16_poly_add(int16_t *r, const int16_t *a, const int16_t *b, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        r[i] = (int16_t)(a[i] + b[i]);
    }
}

void fp16_poly_sub(int16_t *r, const int16_t *a, const int16_t *b, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        r[i] = (int16_t)(a[i] - b[i]);
    }
}

void fp16_poly_barrett(int16_t *r, size_t n, const fp16_ctx *ctx)
{
    for (size_t i = 0; i < n; ++i) {
        r[i] = fp16_barrett_reduce(r[i], ctx);
    }
}

void fp16_poly_caddq(int16_t *r, size_t n, const fp16_ctx *ctx)
{
    for (size_t i = 0; i < n; ++i) {
        r[i] = fp16_caddq(r[i], ctx);
    }
}

void fp16_poly_caddq2(int16_t *r, size_t n, const fp16_ctx *ctx)
{
    for (size_t i = 0; i < n; ++i) {
        r[i] = fp16_caddq2(r[i], ctx);
    }
}

void fp16_poly_tomont(int16_t *r, size_t n, const fp16_ctx *ctx)
{
    for (size_t i = 0; i < n; ++i) {
        r[i] = fp16_montgomery_reduce((int32_t)r[i] * (int32_t)ctx->mont_r2, ctx);
    }
}

void fp16_poly_cbd(int16_t *r, size_t n, const uint8_t *buf, unsigned eta)
{
    unsigned hi = 1;
    unsigned lo;

    if (eta == 0) {
        for (size_t i = 0; i < n; ++i) {
            r[i] = 0;
        }
        return;
    }
    while ((hi << 1) <= eta) {
        hi <<= 1;
    }
    lo = eta - hi;

    cbd_segment(r, n, buf, hi, 0);
    if (lo != 0) {
        cbd_segment(r, n, buf + (n * hi) / 4u, lo, 1);
    }
}

size_t fp16_poly_rej_uniform12(int16_t *out, size_t n, const uint8_t *in,
                               size_t len, const fp16_ctx *ctx)
{
    size_t count = 0;
    while (len >= 3 && count < n) {
        const int16_t a = (int16_t)((in[0] | ((unsigned)in[1] << 8)) & 0xfffu);
        const int16_t b = (int16_t)((in[1] >> 4) | ((unsigned)in[2] << 4));
        if (a < ctx->q)
            out[count++] = a;
        if (b < ctx->q && count < n)
            out[count++] = b;
        in += 3;
        len -= 3;
    }
    return count;
}
