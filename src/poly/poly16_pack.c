#include "poly/poly16.h"

static void encode(uint8_t *out, const int16_t *coeffs, size_t n,
                   unsigned bits, uint16_t q)
{
    const uint32_t mask = (1u << bits) - 1u;
    uint32_t word = 0;
    unsigned pending = 0;

    for (size_t i = 0; i < n; ++i) {
        uint32_t coeff = (uint16_t)coeffs[i];
        if (q == 3329)
            coeff = ((coeff << bits) + 1664u) / 3329u;
        else if (q != 0)
            coeff = ((coeff << bits) + q / 2u) / q;

        word |= (coeff & mask) << pending;
        pending += bits;
        while (pending >= 8) {
            *out++ = (uint8_t)word;
            word >>= 8;
            pending -= 8;
        }
    }
    if (pending != 0)
        *out = (uint8_t)word;
}

static void decode(int16_t *out, const uint8_t *in, size_t n,
                   unsigned bits, uint16_t q)
{
    const uint32_t mask = (1u << bits) - 1u;
    uint32_t word = 0;
    unsigned pending = 0;

    for (size_t i = 0; i < n; ++i) {
        while (pending < bits) {
            word |= (uint32_t)*in++ << pending;
            pending += 8;
        }
        uint32_t coeff = word & mask;
        if (q != 0)
            coeff = (coeff * q + (1u << (bits - 1))) >> bits;
        out[i] = (int16_t)coeff;
        word >>= bits;
        pending -= bits;
    }
}

void fp16_poly_encode(uint8_t *out, const int16_t *in, size_t n, unsigned bits)
{
    encode(out, in, n, bits, 0);
}

void fp16_poly_decode(int16_t *out, const uint8_t *in, size_t n, unsigned bits)
{
    decode(out, in, n, bits, 0);
}

void fp16_poly_compress(uint8_t *out, const int16_t *in, size_t n,
                         unsigned bits, const fp16_ctx *ctx)
{
    encode(out, in, n, bits, (uint16_t)ctx->q);
}

void fp16_poly_decompress(int16_t *out, const uint8_t *in, size_t n,
                           unsigned bits, const fp16_ctx *ctx)
{
    decode(out, in, n, bits, (uint16_t)ctx->q);
}

void fp16_poly_encode_repeated_bits(int16_t *out, const uint8_t *msg,
                                    size_t msg_bytes, size_t chunk_bytes,
                                    unsigned copies, const fp16_ctx *ctx)
{
    const size_t stride = 8 * chunk_bytes;
    const int16_t midpoint = (int16_t)((ctx->q + 1) / 2);
    for (size_t chunk = 0; chunk < msg_bytes; chunk += chunk_bytes) {
        for (size_t byte = 0; byte < chunk_bytes; ++byte) {
            for (unsigned bit = 0; bit < 8; ++bit) {
                const int16_t mask = (int16_t)-(int)((msg[chunk + byte] >> bit) & 1u);
                for (unsigned copy = 0; copy < copies; ++copy)
                    out[copy * stride + 8 * byte + bit] = mask & midpoint;
            }
        }
        out += copies * stride;
    }
}

void fp16_poly_decode_repeated_bits(uint8_t *out, const int16_t *coeffs,
                                    size_t msg_bytes, size_t chunk_bytes,
                                    unsigned copies, const fp16_ctx *ctx)
{
    const size_t stride = 8 * chunk_bytes;
    const int32_t midpoint = ctx->q / 2;
    const uint32_t threshold = copies * (uint32_t)ctx->q / 4u;

    for (size_t chunk = 0; chunk < msg_bytes; chunk += chunk_bytes) {
        for (size_t byte = 0; byte < chunk_bytes; ++byte) {
            uint8_t decoded = 0;
            for (unsigned bit = 0; bit < 8; ++bit) {
                uint32_t distance = 0;
                for (unsigned copy = 0; copy < copies; ++copy) {
                    const int16_t x = coeffs[copy * stride + 8 * byte + bit];
                    const int32_t delta = fp16_caddq(x, ctx) - midpoint;
                    const int32_t mask = -(int32_t)((uint32_t)delta >> 31);
                    distance += (uint32_t)((delta ^ mask) - mask);
                }
                decoded |= (uint8_t)(((distance - threshold) >> 31) << bit);
            }
            out[chunk + byte] = decoded;
        }
        coeffs += copies * stride;
    }
}
