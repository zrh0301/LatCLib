#include "poly/poly16.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

static int test_repeated_bits(const fp16_ctx *ctx, unsigned copies, size_t chunk_bytes)
{
    const uint8_t msg[] = {0x96, 0x45, 0xac, 0x53, 0xff, 0x00, 0x81, 0x7e};
    const size_t stride = 8 * chunk_bytes;
    const int16_t midpoint = ctx->q / 2;
    int16_t coeffs[8 * 4 * sizeof(msg)];
    uint8_t decoded[sizeof(msg)];

    fp16_poly_encode_repeated_bits(coeffs, msg, sizeof(msg), chunk_bytes, copies, ctx);
    for (size_t byte = 0; byte < sizeof(msg); ++byte) {
        const size_t offset = (byte / chunk_bytes) * copies * stride;
        for (unsigned copy = 0; copy < copies; ++copy)
            for (unsigned bit = 0; bit < 8; ++bit) {
                const size_t index = offset + copy * stride + 8 * (byte % chunk_bytes) + bit;
                const int16_t expected = (msg[byte] >> bit) & 1u ? midpoint + 1 : 0;
                CHECK(coeffs[index] == expected);
                /* Add small noise using both centered and nonnegative representatives. */
                coeffs[index] = (int16_t)(coeffs[index] + (int)copy - 1);
            }
    }
    fp16_poly_decode_repeated_bits(decoded, coeffs, sizeof(msg), chunk_bytes, copies, ctx);
    CHECK(memcmp(decoded, msg, sizeof(msg)) == 0);

    /* An exact threshold decodes to 0. One unit below it decodes to 1. */
    for (size_t i = 0; i < copies * stride; ++i)
        coeffs[i] = midpoint;
    const unsigned threshold = copies * (unsigned)ctx->q / 4;
    for (unsigned bit = 0; bit < 2; ++bit) {
        unsigned remaining = threshold - bit;
        for (unsigned copy = 0; copy < copies; ++copy) {
            const unsigned distance = remaining < (unsigned)midpoint
                                      ? remaining : (unsigned)midpoint;
            coeffs[copy * stride + bit] = (int16_t)(midpoint - (int)distance);
            remaining -= distance;
        }
        CHECK(remaining == 0);
    }
    fp16_poly_decode_repeated_bits(decoded, coeffs, chunk_bytes, chunk_bytes, copies, ctx);
    CHECK(decoded[0] == 0xfe);
    return 0;
}

int main(void)
{
    const int16_t coeffs[] = {0, 1, 3328, 1664, 1665, 17, 1500};
    const size_t n = sizeof(coeffs) / sizeof(*coeffs);
    uint8_t encoded[16], expected[16];
    int16_t decoded[7];
    fp16_ctx ctx;
    CHECK(fp16_init(&ctx, 3329) == 0);
    CHECK(test_repeated_bits(&ctx, 2, 8) == 0);
    CHECK(test_repeated_bits(&ctx, 4, 4) == 0);

    for (unsigned bits = 1; bits <= 15; ++bits) {
        const unsigned mask = (1u << bits) - 1u;
        for (unsigned compress = 0; compress <= 1; ++compress) {
            memset(expected, 0, sizeof(expected));
            memset(encoded, 0x55, sizeof(encoded));
            for (size_t i = 0; i < n; ++i) {
                unsigned c = (unsigned)coeffs[i];
                if (compress)
                    c = ((c << bits) + 1664u) / 3329u;
                c &= mask;
                for (size_t j = 0; j < bits; ++j) {
                    size_t pos = i * bits + j;
                    expected[pos / 8] |= (uint8_t)(((c >> j) & 1u) << (pos % 8));
                }
            }
            if (compress) {
                fp16_poly_compress(encoded, coeffs, n, bits, &ctx);
                fp16_poly_decompress(decoded, encoded, n, bits, &ctx);
            } else {
                fp16_poly_encode(encoded, coeffs, n, bits);
                fp16_poly_decode(decoded, encoded, n, bits);
            }
            const size_t bytes = (n * bits + 7) / 8;
            CHECK(memcmp(encoded, expected, bytes) == 0 && encoded[bytes] == 0x55);
            for (size_t i = 0; i < n; ++i) {
                unsigned c = (unsigned)coeffs[i];
                if (compress) {
                    c = (((c << bits) + 1664u) / 3329u) & mask;
                    c = (c * 3329u + (1u << (bits - 1))) >> bits;
                }
                CHECK(decoded[i] == (int16_t)(c & (compress ? 0x7fff : mask)));
            }
        }
    }

    /* Candidates in pairs: (4095, 0), (3328, 3329), (1, 256), plus tail. */
    const uint8_t candidates[] = {0xff, 0x0f, 0, 0, 0x1d, 0xd0, 1, 0, 0x10, 2, 0};
    const int16_t accepted[] = {0, 3328, 1, 256};
    CHECK(fp16_poly_rej_uniform12(decoded, n, candidates, sizeof(candidates), &ctx) == 4);
    CHECK(memcmp(decoded, accepted, sizeof(accepted)) == 0);
    decoded[1] = -1;
    CHECK(fp16_poly_rej_uniform12(decoded, 1, candidates, sizeof(candidates), &ctx) == 1);
    CHECK(decoded[0] == 0 && decoded[1] == -1);
    CHECK(fp16_poly_rej_uniform12(NULL, 0, NULL, 0, &ctx) == 0);
    return 0;
}
