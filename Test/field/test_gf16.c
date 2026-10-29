#include "gf16.h"
#include "gf16_matrix.h"

#include <stdio.h>
#include <string.h>

#define CHECK(expr) do { \
    if (!(expr)) { \
        fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr); \
        return 1; \
    } \
} while (0)

static uint8_t reference_mul(uint8_t a, uint8_t b) {
    uint8_t result = 0;
    for (int i = 0; i < 4; ++i) {
        if (b & 1u) {
            result ^= a;
        }
        b >>= 1;
        a <<= 1;
        if (a & 0x10u) {
            a ^= 0x13u;
        }
    }
    return result;
}

int main(void) {
    for (uint8_t a = 0; a < 16; ++a) {
        CHECK(gf16_add(a, a) == 0);
        CHECK(gf16_sub(a, a) == 0);
        for (uint8_t b = 0; b < 16; ++b) {
            CHECK(gf16_add(a, b) == (uint8_t)(a ^ b));
            CHECK(gf16_sub(a, b) == (uint8_t)(a ^ b));
            CHECK(gf16_mul(a, b) == reference_mul(a, b));
        }
        if (a != 0) {
            CHECK(reference_mul(a, gf16_inv(a)) == 1);
        } else {
            CHECK(gf16_inv(a) == 0);
        }
    }

    {
        const uint8_t a[6] = {1, 2, 3, 4, 5, 6};
        const uint8_t b[6] = {7, 8, 9, 10, 11, 12};
        const uint8_t expected[4] = {
            (uint8_t)(reference_mul(1, 7) ^ reference_mul(2, 9) ^ reference_mul(3, 11)),
            (uint8_t)(reference_mul(1, 8) ^ reference_mul(2, 10) ^ reference_mul(3, 12)),
            (uint8_t)(reference_mul(4, 7) ^ reference_mul(5, 9) ^ reference_mul(6, 11)),
            (uint8_t)(reference_mul(4, 8) ^ reference_mul(5, 10) ^ reference_mul(6, 12))
        };
        uint8_t out[4] = {0};
        CHECK(gf16_matrix_mul(out, a, b, 2, 3, 2) == 0);
        CHECK(memcmp(out, expected, sizeof(out)) == 0);
        CHECK(gf16_matrix_mul_add(out, a, b, 2, 3, 2) == 0);
        CHECK(memcmp(out, (uint8_t[4]){0}, sizeof(out)) == 0);
    }

    for (size_t n = 1; n <= 5; ++n) {
        uint8_t identity[25] = {0};
        uint8_t upper[25] = {0};
        uint8_t det = 0xff;
        uint8_t expected_det = 1;
        for (size_t i = 0; i < n; ++i) {
            identity[i * n + i] = 1;
            upper[i * n + i] = (uint8_t)(i + 2);
            expected_det = reference_mul(expected_det, upper[i * n + i]);
            for (size_t j = i + 1; j < n; ++j) {
                upper[i * n + j] = (uint8_t)(i + j + 1);
            }
        }
        CHECK(gf16_matrix_determinant(&det, identity, n) == 0);
        CHECK(det == 1);
        CHECK(gf16_matrix_determinant(&det, upper, n) == 0);
        CHECK(det == expected_det);
        if (n > 1) {
            memcpy(identity + n, identity, n);
            CHECK(gf16_matrix_determinant(&det, identity, n) == 0);
            CHECK(det == 0);
        }
    }

    {
        uint8_t sentinel = 0x5a;
        CHECK(gf16_matrix_determinant(&sentinel, NULL, 2) == -1);
        CHECK(sentinel == 0x5a);
        CHECK(gf16_matrix_determinant(&sentinel, (uint8_t[36]){0}, 6) == -1);
        CHECK(sentinel == 0x5a);
    }

    puts("GF(16) and matrix operators passed");
    return 0;
}
