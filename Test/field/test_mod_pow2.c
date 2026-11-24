#include "gf/mod_pow2.h"

#include <stdio.h>
#include <stdint.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "mod_pow2 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

int main(void)
{
    /* Reduction is a mask for every supported width. */
    for (unsigned int logq = 1; logq <= 16; ++logq) {
        const uint32_t q = 1u << logq;
        for (uint32_t v = 0; v < 4 * q; ++v)
            CHECK(pqc_mod_pow2(v, logq) == (uint16_t)(v % q));
        CHECK(pqc_mod_pow2(UINT32_MAX, logq) == (uint16_t)(q - 1));
    }

    /* Rounding reproduces round-half-up division by 2^(logq-bits). */
    for (unsigned int logq = 2; logq <= 16; ++logq) {
        const uint32_t q = 1u << logq;
        for (unsigned int bits = 1; bits < logq; ++bits) {
            const uint32_t shift = 1u << (logq - bits);
            for (uint32_t v = 0; v < q; v += 7) {
                const uint32_t want = (v + shift / 2) / shift;
                CHECK(pqc_round_pow2((uint16_t)v, logq, bits) == (int32_t)want);
            }
        }
    }

    /* Scloud+ usage: q = 2^10. Round-half-up can reach 2^bits. */
    CHECK(pqc_mod_pow2(1234, 10) == 210);
    CHECK(pqc_round_pow2(1023, 10, 1) == 2);
    CHECK(pqc_round_pow2(255, 10, 1) == 0);
    CHECK(pqc_round_pow2(256, 10, 1) == 1);

    puts("mod_pow2: masking and rounding passed");
    return 0;
}
