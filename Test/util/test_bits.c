#include "util/bits.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "bits check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

static uint32_t next_random(void)
{
    static uint32_t state = UINT32_C(0x9e3779b9);
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

/* Bit-by-bit oracle: bit i of value goes to bit (offset + i), LSB first. */
static void reference_store(uint8_t *dst, size_t offset, uint32_t value,
                            unsigned int nbits)
{
    for (unsigned int i = 0; i < nbits; ++i) {
        size_t pos = offset + i;
        if ((value >> i) & 1u)
            dst[pos >> 3] |= (uint8_t)(1u << (pos & 7u));
    }
}

static uint32_t reference_load(const uint8_t *src, size_t offset,
                               unsigned int nbits)
{
    uint32_t value = 0;
    for (unsigned int i = 0; i < nbits; ++i) {
        size_t pos = offset + i;
        value |= (uint32_t)((src[pos >> 3] >> (pos & 7u)) & 1u) << i;
    }
    return value;
}

int main(void)
{
    uint8_t got[32], want[32];

    /* Every width against the bitwise oracle, at every byte phase. */
    for (unsigned int nbits = 1; nbits <= 32; ++nbits) {
        for (size_t offset = 0; offset < 24; ++offset) {
            for (int round = 0; round < 8; ++round) {
                uint32_t value = next_random();
                if (!(value & UINT32_C(1)))
                    value = UINT32_MAX - value;

                memset(got, 0, sizeof(got));
                memset(want, 0, sizeof(want));
                pqc_store_bits_le(got, offset, value, nbits);
                reference_store(want, offset, value, nbits);

                CHECK(memcmp(got, want, sizeof(got)) == 0);
                CHECK(pqc_load_bits_le(got, offset, nbits) ==
                      (value & pqc_bit_mask(nbits)));
                CHECK(pqc_load_bits_le(want, offset, nbits) ==
                      reference_load(want, offset, nbits));
            }
        }
    }

    /* The masks used by the candidates: 5-bit field elements and 18-bit chunks. */
    CHECK(pqc_bit_mask(5) == 0x1Fu);
    CHECK(pqc_bit_mask(18) == 0x3FFFFu);
    CHECK(pqc_bit_mask(32) == UINT32_MAX);

    /* Sequential cursor packing is what TSUOV and Scloud+ both do. */
    memset(got, 0, sizeof(got));
    memset(want, 0, sizeof(want));
    size_t got_bit = 0, want_bit = 0;
    for (int i = 0; i < 40; ++i) {
        uint32_t value = (uint32_t)i % 31u;
        pqc_store_bits_le(got, got_bit, value, 5);
        reference_store(want, want_bit, value, 5);
        got_bit += 5;
        want_bit += 5;
    }
    CHECK(memcmp(got, want, sizeof(got)) == 0);
    for (int i = 0; i < 40; ++i)
        CHECK(pqc_load_bits_le(got, (size_t)i * 5, 5) == (uint32_t)i % 31u);

    /* Adjacent writes that share a byte must not clobber each other. */
    memset(got, 0, sizeof(got));
    pqc_store_bits_le(got, 0, 0x1F, 5);
    pqc_store_bits_le(got, 5, 0x07, 3);
    CHECK(got[0] == 0xFF);
    CHECK(pqc_load_bits_le(got, 0, 5) == 0x1F);
    CHECK(pqc_load_bits_le(got, 5, 3) == 0x07);

    puts("bits: LSB-first packing matches the bitwise oracle");
    return 0;
}
