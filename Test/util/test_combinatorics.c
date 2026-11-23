#include "util/combinatorics.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

enum { MAX_N = 12 };
static uint64_t pascal[MAX_N + 1][MAX_N + 1];
static uint32_t offsets[MAX_N + 1][MAX_N + 2];
static uint64_t data[MAX_N + 1][MAX_N + 1];
static pqc_binomial_table tables[MAX_N + 1];

static int enumerate(uint32_t *expected, uint32_t position, uint32_t start,
                      uint32_t n, uint32_t k, uint32_t *rank)
{
    if (position == k) {
        pqc_u256 value;
        uint32_t decoded[MAX_N];
        pqc_u256_set(&value, *rank);
        CHECK(pqc_combination_unrank(decoded, n, k, &value, NULL) == 0);
        CHECK(memcmp(decoded, expected, k * sizeof(*expected)) == 0);
        CHECK(pqc_combination_unrank(decoded, n, k, &value, tables) == 0);
        CHECK(memcmp(decoded, expected, k * sizeof(*expected)) == 0);
        ++*rank;
        return 0;
    }
    for (uint32_t i = start; i <= n - (k - position); ++i) {
        expected[position] = i;
        CHECK(enumerate(expected, position + 1, i + 1, n, k, rank) == 0);
    }
    return 0;
}

int main(void)
{
    for (unsigned n = 0; n <= MAX_N; ++n) {
        pascal[n][0] = pascal[n][n] = 1;
        for (unsigned k = 1; k < n; ++k)
            pascal[n][k] = pascal[n - 1][k - 1] + pascal[n - 1][k];
    }
    for (unsigned k = 2; k <= MAX_N; ++k) {
        tables[k] = (pqc_binomial_table){offsets[k], data[k]};
        for (unsigned n = k; n <= MAX_N; ++n) {
            offsets[k][n - k] = n - k;
            data[k][n - k] = pascal[n][k];
        }
        offsets[k][MAX_N - k + 1] = MAX_N - k + 1;
    }
    for (uint32_t n = 0; n <= MAX_N; ++n) {
        for (uint32_t k = 0; k <= n; ++k) {
            pqc_u256 binomial;
            CHECK(pqc_u256_binomial(&binomial, n, k) == 0);
            CHECK(binomial.limb[0] == pascal[n][k]);
            CHECK(binomial.limb[1] == 0 && binomial.limb[2] == 0 && binomial.limb[3] == 0);
            uint32_t expected[MAX_N], rank = 0;
            CHECK(enumerate(expected, 0, 0, n, k, &rank) == 0);
            CHECK(rank == pascal[n][k]);
            pqc_u256_set(&binomial, rank);
            CHECK(pqc_combination_unrank(expected, n, k, &binomial, NULL) == -1);
        }
    }
    const pqc_u256 known = {{UINT64_C(0x2114f400aad36330), UINT64_C(0x5b1df70a0a6e5924),
                             UINT64_C(0xbf41a3ef5a038323), UINT64_C(0x1386ceab9061)}};
    pqc_u256 value;
    CHECK(pqc_u256_binomial(&value, 4721, 27) == 0);
    CHECK(pqc_u256_cmp(&value, &known) == 0);
    CHECK(pqc_u256_binomial(&value, 300, 150) == -1);
    CHECK(pqc_u256_cmp(&value, &known) == 0);
    CHECK(pqc_u256_binomial(&value, 3, 4) == 0 && pqc_u256_is_zero(&value));

    uint32_t subset[27];
    pqc_u256_set(&value, 0);
    CHECK(pqc_combination_unrank(subset, 4721, 27, &value, NULL) == 0);
    for (unsigned i = 0; i < 27; ++i)
        CHECK(subset[i] == i);
    value = known;
    pqc_u256 one;
    pqc_u256_set(&one, 1);
    pqc_u256_sub_assign(&value, &one);
    CHECK(pqc_combination_unrank(subset, 4721, 27, &value, NULL) == 0);
    for (unsigned i = 0; i < 27; ++i)
        CHECK(subset[i] == 4721 - 27 + i);
    CHECK(pqc_combination_unrank(NULL, 0, 0, &one, NULL) == -1);
    pqc_u256_set(&value, 0);
    CHECK(pqc_combination_unrank(NULL, 0, 0, &value, NULL) == 0);

    /* Four-byte big-endian numbers reduced modulo 65521. */
    pqc_u256 residues[4][256], modulus;
    pqc_u256_set(&modulus, 65521);
    uint64_t weight = 1;
    for (size_t position = 4; position-- > 0;) {
        for (unsigned b = 0; b < 256; ++b)
            pqc_u256_set(&residues[position][b], weight * b % 65521);
        weight = weight * 256 % 65521;
    }
    uint32_t random = 0x654321fe;
    for (unsigned trial = 0; trial < 1024; ++trial) {
        random = random * 1664525u + 1013904223u;
        uint8_t bytes[] = {(uint8_t)(random >> 24), (uint8_t)(random >> 16),
                           (uint8_t)(random >> 8), (uint8_t)random};
        pqc_u256_reduce_bytes_table(&value, bytes, sizeof(bytes), residues, &modulus);
        CHECK(value.limb[0] == random % 65521 && value.limb[1] == 0);
    }
    pqc_u256_reduce_bytes_table(&value, NULL, 0, NULL, &modulus);
    CHECK(pqc_u256_is_zero(&value));
    return 0;
}
