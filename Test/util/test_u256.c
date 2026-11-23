#include "util/bytes.h"
#include "util/u256.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

int main(void)
{
    pqc_u256 x, y, term, mod;
    const uint64_t lo[] = {UINT64_MAX, UINT64_MAX};
    pqc_u256_load_limbs(&x, lo, 2);
    pqc_u256_set(&y, 1);
    CHECK(pqc_u256_add_assign(&x, &y) == 0);
    CHECK(x.limb[2] == 1 && x.limb[1] == 0 && x.limb[0] == 0);
    CHECK(pqc_u256_significant_limbs(&x) == 3);
    CHECK(pqc_u256_cmp_limbs(lo, 2, &x, 3) < 0);
    CHECK(pqc_u256_sub(&x, &x, &y) == 0);
    CHECK(pqc_u256_cmp_limbs(lo, 2, &x, 2) == 0);
    uint32_t v32;
    uint64_t v64;
    CHECK(pqc_u256_to_u32_checked(&v32, &x) != 0);
    CHECK(pqc_u256_to_u64_checked(&v64, &x) != 0);
    CHECK(pqc_u256_to_u32_checked(&v32, &y) == 0 && v32 == 1);
    CHECK(pqc_u256_to_u64_checked(&v64, &y) == 0 && v64 == 1);

    pqc_u256_zero(&term);
    CHECK(pqc_u256_sub(&mod, &term, &y) == 1);
    CHECK(pqc_u256_sub(&x, &mod, &y) == 0);
    pqc_u256_add_mod_assign(&x, &x, &mod);
    CHECK(x.limb[0] == UINT64_MAX - 2 && x.limb[3] == UINT64_MAX);
    pqc_u256_set(&mod, 17);
    pqc_u256_set(&x, 16);
    pqc_u256_add_mod_assign(&x, &x, &mod);
    CHECK(x.limb[0] == 15);
    pqc_u256_set(&x, 7);
    pqc_u256_add_mod_assign(&x, &x, &mod);
    CHECK(x.limb[0] == 14);

    const uint64_t limbs[] = {UINT64_MAX, 0x8e88b2ab3a37c394ULL, 5};
    const uint32_t factors[] = {1, 2, 3, 6, 0x80000000u, UINT32_MAX};
    pqc_u256_load_limbs(&term, limbs, 3);
    for (size_t i = 0; i < sizeof(factors) / sizeof(*factors); ++i) {
        x = term;
        pqc_u256_mul_u32(&x, factors[i]);
        pqc_u256_div_u32_exact(&x, factors[i]);
        CHECK(pqc_u256_cmp(&x, &term) == 0);
    }
    pqc_u256_set(&x, UINT64_MAX);
    pqc_u256_mul_u32(&x, UINT32_MAX);
    CHECK(x.limb[0] == UINT64_MAX - UINT32_MAX + 1u);
    CHECK(x.limb[1] == UINT32_MAX - 1u);
    pqc_u256_mul_u32(&x, 0);
    CHECK(pqc_u256_is_zero(&x));
    CHECK(pqc_u32_gcd(18, 12) == 6 && pqc_u32_gcd(0, 7) == 7);

    const uint8_t expected[] = {0, 0, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0};
    uint8_t bytes[sizeof(expected)];
    pqc_store_be(bytes, sizeof(bytes), UINT64_C(0x123456789abcdef0));
    CHECK(memcmp(bytes, expected, sizeof(bytes)) == 0);
    CHECK(pqc_load_be(bytes + 2, 8) == UINT64_C(0x123456789abcdef0));
    pqc_store_be(bytes, 3, UINT64_C(0x123456789abcdef0));
    CHECK(pqc_load_be(bytes, 3) == UINT64_C(0xbcdef0));
    pqc_store_be(NULL, 0, 1);
    CHECK(pqc_load_be(NULL, 0) == 0);
    return 0;
}
