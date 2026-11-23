#include "util/combinatorics.h"

int pqc_u256_binomial(pqc_u256 *out, uint32_t n, uint32_t k)
{
    pqc_u256 value;
    if (k > n) {
        pqc_u256_zero(out);
        return 0;
    }
    if (k > n - k)
        k = n - k;
    pqc_u256_set(&value, 1);
    for (uint32_t i = 1; i <= k; ++i) {
        uint32_t numerator = n - k + i;
        const uint32_t common = pqc_u32_gcd(numerator, i);
        numerator /= common;
        pqc_u256_div_u32_exact(&value, i / common);
        uint64_t carry = 0;
        for (size_t j = 0; j < PQC_U256_LIMBS; ++j)
            value.limb[j] = pqc_u64_mul_u32_carry(value.limb[j], numerator, &carry);
        if (carry != 0)
            return -1;
    }
    *out = value;
    return 0;
}

int pqc_u256_binomial_get(pqc_u256 *out, uint32_t n, uint32_t k,
                          const pqc_binomial_table *tables)
{
    if (tables == NULL)
        return pqc_u256_binomial(out, n, k);
    if (k > n) {
        pqc_u256_zero(out);
    } else if (k == 0) {
        pqc_u256_set(out, 1);
    } else if (k == 1) {
        pqc_u256_set(out, n);
    } else {
        const pqc_binomial_table *table = tables + k;
        const uint32_t offset = table->offsets[n - k];
        const uint32_t len = table->offsets[n - k + 1] - offset;
        pqc_u256_load_limbs(out, table->data + offset, len);
    }
    return 0;
}

int pqc_combination_unrank(uint32_t *out, uint32_t n, uint32_t k,
                           const pqc_u256 *rank, const pqc_binomial_table *tables)
{
    pqc_u256 remaining, base, target, binomial;
    if (!rank || (k != 0 && !out) || k > n ||
        pqc_u256_binomial_get(&base, n, k, tables) != 0 ||
        pqc_u256_cmp(rank, &base) >= 0)
        return -1;

    remaining = *rank;
    uint32_t previous = 0;
    for (uint32_t i = 0; i < k; ++i) {
        const uint32_t count = k - i;
        uint32_t low = count;
        uint32_t high = n - previous;
        if (pqc_u256_binomial_get(&base, high, count, tables) != 0)
            return -1;
        pqc_u256_sub(&target, &base, &remaining);

        while (low < high) {
            const uint32_t mid = low + (high - low) / 2;
            if (pqc_u256_binomial_get(&binomial, mid, count, tables) != 0)
                return -1;
            if (pqc_u256_cmp(&binomial, &target) >= 0)
                high = mid;
            else
                low = mid + 1;
        }
        if (pqc_u256_binomial_get(&binomial, low, count, tables) != 0)
            return -1;
        pqc_u256_sub(&remaining, &binomial, &target);
        out[i] = n - low;
        previous = out[i] + 1;
    }
    return pqc_u256_is_zero(&remaining) ? 0 : -1;
}

void pqc_u256_reduce_bytes_table(pqc_u256 *out, const uint8_t *in, size_t len,
                                 const pqc_u256 table[][256], const pqc_u256 *modulus)
{
    pqc_u256_zero(out);
    for (size_t i = 0; i < len; ++i) {
        if (in[i] != 0)
            pqc_u256_add_mod_assign(out, &table[i][in[i]], modulus);
    }
}
