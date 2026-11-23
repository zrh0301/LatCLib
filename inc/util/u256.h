#ifndef PQC_U256_H
#define PQC_U256_H

#include <stddef.h>
#include <stdint.h>

#define PQC_U256_LIMBS 4u

/* Unsigned arithmetic on four little-endian 64-bit limbs.
 * Comparisons and exact division are intended for public values. */
typedef struct {
    uint64_t limb[PQC_U256_LIMBS];
} pqc_u256;

static inline void pqc_u256_zero(pqc_u256 *x)
{
    size_t i;

    for (i = 0; i < PQC_U256_LIMBS; i++) {
        x->limb[i] = 0;
    }
}

static inline void pqc_u256_set(pqc_u256 *x, uint64_t v)
{
    size_t i;

    x->limb[0] = v;
    for (i = 1; i < PQC_U256_LIMBS; i++) {
        x->limb[i] = 0;
    }
}

/* src_limbs <= PQC_U256_LIMBS. */
static inline void pqc_u256_load_limbs(pqc_u256 *out,
                                         const uint64_t *src,
                                         uint32_t src_limbs)
{
    uint32_t i;

    for (i = 0; i < src_limbs; i++) {
        out->limb[i] = src[i];
    }
    for (; i < PQC_U256_LIMBS; i++) {
        out->limb[i] = 0;
    }
}

static inline int pqc_u256_is_zero(const pqc_u256 *x)
{
    size_t i;
    uint64_t acc = 0;

    for (i = 0; i < PQC_U256_LIMBS; i++) {
        acc |= x->limb[i];
    }

    return acc == 0;
}

static inline int pqc_u256_to_u32_checked(uint32_t *out, const pqc_u256 *x)
{
    size_t i;

    if (x->limb[0] > UINT32_MAX) {
        return -1;
    }
    for (i = 1; i < PQC_U256_LIMBS; i++) {
        if (x->limb[i] != 0) {
            return -1;
        }
    }

    *out = (uint32_t)x->limb[0];
    return 0;
}

static inline int pqc_u256_to_u64_checked(uint64_t *out, const pqc_u256 *x)
{
    size_t i;

    for (i = 1; i < PQC_U256_LIMBS; i++) {
        if (x->limb[i] != 0) {
            return -1;
        }
    }

    *out = x->limb[0];
    return 0;
}

static inline int pqc_u256_cmp(const pqc_u256 *a, const pqc_u256 *b)
{
    size_t i = PQC_U256_LIMBS;

    while (i > 0) {
        i--;
        if (a->limb[i] < b->limb[i]) {
            return -1;
        }
        if (a->limb[i] > b->limb[i]) {
            return 1;
        }
    }

    return 0;
}

static inline uint64_t pqc_u256_add_assign(pqc_u256 *a, const pqc_u256 *b)
{
    size_t i;
    uint64_t carry = 0;

    for (i = 0; i < PQC_U256_LIMBS; i++) {
        const uint64_t ai = a->limb[i];
        const uint64_t bi = b->limb[i];
        const uint64_t tmp = ai + carry;
        const uint64_t carry1 = (tmp < ai);
        const uint64_t res = tmp + bi;
        const uint64_t carry2 = (res < tmp);

        a->limb[i] = res;
        carry = carry1 | carry2;
    }

    return carry;
}

static inline uint64_t pqc_u256_sub(pqc_u256 *out,
                                      const pqc_u256 *a,
                                      const pqc_u256 *b)
{
    size_t i;
    uint64_t borrow = 0;

    for (i = 0; i < PQC_U256_LIMBS; i++) {
        const uint64_t ai = a->limb[i];
        const uint64_t bi = b->limb[i];
        const uint64_t tmp = ai - borrow;
        const uint64_t borrow1 = (ai < borrow);
        const uint64_t res = tmp - bi;
        const uint64_t borrow2 = (tmp < bi);

        out->limb[i] = res;
        borrow = borrow1 | borrow2;
    }

    return borrow;
}

static inline uint64_t pqc_u256_sub_assign(pqc_u256 *a, const pqc_u256 *b)
{
    return pqc_u256_sub(a, a, b);
}

/* acc = (acc + term) mod mod. Requires 0 <= acc, term < mod and mod < 2^(64*limbs). */
static inline void pqc_u256_add_mod_assign(pqc_u256 *acc,
                                             const pqc_u256 *term,
                                             const pqc_u256 *mod)
{
    const uint64_t carry = pqc_u256_add_assign(acc, term);

    if (carry != 0u || pqc_u256_cmp(acc, mod) >= 0) {
        (void)pqc_u256_sub_assign(acc, mod);
    }
}

/* Low 64 bits and carry of x*m + carry, with carry < 2^32. */
static inline uint64_t pqc_u64_mul_u32_carry(uint64_t x, uint32_t m,
                                              uint64_t *carry)
{
    const uint64_t low = (uint64_t)(uint32_t)x * m + *carry;
    const uint64_t high = (x >> 32) * m + (low >> 32);
    *carry = high >> 32;
    return (high << 32) | (uint32_t)low;
}

/* Truncates the product modulo 2^256. */
static inline void pqc_u256_mul_u32(pqc_u256 *x, uint32_t m)
{
    uint64_t carry = 0;
    size_t i;

    for (i = 0; i < PQC_U256_LIMBS; i++) {
        x->limb[i] = pqc_u64_mul_u32_carry(x->limb[i], m, &carry);
    }
}

static inline uint32_t pqc_u32_gcd(uint32_t a, uint32_t b)
{
    while (b != 0u) {
        const uint32_t t = a % b;
        a = b;
        b = t;
    }

    return a;
}

/* x must be nonzero. */
static inline uint32_t pqc_u32_ctz(uint32_t x)
{
    uint32_t n = 0;

    while ((x & 1u) == 0u) {
        x >>= 1;
        n++;
    }

    return n;
}

/* Shift by 0 through 63 bits. */
static inline void pqc_u256_rshift_bits(pqc_u256 *x, uint32_t s)
{
    uint64_t carry = 0;
    size_t i = PQC_U256_LIMBS;

    if (s == 0u) {
        return;
    }

    while (i > 0) {
        const uint64_t limb = x->limb[i - 1u];
        x->limb[i - 1u] = (limb >> s) | carry;
        carry = limb << (64u - s);
        i--;
    }
}

static inline uint64_t pqc_inv_odd_u64(uint32_t d)
{
    uint64_t inv = d;

    inv *= 2u - (uint64_t)d * inv;
    inv *= 2u - (uint64_t)d * inv;
    inv *= 2u - (uint64_t)d * inv;
    inv *= 2u - (uint64_t)d * inv;
    inv *= 2u - (uint64_t)d * inv;
    inv *= 2u - (uint64_t)d * inv;

    return inv;
}

static inline void pqc_u256_div_odd_u32_exact(pqc_u256 *x, uint32_t d)
{
    const uint64_t inv = pqc_inv_odd_u64(d);
    uint64_t carry = 0;
    size_t i;

    for (i = 0; i < PQC_U256_LIMBS; i++) {
        const uint64_t ai = x->limb[i];
        const uint64_t qi = (ai - carry) * inv;
        (void)pqc_u64_mul_u32_carry(qi, d, &carry);
        x->limb[i] = qi;
    }
}

/* Exact division by a non-zero 32-bit integer. The input must be divisible by d. */
static inline void pqc_u256_div_u32_exact(pqc_u256 *x, uint32_t d)
{
    uint32_t shift;

    if (d == 1u) {
        return;
    }

    shift = pqc_u32_ctz(d);
    if (shift != 0u) {
        pqc_u256_rshift_bits(x, shift);
        d >>= shift;
    }

    if (d != 1u) {
        pqc_u256_div_odd_u32_exact(x, d);
    }
}

static inline uint32_t pqc_u256_significant_limbs(const pqc_u256 *x)
{
    uint32_t i = PQC_U256_LIMBS;

    while (i > 0u && x->limb[i - 1u] == 0u) {
        i--;
    }

    return i;
}

static inline int pqc_u256_cmp_limbs(const uint64_t *data,
                                                  uint32_t len,
                                                  const pqc_u256 *target,
                                                  uint32_t target_len)
{
    if (len < target_len) {
        return -1;
    }
    if (len > target_len) {
        return 1;
    }

    while (len > 0u) {
        const uint64_t ai = data[len - 1u];
        const uint64_t bi = target->limb[len - 1u];

        if (ai < bi) {
            return -1;
        }
        if (ai > bi) {
            return 1;
        }
        len--;
    }

    return 0;
}

#endif
