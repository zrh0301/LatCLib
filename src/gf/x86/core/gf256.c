#include "gfext.h"

#include <string.h>

/* GF(256) = GF(2)[x] / (x^8 + x^4 + x^3 + x + 1). */
uint8_t pqc_kernel_gf256_mul_ref(uint8_t a, uint8_t b) {
    uint8_t result = 0;
    for (unsigned i = 0; i < 8; ++i) {
        const uint8_t select = (uint8_t)-(int)(b & 1u);
        const uint8_t carry = (uint8_t)(a >> 7);
        result ^= (uint8_t)(a & select);
        a = (uint8_t)(a << 1);
        a ^= (uint8_t)(0x1bu & (uint8_t)-(int)carry);
        b >>= 1;
    }
    return result;
}

/* The MQOM reference implementation evaluates four independent GF(256)
 * products in one 32-bit word.  Each byte lane is reduced independently; the
 * arithmetic is branch-free and preserves the 0x11b field. */
static inline uint32_t pqc_kernel_gf256_lane_xtime_x4(uint32_t v) {
    const uint32_t msb = v & 0x80808080u;
    const uint32_t shifted = (v & 0x7f7f7f7fu) << 1;
    const uint32_t reduction = (msb - (msb >> 7)) & 0x1b1b1b1bu;
    return shifted ^ reduction;
}

static inline uint32_t pqc_kernel_gf256_mul4_swar(uint32_t x, uint32_t y) {
    uint32_t acc = 0;
    uint32_t a = x;
#define GF256_SWAR_STEP(bit) do { \
        uint32_t mask = (y >> (bit)) & 0x01010101u; \
        mask |= (mask << 1) & 0xfefefefeu; \
        mask |= (mask << 2) & 0xfcfcfcfcu; \
        mask |= (mask << 4) & 0xf0f0f0f0u; \
        acc ^= a & mask; \
        a = pqc_kernel_gf256_lane_xtime_x4(a); \
    } while (0)
    GF256_SWAR_STEP(0);
    GF256_SWAR_STEP(1);
    GF256_SWAR_STEP(2);
    GF256_SWAR_STEP(3);
    GF256_SWAR_STEP(4);
    GF256_SWAR_STEP(5);
    GF256_SWAR_STEP(6);
    GF256_SWAR_STEP(7);
#undef GF256_SWAR_STEP
    return acc;
}

uint8_t pqc_kernel_gf256_mul_swar_ref(uint8_t a, uint8_t b) {
    return (uint8_t)pqc_kernel_gf256_mul4_swar((uint32_t)a, (uint32_t)b);
}

uint32_t pqc_kernel_gf256_mul4_swar_word(uint32_t x, uint32_t y) {
    return pqc_kernel_gf256_mul4_swar(x, y);
}

void pqc_kernel_gf256_mul_const_swar_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    const uint32_t kv = (uint32_t)k * 0x01010101u;
    size_t i = 0;
    for (; i + 4u <= n; i += 4u) {
        uint32_t av, product;
        memcpy(&av, a + i, sizeof(av));
        product = pqc_kernel_gf256_mul4_swar(av, kv);
        memcpy(out + i, &product, sizeof(product));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_swar_ref(a[i], k);
}

void pqc_kernel_gf256_mulacc_swar_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + 4u <= n; i += 4u) {
        uint32_t av, bv, cv, product;
        memcpy(&av, a + i, sizeof(av));
        memcpy(&bv, b + i, sizeof(bv));
        memcpy(&cv, acc + i, sizeof(cv));
        product = pqc_kernel_gf256_mul4_swar(av, bv) ^ cv;
        memcpy(acc + i, &product, sizeof(product));
    }
    for (; i < n; ++i) acc[i] ^= pqc_kernel_gf256_mul_swar_ref(a[i], b[i]);
}

uint8_t pqc_kernel_gf256_dot_swar_array(const uint8_t *a, const uint8_t *b, size_t n) {
    uint32_t acc = 0;
    size_t i = 0;
    for (; i + 4u <= n; i += 4u) {
        uint32_t av, bv;
        memcpy(&av, a + i, sizeof(av));
        memcpy(&bv, b + i, sizeof(bv));
        acc ^= pqc_kernel_gf256_mul4_swar(av, bv);
    }
    uint8_t r = (uint8_t)(acc ^ (acc >> 8) ^ (acc >> 16) ^ (acc >> 24));
    for (; i < n; ++i) r ^= pqc_kernel_gf256_mul_swar_ref(a[i], b[i]);
    return r;
}

/* ------------------------------------------------- 求逆、平方、幂、异或的原语 */

/* 256 项 log/exp 表：GF(256) 的非零元构成 255 阶循环群（0x03 是本原元）。 */
static uint8_t pqc_kernel_gf256_log_table[256];
static uint8_t pqc_kernel_gf256_exp_table[256];
static int pqc_kernel_gf256_logexp_ready = 0;

static void pqc_kernel_gf256_logexp_init(void) {
    if (pqc_kernel_gf256_logexp_ready) return;
    uint8_t x = 1;
    for (unsigned e = 0; e < 255; ++e) {
        pqc_kernel_gf256_exp_table[e] = x;
        pqc_kernel_gf256_log_table[x] = (uint8_t)e;
        x = pqc_kernel_gf256_mul_ref(x, 0x03u);
    }
    pqc_kernel_gf256_log_table[0] = 0;
    pqc_kernel_gf256_exp_table[255] = pqc_kernel_gf256_exp_table[0];
    pqc_kernel_gf256_logexp_ready = 1;
}

/* a^-1 = exp[(255 - log a) mod 255]，约定 inv(0) = 0。 */
uint8_t pqc_kernel_gf256_inv_logexp_ref(uint8_t a) {
    pqc_kernel_gf256_logexp_init();
    if (a == 0) return 0;
    return pqc_kernel_gf256_exp_table[(255u - pqc_kernel_gf256_log_table[a]) % 255u];
}

void pqc_kernel_gf256_inv_logexp_array(uint8_t *out, const uint8_t *a, size_t n) {
    pqc_kernel_gf256_logexp_init();
    for (size_t i = 0; i < n; ++i) {
        const uint8_t v = a[i];
        out[i] = v ? pqc_kernel_gf256_exp_table[(255u - pqc_kernel_gf256_log_table[v]) % 255u] : 0;
    }
}

uint8_t pqc_kernel_gf256_pow_ref(uint8_t a, uint32_t e) {
    uint8_t result = 1;
    uint8_t base = a;
    while (e) {
        if (e & 1u) result = pqc_kernel_gf256_mul_ref(result, base);
        e >>= 1;
        if (e) base = pqc_kernel_gf256_mul_ref(base, base);
    }
    return result;
}

void pqc_kernel_gf256_pow_scalar_array(uint8_t *out, const uint8_t *a, uint32_t e, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = pqc_kernel_gf256_pow_ref(a[i], e);
}

void pqc_kernel_gf256_square_ref_array(uint8_t *out, const uint8_t *a, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(a[i], a[i]);
}

void pqc_kernel_gf256_xor_scalar_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = (uint8_t)(a[i] ^ b[i]);
}

void pqc_kernel_gf256_mul_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = pqc_kernel_gf256_mul_ref(a[i], b[i]);
    }
}

void pqc_kernel_gf256_mul_swar_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + 4u <= n; i += 4u) {
        uint32_t av, bv;
        memcpy(&av, a + i, sizeof(av));
        memcpy(&bv, b + i, sizeof(bv));
        const uint32_t product = pqc_kernel_gf256_mul4_swar(av, bv);
        memcpy(out + i, &product, sizeof(product));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_swar_ref(a[i], b[i]);
}

void pqc_kernel_gf256_mul_const_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = pqc_kernel_gf256_mul_ref(a[i], k);
    }
}

void pqc_kernel_gf256_mulacc_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        acc[i] ^= pqc_kernel_gf256_mul_ref(a[i], b[i]);
    }
}

/* 向量异或：GF(256) 的加法就是逐字节异或，是常数乘之外的“零成本”底座。 */
void pqc_kernel_gf256_xor_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i)));
    }
    for (; i < n; ++i) out[i] = (uint8_t)(a[i] ^ b[i]);
}
