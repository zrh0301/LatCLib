#include "gfext.h"
#include "gfni_field_map.h"

#include <string.h>

uint8_t pqc_kernel_gf2_parity_u64(uint64_t x) {
    x ^= x >> 32; x ^= x >> 16; x ^= x >> 8;
    x ^= x >> 4; x ^= x >> 2; x ^= x >> 1;
    return (uint8_t)(x & 1u);
}

uint8_t pqc_kernel_gf2_dot_u64(uint64_t a, uint64_t b) { return pqc_kernel_gf2_parity_u64(a & b); }

void pqc_kernel_gf2_xor_words(uint64_t *dst, const uint64_t *src, size_t words) {
    for (size_t i = 0; i < words; ++i) dst[i] ^= src[i];
}

uint8_t pqc_kernel_gf2_dot_words(const uint64_t *a, const uint64_t *b, size_t words) {
    uint8_t r = 0;
    for (size_t i = 0; i < words; ++i) r ^= pqc_kernel_gf2_dot_u64(a[i], b[i]);
    return r;
}

void pqc_kernel_gf2_mat_vec_u64(uint64_t *out, const uint64_t *rows, const uint64_t *v,
                     size_t row_count, size_t words_per_row) {
    const size_t out_words = (row_count + 63u) / 64u;
    for (size_t i = 0; i < out_words; ++i) out[i] = 0;
    for (size_t row = 0; row < row_count; ++row) {
        const uint8_t bit = pqc_kernel_gf2_dot_words(rows + row * words_per_row, v, words_per_row);
        out[row >> 6] |= (uint64_t)bit << (row & 63u);
    }
}

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

static inline uint32_t pqc_kernel_gf256_lane_xtime_x4(uint32_t v) {
    const uint32_t msb = v & 0x80808080u;
    const uint32_t shifted = (v & 0x7f7f7f7fu) << 1;
    return shifted ^ ((msb - (msb >> 7)) & 0x1b1b1b1bu);
}

static inline uint32_t pqc_kernel_gf256_mul4_swar(uint32_t x, uint32_t y) {
    uint32_t acc = 0;
    uint32_t a = x;
#define GF256_SWAR_STEP(bit) do { \
        uint32_t mask = (y >> (bit)) & 0x01010101u; \
        mask |= (mask << 1) & 0xfefefefeu; \
        mask |= (mask << 2) & 0xfcfcfcfcu; \
        mask |= (mask << 4) & 0xf0f0f0f0u; \
        acc ^= a & mask; a = pqc_kernel_gf256_lane_xtime_x4(a); \
    } while (0)
    GF256_SWAR_STEP(0); GF256_SWAR_STEP(1); GF256_SWAR_STEP(2); GF256_SWAR_STEP(3);
    GF256_SWAR_STEP(4); GF256_SWAR_STEP(5); GF256_SWAR_STEP(6); GF256_SWAR_STEP(7);
#undef GF256_SWAR_STEP
    return acc;
}

uint8_t pqc_kernel_gf256_mul_swar_ref(uint8_t a, uint8_t b) {
    return (uint8_t)pqc_kernel_gf256_mul4_swar((uint32_t)a, (uint32_t)b);
}

void pqc_kernel_gf256_mul_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(a[i], b[i]);
}

void pqc_kernel_gf256_mul_swar_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + 4u <= n; i += 4u) {
        uint32_t av, bv, product;
        memcpy(&av, a + i, sizeof(av)); memcpy(&bv, b + i, sizeof(bv));
        product = pqc_kernel_gf256_mul4_swar(av, bv); memcpy(out + i, &product, sizeof(product));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_swar_ref(a[i], b[i]);
}

void pqc_kernel_gf256_mul_const_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(a[i], k);
}

void pqc_kernel_gf256_mulacc_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) acc[i] ^= pqc_kernel_gf256_mul_ref(a[i], b[i]);
}

pqc_kernel_gf256_vec pqc_kernel_gf256_vec_mul_neon(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b) {
    pqc_kernel_gf256_vec result = vdupq_n_u8(0);
    const pqc_kernel_gf256_vec one = vdupq_n_u8(1);
    const pqc_kernel_gf256_vec reduction = vdupq_n_u8(0x1b);
    for (unsigned i = 0; i < 8; ++i) {
        const pqc_kernel_gf256_vec take = vceqq_u8(vandq_u8(b, one), one);
        result = veorq_u8(result, vandq_u8(a, take));
        const pqc_kernel_gf256_vec high = vshrq_n_u8(a, 7);
        const pqc_kernel_gf256_vec reduce = vceqq_u8(high, one);
        a = vshlq_n_u8(a, 1);
        a = veorq_u8(a, vandq_u8(reduction, reduce));
        b = vshrq_n_u8(b, 1);
    }
    return result;
}

void pqc_kernel_gf256_mul_neon_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf256_vec_mul_neon(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i)));
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(a[i], b[i]);
}

void pqc_kernel_gf256_mulacc_neon_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec c = pqc_kernel_gf16_vec_load(acc + i);
        const pqc_kernel_gf256_vec p = pqc_kernel_gf256_vec_mul_neon(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i));
        pqc_kernel_gf16_vec_store(acc + i, veorq_u8(c, p));
    }
    for (; i < n; ++i) acc[i] ^= pqc_kernel_gf256_mul_ref(a[i], b[i]);
}

void pqc_kernel_gf256_mul_const_neon_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    const pqc_kernel_gf256_vec kv = vdupq_n_u8(k);
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf256_vec_mul_neon(pqc_kernel_gf16_vec_load(a + i), kv));
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(a[i], k);
}

static uint8_t reduce_vec(pqc_kernel_gf256_vec v) {
    uint8_t tmp[GF16_VEC_BYTES];
    pqc_kernel_gf16_vec_store(tmp, v);
    uint8_t r = 0;
    for (size_t i = 0; i < GF16_VEC_BYTES; ++i) r ^= tmp[i];
    return r;
}

uint8_t pqc_kernel_gf256_dot_neon_array(const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    pqc_kernel_gf256_vec acc = vdupq_n_u8(0);
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        acc = veorq_u8(acc, pqc_kernel_gf256_vec_mul_neon(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i)));
    uint8_t r = reduce_vec(acc);
    for (; i < n; ++i) r ^= pqc_kernel_gf256_mul_ref(a[i], b[i]);
    return r;
}

static pqc_kernel_gf256_vec embed_vec(pqc_kernel_gf256_vec v, const uint8_t *table) {
    return vqtbl1q_u8(vld1q_u8(table), vandq_u8(v, vdupq_n_u8(0x0f)));
}

uint8_t pqc_kernel_gf16_embed_gf256(uint8_t a) { return gfni_gf16_embed_scalar_table[a & 0x0f]; }

void pqc_kernel_gf16_to_gf256_array(uint8_t *out, const uint8_t *in, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(out + i, embed_vec(pqc_kernel_gf16_vec_load(in + i), gfni_gf16_embed_scalar_table));
    for (; i < n; ++i) out[i] = pqc_kernel_gf16_embed_gf256(in[i]);
}

void pqc_kernel_gf16_gf256_mul_array(uint8_t *out, const uint8_t *gf16, const uint8_t *gf256, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec a = embed_vec(pqc_kernel_gf16_vec_load(gf16 + i), gfni_gf16_embed_scalar_table);
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf256_vec_mul_neon(a, pqc_kernel_gf16_vec_load(gf256 + i)));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(pqc_kernel_gf16_embed_gf256(gf16[i]), gf256[i]);
}

void pqc_kernel_gf16_gf256_mulacc_array(uint8_t *acc, const uint8_t *gf16, const uint8_t *gf256, size_t n) {
    for (size_t i = 0; i < n; ++i) acc[i] ^= pqc_kernel_gf256_mul_ref(pqc_kernel_gf16_embed_gf256(gf16[i]), gf256[i]);
}

void pqc_kernel_gf16_gf256_mulacc_neon_array(uint8_t *acc, const uint8_t *gf16,
                                  const uint8_t *gf256, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec a = embed_vec(pqc_kernel_gf16_vec_load(gf16 + i), gfni_gf16_embed_scalar_table);
        const pqc_kernel_gf256_vec p = pqc_kernel_gf256_vec_mul_neon(a, pqc_kernel_gf16_vec_load(gf256 + i));
        pqc_kernel_gf16_vec_store(acc + i, veorq_u8(pqc_kernel_gf16_vec_load(acc + i), p));
    }
    for (; i < n; ++i) acc[i] ^= pqc_kernel_gf256_mul_ref(pqc_kernel_gf16_embed_gf256(gf16[i]), gf256[i]);
}

uint8_t pqc_kernel_mqom_gf16_embed_gf256(uint8_t a) { return pqc_kernel_mqom_gf16_embed_scalar_table[a & 0x0f]; }

static void mqom_lift_block(uint8_t *out, const uint8_t *packed) {
    uint8x8_t p = vld1_u8(packed);
    uint8x8_t lo = vand_u8(p, vdup_n_u8(0x0f));
    uint8x8_t hi = vand_u8(vshr_n_u8(p, 4), vdup_n_u8(0x0f));
    uint8x16_t table = vld1q_u8(pqc_kernel_mqom_gf16_embed_scalar_table);
    uint8x8_t l = vget_low_u8(vqtbl1q_u8(table, vcombine_u8(lo, vdup_n_u8(0))));
    uint8x8_t h = vget_low_u8(vqtbl1q_u8(table, vcombine_u8(hi, vdup_n_u8(0))));
    const uint8x8_t inter0 = vzip1_u8(l, h);
    const uint8x8_t inter1 = vzip2_u8(l, h);
    vst1_u8(out, inter0); vst1_u8(out + 8, inter1);
}

void pqc_kernel_mqom_gf16_to_gf256_packed_array(uint8_t *out, const uint8_t *packed, size_t n) {
    size_t i = 0;
    for (; i + 16u <= n; i += 16u) mqom_lift_block(out + i, packed + i / 2u);
    if (i < n) {
        /* 尾部补齐：剩余 lane 补零成一个完整向量，只写回 rem 个字节。 */
        const size_t rem = n - i;
        uint8_t pbuf[8], lbuf[16];
        memset(pbuf, 0, sizeof(pbuf));
        memcpy(pbuf, packed + i / 2u, (rem + 1u) / 2u);
        mqom_lift_block(lbuf, pbuf);
        memcpy(out + i, lbuf, rem);
    }
}

void pqc_kernel_mqom_gf16_gf256_const_mul_array(uint8_t *out, uint8_t pqc_kernel_gf16_const,
                                     const uint8_t *gf256, size_t n) {
    pqc_kernel_gf256_mul_const_neon_array(out, gf256, pqc_kernel_mqom_gf16_embed_gf256(pqc_kernel_gf16_const), n);
}

void pqc_kernel_mqom_gf16_gf256_mul_packed_array(uint8_t *out, const uint8_t *packed,
                                      const uint8_t *gf256, size_t n) {
    uint8_t lifted[16];
    size_t i = 0;
    for (; i + 16u <= n; i += 16u) {
        mqom_lift_block(lifted, packed + i / 2u);
        pqc_kernel_gf256_mul_neon_array(out + i, lifted, gf256 + i, 16u);
    }
    if (i < n) {
        const size_t rem = n - i;
        uint8_t pbuf[8], bbuf[16], obuf[16];
        memset(pbuf, 0, sizeof(pbuf));
        memset(bbuf, 0, sizeof(bbuf));
        memcpy(pbuf, packed + i / 2u, (rem + 1u) / 2u);
        memcpy(bbuf, gf256 + i, rem);
        mqom_lift_block(lifted, pbuf);
        pqc_kernel_gf256_mul_neon_array(obuf, lifted, bbuf, 16u);
        memcpy(out + i, obuf, rem);
    }
}

void pqc_kernel_mqom_gf16_gf256_mulacc_packed_array(uint8_t *acc, const uint8_t *packed,
                                         const uint8_t *gf256, size_t n) {
    uint8_t lifted[16];
    size_t i = 0;
    for (; i + 16u <= n; i += 16u) {
        mqom_lift_block(lifted, packed + i / 2u);
        pqc_kernel_gf256_mulacc_neon_array(acc + i, lifted, gf256 + i, 16u);
    }
    if (i < n) {
        const size_t rem = n - i;
        uint8_t pbuf[8], bbuf[16], cbuf[16];
        memset(pbuf, 0, sizeof(pbuf));
        memset(bbuf, 0, sizeof(bbuf));
        memset(cbuf, 0, sizeof(cbuf));
        memcpy(pbuf, packed + i / 2u, (rem + 1u) / 2u);
        memcpy(bbuf, gf256 + i, rem);
        memcpy(cbuf, acc + i, rem);
        mqom_lift_block(lifted, pbuf);
        pqc_kernel_gf256_mulacc_neon_array(cbuf, lifted, bbuf, 16u);
        memcpy(acc + i, cbuf, rem);
    }
}

uint8_t pqc_kernel_mqom_gf16_gf256_dot_packed(const uint8_t *packed, const uint8_t *gf256, size_t n) {
    uint8_t lifted[16];
    uint8_t result = 0;
    size_t i = 0;
    for (; i + 16u <= n; i += 16u) {
        mqom_lift_block(lifted, packed + i / 2u);
        result ^= pqc_kernel_gf256_dot_neon_array(lifted, gf256 + i, 16u);
    }
    if (i < n) {
        const size_t rem = n - i;
        uint8_t pbuf[8], bbuf[16];
        memset(pbuf, 0, sizeof(pbuf));
        memset(bbuf, 0, sizeof(bbuf));
        memcpy(pbuf, packed + i / 2u, (rem + 1u) / 2u);
        memcpy(bbuf, gf256 + i, rem);
        mqom_lift_block(lifted, pbuf);
        result ^= pqc_kernel_gf256_dot_neon_array(lifted, bbuf, 16u);
    }
    return result;
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

void pqc_kernel_gf256_xor_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i)));
    }
    for (; i < n; ++i) out[i] = (uint8_t)(a[i] ^ b[i]);
}
