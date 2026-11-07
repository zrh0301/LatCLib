#include "gfext.h"
#include "gfni_field_map.h"

#include <string.h>

uint8_t pqc_kernel_gf16_embed_gf256(uint8_t a) {
    return gfni_gf16_embed_scalar_table[a & 0x0f];
}

static pqc_kernel_gf256_vec pqc_kernel_gf16_embed_vec(pqc_kernel_gf256_vec value) {
    const pqc_kernel_gf256_vec table = _mm256_load_si256((const __m256i *)gfni_gf16_embed_table);
    return pqc_kernel_gf16_vec_shuffle(table, pqc_kernel_gf16_vec_mask_u8(value, 0x0f));
}

static pqc_kernel_gf256_vec pqc_kernel_gf16_gf256_mul_vec(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b, int use_gfni) {
#if defined(GF16_ENABLE_GFNI)
    if (use_gfni) {
        return _mm256_gf2p8mul_epi8(a, b);
    }
#else
    (void)use_gfni;
#endif
    return pqc_kernel_gf256_vec_mul_avx2(a, b);
}

/* The MQOM flow calls the packed dot/mulacc primitive once per matrix row.
 * Cache the feature probe so CPUID is paid once at process start, matching the
 * MQOM implementation's compile/runtime dispatch boundary. */
static int pqc_kernel_mqom_gfni_cached = -1;

static int mqom_use_gfni(void) {
    if (pqc_kernel_mqom_gfni_cached < 0) {
#if defined(GF16_ENABLE_GFNI)
        pqc_kernel_mqom_gfni_cached = pqc_kernel_gf256_gfni_available();
#else
        pqc_kernel_mqom_gfni_cached = 0;
#endif
    }
    return pqc_kernel_mqom_gfni_cached;
}

void pqc_kernel_gf16_to_gf256_array(uint8_t *out, const uint8_t *in, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(in + i)));
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_embed_gf256(in[i]);
    }
}

void pqc_kernel_gf16_gf256_mul_array(uint8_t *out, const uint8_t *gf16,
                          const uint8_t *gf256, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec a = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i));
        const pqc_kernel_gf256_vec b = pqc_kernel_gf16_vec_load(gf256 + i);
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf256_vec_mul_avx2(a, b));
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf256_mul_ref(pqc_kernel_gf16_embed_gf256(gf16[i]), gf256[i]);
    }
}

void pqc_kernel_gf16_gf256_mulacc_array(uint8_t *acc, const uint8_t *gf16,
                             const uint8_t *gf256, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec a = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i));
        const pqc_kernel_gf256_vec b = pqc_kernel_gf16_vec_load(gf256 + i);
        const pqc_kernel_gf256_vec c = pqc_kernel_gf16_vec_load(acc + i);
        pqc_kernel_gf16_vec_store(acc + i, _mm256_xor_si256(c, pqc_kernel_gf16_gf256_mul_vec(a, b, 0)));
    }
    for (; i < n; ++i) {
        acc[i] ^= pqc_kernel_gf256_mul_ref(pqc_kernel_gf16_embed_gf256(gf16[i]), gf256[i]);
    }
}

void pqc_kernel_gf16_gf256_mulacc_gfni_array(uint8_t *acc, const uint8_t *gf16,
                                  const uint8_t *gf256, size_t n) {
    size_t i = 0;
#if defined(GF16_ENABLE_GFNI)
    const int use_gfni = pqc_kernel_gf256_gfni_available();
#else
    const int use_gfni = 0;
#endif
    /* Four independent accumulators hide the GFNI/shuffle latency on long rows. */
    for (; i + 4u * GF16_VEC_BYTES <= n; i += 4u * GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec a0 = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i));
        const pqc_kernel_gf256_vec a1 = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i + GF16_VEC_BYTES));
        const pqc_kernel_gf256_vec a2 = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i + 2u * GF16_VEC_BYTES));
        const pqc_kernel_gf256_vec a3 = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i + 3u * GF16_VEC_BYTES));
        const pqc_kernel_gf256_vec b0 = pqc_kernel_gf16_vec_load(gf256 + i);
        const pqc_kernel_gf256_vec b1 = pqc_kernel_gf16_vec_load(gf256 + i + GF16_VEC_BYTES);
        const pqc_kernel_gf256_vec b2 = pqc_kernel_gf16_vec_load(gf256 + i + 2u * GF16_VEC_BYTES);
        const pqc_kernel_gf256_vec b3 = pqc_kernel_gf16_vec_load(gf256 + i + 3u * GF16_VEC_BYTES);
        const pqc_kernel_gf256_vec c0 = pqc_kernel_gf16_vec_load(acc + i);
        const pqc_kernel_gf256_vec c1 = pqc_kernel_gf16_vec_load(acc + i + GF16_VEC_BYTES);
        const pqc_kernel_gf256_vec c2 = pqc_kernel_gf16_vec_load(acc + i + 2u * GF16_VEC_BYTES);
        const pqc_kernel_gf256_vec c3 = pqc_kernel_gf16_vec_load(acc + i + 3u * GF16_VEC_BYTES);
        pqc_kernel_gf16_vec_store(acc + i, _mm256_xor_si256(c0, pqc_kernel_gf16_gf256_mul_vec(a0, b0, use_gfni)));
        pqc_kernel_gf16_vec_store(acc + i + GF16_VEC_BYTES,
                       _mm256_xor_si256(c1, pqc_kernel_gf16_gf256_mul_vec(a1, b1, use_gfni)));
        pqc_kernel_gf16_vec_store(acc + i + 2u * GF16_VEC_BYTES,
                       _mm256_xor_si256(c2, pqc_kernel_gf16_gf256_mul_vec(a2, b2, use_gfni)));
        pqc_kernel_gf16_vec_store(acc + i + 3u * GF16_VEC_BYTES,
                       _mm256_xor_si256(c3, pqc_kernel_gf16_gf256_mul_vec(a3, b3, use_gfni)));
    }
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec a = pqc_kernel_gf16_embed_vec(pqc_kernel_gf16_vec_load(gf16 + i));
        const pqc_kernel_gf256_vec b = pqc_kernel_gf16_vec_load(gf256 + i);
        const pqc_kernel_gf256_vec c = pqc_kernel_gf16_vec_load(acc + i);
        const pqc_kernel_gf256_vec product = pqc_kernel_gf16_gf256_mul_vec(a, b, use_gfni);
        pqc_kernel_gf16_vec_store(acc + i, _mm256_xor_si256(c, product));
    }
    for (; i < n; ++i) {
        acc[i] ^= pqc_kernel_gf256_mul_ref(pqc_kernel_gf16_embed_gf256(gf16[i]), gf256[i]);
    }
}

/* === MQOM v2.1.1 packed GF(16)-GF(256) compatibility path ===
 * MQOM stores two GF(16) coefficients per byte (low nibble first).  Its
 * AVX2 lift expands 16 packed bytes into 32 GF(256) lanes with two shuffles:
 * one separates the nibbles and one applies the subfield map. */
static inline pqc_kernel_gf256_vec mqom_lift_packed_vec(const uint8_t *packed,
                                             size_t packed_bytes) {
    const pqc_kernel_gf256_vec nib_mask = _mm256_set1_epi8(0x0f);
    const pqc_kernel_gf256_vec shuff_even = _mm256_set_epi8(
        -1, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10, -1, 9, -1, 8,
        -1, 7, -1, 6, -1, 5, -1, 4, -1, 3, -1, 2, -1, 1, -1, 0);
    const pqc_kernel_gf256_vec shuff_odd = _mm256_set_epi8(
        15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10, -1, 9, -1, 8, -1,
        7, -1, 6, -1, 5, -1, 4, -1, 3, -1, 2, -1, 1, -1, 0, -1);
    const pqc_kernel_gf256_vec lift_table = _mm256_load_si256(
        (const __m256i *)pqc_kernel_mqom_gf16_embed_table);
    uint8_t local[16] = {0};
    if (packed_bytes > sizeof(local)) packed_bytes = sizeof(local);
    memcpy(local, packed, packed_bytes);
    const __m128i low = _mm_loadu_si128((const __m128i *)local);
    const pqc_kernel_gf256_vec packed_vec = _mm256_broadcastsi128_si256(low);
    const pqc_kernel_gf256_vec lo = _mm256_and_si256(packed_vec, nib_mask);
    const pqc_kernel_gf256_vec hi = _mm256_and_si256(_mm256_srli_epi64(packed_vec, 4), nib_mask);
    const pqc_kernel_gf256_vec separated = _mm256_or_si256(
        _mm256_shuffle_epi8(lo, shuff_even),
        _mm256_shuffle_epi8(hi, shuff_odd));
    return _mm256_shuffle_epi8(lift_table, separated);
}

static inline pqc_kernel_gf256_vec pqc_kernel_mqom_gf16_gf256_mul_vec(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b,
                                                int use_gfni) {
#if defined(GF16_ENABLE_GFNI)
    if (use_gfni) return _mm256_gf2p8mul_epi8(a, b);
#else
    (void)use_gfni;
#endif
    return pqc_kernel_gf256_vec_mul_avx2(a, b);
}

static uint8_t mqom_reduce_vec(pqc_kernel_gf256_vec v) {
    uint64_t x = (uint64_t)_mm256_extract_epi64(v, 0) ^
                 (uint64_t)_mm256_extract_epi64(v, 1) ^
                 (uint64_t)_mm256_extract_epi64(v, 2) ^
                 (uint64_t)_mm256_extract_epi64(v, 3);
    x ^= x >> 32;
    x ^= x >> 16;
    x ^= x >> 8;
    return (uint8_t)x;
}

uint8_t pqc_kernel_mqom_gf16_embed_gf256(uint8_t a) {
    return pqc_kernel_mqom_gf16_embed_scalar_table[a & 0x0f];
}

void pqc_kernel_mqom_gf16_to_gf256_packed_array(uint8_t *out, const uint8_t *packed,
                                     size_t n) {
    size_t i = 0;
    for (; i + 32u <= n; i += 32u) {
        pqc_kernel_gf16_vec_store(out + i, mqom_lift_packed_vec(packed + i / 2u, 16));
    }
    if (i < n) {
        const size_t rem = n - i;
        pqc_kernel_gf16_vec lifted = mqom_lift_packed_vec(packed + i / 2u, (rem + 1u) / 2u);
        uint8_t local[32];
        pqc_kernel_gf16_vec_store(local, lifted);
        memcpy(out + i, local, rem);
    }
}

void pqc_kernel_mqom_gf16_gf256_const_mul_array(uint8_t *out, uint8_t pqc_kernel_gf16_const,
                                     const uint8_t *gf256, size_t n) {
    pqc_kernel_gf256_mul_const_gfni_array(out, gf256,
                               pqc_kernel_mqom_gf16_embed_gf256(pqc_kernel_gf16_const), n);
}

void pqc_kernel_mqom_gf16_gf256_mul_packed_array(uint8_t *out,
                                      const uint8_t *pqc_kernel_gf16_packed,
                                      const uint8_t *gf256, size_t n) {
    size_t i = 0;
    const int use_gfni = mqom_use_gfni();
    for (; i + 32u <= n; i += 32u) {
        const pqc_kernel_gf256_vec a = mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u, 16);
        const pqc_kernel_gf256_vec b = pqc_kernel_gf16_vec_load(gf256 + i);
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_mqom_gf16_gf256_mul_vec(a, b, use_gfni));
    }
    for (; i < n; ++i) {
        const uint8_t nibble = (uint8_t)((pqc_kernel_gf16_packed[i / 2u] >> (4u * (i & 1u))) & 0x0f);
        out[i] = pqc_kernel_gf256_mul_ref(pqc_kernel_mqom_gf16_embed_gf256(nibble), gf256[i]);
    }
}

void pqc_kernel_mqom_gf16_gf256_mulacc_packed_array(uint8_t *acc,
                                         const uint8_t *pqc_kernel_gf16_packed,
                                         const uint8_t *gf256, size_t n) {
    size_t i = 0;
    const int use_gfni = mqom_use_gfni();
    for (; i + 4u * 32u <= n; i += 4u * 32u) {
        const pqc_kernel_gf256_vec a0 = mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u, 16);
        const pqc_kernel_gf256_vec a1 = mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u + 16u, 16);
        const pqc_kernel_gf256_vec a2 = mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u + 32u, 16);
        const pqc_kernel_gf256_vec a3 = mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u + 48u, 16);
        const pqc_kernel_gf256_vec b0 = pqc_kernel_gf16_vec_load(gf256 + i);
        const pqc_kernel_gf256_vec b1 = pqc_kernel_gf16_vec_load(gf256 + i + 32u);
        const pqc_kernel_gf256_vec b2 = pqc_kernel_gf16_vec_load(gf256 + i + 64u);
        const pqc_kernel_gf256_vec b3 = pqc_kernel_gf16_vec_load(gf256 + i + 96u);
        const pqc_kernel_gf256_vec c0 = pqc_kernel_gf16_vec_load(acc + i);
        const pqc_kernel_gf256_vec c1 = pqc_kernel_gf16_vec_load(acc + i + 32u);
        const pqc_kernel_gf256_vec c2 = pqc_kernel_gf16_vec_load(acc + i + 64u);
        const pqc_kernel_gf256_vec c3 = pqc_kernel_gf16_vec_load(acc + i + 96u);
        pqc_kernel_gf16_vec_store(acc + i, _mm256_xor_si256(c0, pqc_kernel_mqom_gf16_gf256_mul_vec(a0, b0, use_gfni)));
        pqc_kernel_gf16_vec_store(acc + i + 32u, _mm256_xor_si256(c1, pqc_kernel_mqom_gf16_gf256_mul_vec(a1, b1, use_gfni)));
        pqc_kernel_gf16_vec_store(acc + i + 64u, _mm256_xor_si256(c2, pqc_kernel_mqom_gf16_gf256_mul_vec(a2, b2, use_gfni)));
        pqc_kernel_gf16_vec_store(acc + i + 96u, _mm256_xor_si256(c3, pqc_kernel_mqom_gf16_gf256_mul_vec(a3, b3, use_gfni)));
    }
    for (; i + 32u <= n; i += 32u) {
        const pqc_kernel_gf256_vec a = mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u, 16);
        const pqc_kernel_gf256_vec b = pqc_kernel_gf16_vec_load(gf256 + i);
        const pqc_kernel_gf256_vec c = pqc_kernel_gf16_vec_load(acc + i);
        pqc_kernel_gf16_vec_store(acc + i, _mm256_xor_si256(c, pqc_kernel_mqom_gf16_gf256_mul_vec(a, b, use_gfni)));
    }
    if (i < n) {
        /* 尾部补齐：n=116 这类长度剩余的 20 个字节若走标量，会成为固定开销。
         * 把剩余 lane 补零成一个完整向量，只写回 rem 个字节。 */
        const size_t rem = n - i;
        uint8_t pbuf[16], bbuf[32], cbuf[32];
        memset(pbuf, 0, sizeof(pbuf));
        memset(bbuf, 0, sizeof(bbuf));
        memset(cbuf, 0, sizeof(cbuf));
        memcpy(pbuf, pqc_kernel_gf16_packed + i / 2u, (rem + 1u) / 2u);
        memcpy(bbuf, gf256 + i, rem);
        memcpy(cbuf, acc + i, rem);
        pqc_kernel_gf16_vec_store(cbuf, _mm256_xor_si256(
            pqc_kernel_gf16_vec_load(cbuf),
            pqc_kernel_mqom_gf16_gf256_mul_vec(mqom_lift_packed_vec(pbuf, 16),
                                    pqc_kernel_gf16_vec_load(bbuf), use_gfni)));
        memcpy(acc + i, cbuf, rem);
    }
}

uint8_t pqc_kernel_mqom_gf16_gf256_dot_packed(const uint8_t *pqc_kernel_gf16_packed,
                                   const uint8_t *gf256, size_t n) {
    size_t i = 0;
    pqc_kernel_gf256_vec acc0 = _mm256_setzero_si256();
    pqc_kernel_gf256_vec acc1 = _mm256_setzero_si256();
    pqc_kernel_gf256_vec acc2 = _mm256_setzero_si256();
    pqc_kernel_gf256_vec acc3 = _mm256_setzero_si256();
    const int use_gfni = mqom_use_gfni();
    for (; i + 4u * 32u <= n; i += 4u * 32u) {
        acc0 = _mm256_xor_si256(acc0, pqc_kernel_mqom_gf16_gf256_mul_vec(
            mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u, 16),
            pqc_kernel_gf16_vec_load(gf256 + i), use_gfni));
        acc1 = _mm256_xor_si256(acc1, pqc_kernel_mqom_gf16_gf256_mul_vec(
            mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u + 16u, 16),
            pqc_kernel_gf16_vec_load(gf256 + i + 32u), use_gfni));
        acc2 = _mm256_xor_si256(acc2, pqc_kernel_mqom_gf16_gf256_mul_vec(
            mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u + 32u, 16),
            pqc_kernel_gf16_vec_load(gf256 + i + 64u), use_gfni));
        acc3 = _mm256_xor_si256(acc3, pqc_kernel_mqom_gf16_gf256_mul_vec(
            mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u + 48u, 16),
            pqc_kernel_gf16_vec_load(gf256 + i + 96u), use_gfni));
    }
    for (; i + 32u <= n; i += 32u) {
        acc0 = _mm256_xor_si256(acc0, pqc_kernel_mqom_gf16_gf256_mul_vec(
            mqom_lift_packed_vec(pqc_kernel_gf16_packed + i / 2u, 16),
            pqc_kernel_gf16_vec_load(gf256 + i), use_gfni));
    }
    if (i < n) {
        /* 尾部补齐，做法与 mulacc 一致。 */
        const size_t rem = n - i;
        uint8_t pbuf[16], bbuf[32];
        memset(pbuf, 0, sizeof(pbuf));
        memset(bbuf, 0, sizeof(bbuf));
        memcpy(pbuf, pqc_kernel_gf16_packed + i / 2u, (rem + 1u) / 2u);
        memcpy(bbuf, gf256 + i, rem);
        acc0 = _mm256_xor_si256(acc0, pqc_kernel_mqom_gf16_gf256_mul_vec(
            mqom_lift_packed_vec(pbuf, 16), pqc_kernel_gf16_vec_load(bbuf), use_gfni));
    }
    return mqom_reduce_vec(_mm256_xor_si256(
        _mm256_xor_si256(acc0, acc1), _mm256_xor_si256(acc2, acc3)));
}
