#include "gfext.h"

/*
 * AVX2-only fallback. Each byte lane is shifted independently and reduced by
 * 0x1b, so this path is constant-time and does not require GFNI.
 */
pqc_kernel_gf256_vec pqc_kernel_gf256_vec_mul_avx2(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b) {
    pqc_kernel_gf256_vec result = _mm256_setzero_si256();
    const pqc_kernel_gf256_vec zero = _mm256_setzero_si256();
    const pqc_kernel_gf256_vec all_ones = _mm256_set1_epi8((char)0xff);
    const pqc_kernel_gf256_vec reduction = _mm256_set1_epi8((char)0x1b);

    for (unsigned i = 0; i < 8; ++i) {
        const pqc_kernel_gf256_vec b_lsb = _mm256_and_si256(b, _mm256_set1_epi8(1));
        const pqc_kernel_gf256_vec take = _mm256_xor_si256(
            _mm256_cmpeq_epi8(b_lsb, zero), all_ones);
        result = _mm256_xor_si256(result, _mm256_and_si256(a, take));

        const pqc_kernel_gf256_vec high = _mm256_and_si256(a, _mm256_set1_epi8((char)0x80));
        const pqc_kernel_gf256_vec reduce = _mm256_xor_si256(
            _mm256_cmpeq_epi8(high, zero), all_ones);
        a = _mm256_and_si256(_mm256_slli_epi16(a, 1),
                             _mm256_set1_epi8((char)0xfe));
        a = _mm256_xor_si256(a, _mm256_and_si256(reduction, reduce));
        b = _mm256_and_si256(_mm256_srli_epi16(b, 1),
                             _mm256_set1_epi8(0x7f));
    }
    return result;
}

int pqc_kernel_gf256_gfni_available(void) {
#if defined(GF16_ENABLE_GFNI)
    return pqc_kernel_gf16_gfni_available();
#else
    return 0;
#endif
}

pqc_kernel_gf256_vec pqc_kernel_gf256_vec_mul_gfni(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b) {
#if defined(GF16_ENABLE_GFNI)
    return _mm256_gf2p8mul_epi8(a, b);
#else
    return pqc_kernel_gf256_vec_mul_avx2(a, b);
#endif
}

void pqc_kernel_gf256_mul_gfni_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
#if defined(GF16_ENABLE_GFNI)
    const int use_gfni = pqc_kernel_gf256_gfni_available();
#else
    const int use_gfni = 0;
#endif
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec av = pqc_kernel_gf16_vec_load(a + i);
        const pqc_kernel_gf256_vec bv = pqc_kernel_gf16_vec_load(b + i);
        const pqc_kernel_gf256_vec product = use_gfni
            ? pqc_kernel_gf256_vec_mul_gfni(av, bv)
            : pqc_kernel_gf256_vec_mul_avx2(av, bv);
        pqc_kernel_gf16_vec_store(out + i, product);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf256_mul_ref(a[i], b[i]);
    }
}

void pqc_kernel_gf256_mulacc_gfni_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
#if defined(GF16_ENABLE_GFNI)
    const int use_gfni = pqc_kernel_gf256_gfni_available();
#else
    const int use_gfni = 0;
#endif
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec av = pqc_kernel_gf16_vec_load(a + i);
        const pqc_kernel_gf256_vec bv = pqc_kernel_gf16_vec_load(b + i);
        const pqc_kernel_gf256_vec cv = pqc_kernel_gf16_vec_load(acc + i);
        const pqc_kernel_gf256_vec product = use_gfni
            ? pqc_kernel_gf256_vec_mul_gfni(av, bv)
            : pqc_kernel_gf256_vec_mul_avx2(av, bv);
        pqc_kernel_gf16_vec_store(acc + i, _mm256_xor_si256(cv, product));
    }
    for (; i < n; ++i) {
        acc[i] ^= pqc_kernel_gf256_mul_ref(a[i], b[i]);
    }
}

static inline pqc_kernel_gf256_vec pqc_kernel_gf256_mul_runtime(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b, int use_gfni) {
#if defined(GF16_ENABLE_GFNI)
    if (use_gfni) return _mm256_gf2p8mul_epi8(a, b);
#else
    (void)use_gfni;
#endif
    return pqc_kernel_gf256_vec_mul_avx2(a, b);
}

static uint8_t pqc_kernel_gf256_reduce_vec(pqc_kernel_gf256_vec v) {
    uint64_t x = (uint64_t)_mm256_extract_epi64(v, 0) ^
                 (uint64_t)_mm256_extract_epi64(v, 1) ^
                 (uint64_t)_mm256_extract_epi64(v, 2) ^
                 (uint64_t)_mm256_extract_epi64(v, 3);
    x ^= x >> 32;
    x ^= x >> 16;
    x ^= x >> 8;
    return (uint8_t)x;
}

/* MQOM's GF(256) constant-vector primitive, exposed as a reusable GF_cal
 * operation.  The constant is broadcast once and the tail is scalar-safe. */
void pqc_kernel_gf256_mul_const_gfni_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    size_t i = 0;
#if defined(GF16_ENABLE_GFNI)
    const int use_gfni = pqc_kernel_gf256_gfni_available();
#else
    const int use_gfni = 0;
#endif
    const pqc_kernel_gf256_vec kv = _mm256_set1_epi8((char)k);
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf256_vec av = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf256_mul_runtime(av, kv, use_gfni));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf256_mul_ref(a[i], k);
}

/* GF(256) vector dot product.  Four independent accumulators follow the
 * instruction-level parallelism used by the mixed MQOM path, while retaining
 * a scalar tail for arbitrary lengths. */
uint8_t pqc_kernel_gf256_dot_gfni_array(const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    pqc_kernel_gf256_vec acc0 = _mm256_setzero_si256();
    pqc_kernel_gf256_vec acc1 = _mm256_setzero_si256();
    pqc_kernel_gf256_vec acc2 = _mm256_setzero_si256();
    pqc_kernel_gf256_vec acc3 = _mm256_setzero_si256();
#if defined(GF16_ENABLE_GFNI)
    const int use_gfni = pqc_kernel_gf256_gfni_available();
#else
    const int use_gfni = 0;
#endif
    for (; i + 4u * GF16_VEC_BYTES <= n; i += 4u * GF16_VEC_BYTES) {
        acc0 = _mm256_xor_si256(acc0, pqc_kernel_gf256_mul_runtime(
            pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i), use_gfni));
        acc1 = _mm256_xor_si256(acc1, pqc_kernel_gf256_mul_runtime(
            pqc_kernel_gf16_vec_load(a + i + GF16_VEC_BYTES),
            pqc_kernel_gf16_vec_load(b + i + GF16_VEC_BYTES), use_gfni));
        acc2 = _mm256_xor_si256(acc2, pqc_kernel_gf256_mul_runtime(
            pqc_kernel_gf16_vec_load(a + i + 2u * GF16_VEC_BYTES),
            pqc_kernel_gf16_vec_load(b + i + 2u * GF16_VEC_BYTES), use_gfni));
        acc3 = _mm256_xor_si256(acc3, pqc_kernel_gf256_mul_runtime(
            pqc_kernel_gf16_vec_load(a + i + 3u * GF16_VEC_BYTES),
            pqc_kernel_gf16_vec_load(b + i + 3u * GF16_VEC_BYTES), use_gfni));
    }
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        acc0 = _mm256_xor_si256(acc0, pqc_kernel_gf256_mul_runtime(
            pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i), use_gfni));
    }
    uint8_t result = pqc_kernel_gf256_reduce_vec(_mm256_xor_si256(
        _mm256_xor_si256(acc0, acc1), _mm256_xor_si256(acc2, acc3)));
    for (; i < n; ++i) result ^= pqc_kernel_gf256_mul_ref(a[i], b[i]);
    return result;
}
