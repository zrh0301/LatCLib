#include "gf16.h"
#include "gfni_field_map.h"

#include <stdlib.h>
#include <string.h>

/*
 * GF(16) = GF(2)[x]/(x^4+x+1) is embedded as a subfield of the byte field
 * used by GFNI. The constants are the same embedding used by the MAYO x86
 * experiment: multiplication is performed by gf2p8mul and decoded by the
 * inverse affine map below.
 */
#if defined(GF16_ENABLE_GFNI)

#if defined(__GNUC__) || defined(__clang__)
#include <cpuid.h>
#elif defined(_MSC_VER)
#include <intrin.h>
#endif

static inline pqc_kernel_gf16_vec pqc_kernel_gf16_gfni_embed(pqc_kernel_gf16_vec value) {
    const pqc_kernel_gf16_vec table = _mm256_load_si256((const __m256i *)gfni_gf16_embed_table);
    return _mm256_shuffle_epi8(table, _mm256_and_si256(value, _mm256_set1_epi8(0x0f)));
}

static inline pqc_kernel_gf16_vec pqc_kernel_gf16_gfni_decode(pqc_kernel_gf16_vec value) {
    const pqc_kernel_gf16_vec inverse = _mm256_set1_epi64x(0x0104201400000000ULL);
    return _mm256_and_si256(
        _mm256_gf2p8affine_epi64_epi8(value, inverse, 0),
        _mm256_set1_epi8(0x0f));
}

static inline pqc_kernel_gf16_vec pqc_kernel_gf16_gfni_mul_embedded(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) {
    return pqc_kernel_gf16_gfni_decode(_mm256_gf2p8mul_epi8(a, b));
}

static inline uint8_t pqc_kernel_gf16_gfni_embedded_scalar(uint8_t value) {
    return gfni_gf16_embed_scalar_table[value & 0x0f];
}

#endif

int pqc_kernel_gf16_gfni_available(void) {
#if defined(GF16_ENABLE_GFNI)
#if defined(__GNUC__) || defined(__clang__)
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        return (ecx & (1u << 8)) != 0;
    }
    return 0;
#elif defined(_MSC_VER)
    int cpu_info[4] = {0, 0, 0, 0};
    __cpuidex(cpu_info, 7, 0);
    return (cpu_info[2] & (1 << 8)) != 0;
#else
    return 1;
#endif
#else
    return 0;
#endif
}

pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_gfni(pqc_kernel_gf16_vec a, uint8_t k) {
#if defined(GF16_ENABLE_GFNI)
    const pqc_kernel_gf16_vec embedded_a = pqc_kernel_gf16_gfni_embed(a);
    const pqc_kernel_gf16_vec embedded_k = _mm256_set1_epi8((char)pqc_kernel_gf16_gfni_embedded_scalar(k));
    return pqc_kernel_gf16_gfni_mul_embedded(embedded_a, embedded_k);
#else
    return pqc_kernel_gf16_vec_mul_const_shuffle(a, k);
#endif
}

pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_gfni(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) {
#if defined(GF16_ENABLE_GFNI)
    return pqc_kernel_gf16_gfni_mul_embedded(pqc_kernel_gf16_gfni_embed(a), pqc_kernel_gf16_gfni_embed(b));
#else
    return pqc_kernel_gf16_vec_mul_bitplane(a, b);
#endif
}

void pqc_kernel_gf16_gfni_build_multab(pqc_kernel_gf16_vec *table, uint8_t k) {
    if (!table) {
        return;
    }
#if defined(GF16_ENABLE_GFNI)
    const pqc_kernel_gf16_vec domain = _mm256_load_si256((const __m256i *)gfni_gf16_domain_table);
    const pqc_kernel_gf16_vec embedded_k = _mm256_set1_epi8((char)pqc_kernel_gf16_gfni_embedded_scalar(k));
    *table = pqc_kernel_gf16_gfni_mul_embedded(pqc_kernel_gf16_gfni_embed(domain), embedded_k);
#else
    *table = pqc_kernel_gf16_mtk_uniform[k & 0x0f];
#endif
}

pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_gfni_table(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec table) {
    return pqc_kernel_gf16_vec_shuffle(table, pqc_kernel_gf16_vec_mask_u8(a, 0x0f));
}

void pqc_kernel_gf16_mul_const_gfni_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    size_t i = 0;
#if defined(GF16_ENABLE_GFNI)
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_const_gfni(pqc_kernel_gf16_vec_load(a + i), k);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
#else
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_const_shuffle(pqc_kernel_gf16_vec_load(a + i), k);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
#endif
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_table[a[i] & 0x0f][k & 0x0f];
    }
}

void pqc_kernel_gf16_mul_gfni_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
#if defined(GF16_ENABLE_GFNI)
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_gfni(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i));
        pqc_kernel_gf16_vec_store(out + i, r);
    }
#else
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i));
        pqc_kernel_gf16_vec_store(out + i, r);
    }
#endif
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_table[a[i] & 0x0f][b[i] & 0x0f];
    }
}

void pqc_kernel_gf16_mul_const_gfni_table_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    pqc_kernel_gf16_vec table;
    pqc_kernel_gf16_gfni_build_multab(&table, k);
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_const_gfni_table(pqc_kernel_gf16_vec_load(a + i), table);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_table[a[i] & 0x0f][k & 0x0f];
    }
}

void pqc_kernel_gf16_snova_ps_gfni(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                        size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc,
                    pqc_kernel_gf16_vec_mul_const_gfni(pqc_kernel_gf16_vec_load(sig + c * lane_count + i),
                                             coeff[r * cols + c]));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[sig[c * lane_count + i] & 0x0f]
                                      [coeff[r * cols + c] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

void pqc_kernel_gf16_snova_ps_gfni_table(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                              size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    pqc_kernel_gf16_vec tables[16];
    for (unsigned k = 0; k < 16; ++k) {
        pqc_kernel_gf16_gfni_build_multab(&tables[k], (uint8_t)k);
    }
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc,
                    pqc_kernel_gf16_vec_mul_const_gfni_table(
                        pqc_kernel_gf16_vec_load(sig + c * lane_count + i),
                        tables[coeff[r * cols + c] & 0x0f]));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[sig[c * lane_count + i] & 0x0f]
                                      [coeff[r * cols + c] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

void pqc_kernel_gf16_snova_sts_gfni(uint8_t *out, const uint8_t *left, const uint8_t *right,
                         size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
#if defined(GF16_ENABLE_GFNI)
    pqc_kernel_gf16_vec *acc = (pqc_kernel_gf16_vec *)malloc(rows * sizeof(*acc));
    if (!acc) {
        /* Keep the public API usable under allocation pressure. */
        pqc_kernel_gf16_snova_sts_bitplane(out, left, right, rows, cols, lane_count);
        return;
    }
    for (size_t i = 0; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
        for (size_t r = 0; r < rows; ++r) {
            acc[r] = pqc_kernel_gf16_vec_zero();
        }
        for (size_t c = 0; c < cols; ++c) {
            const pqc_kernel_gf16_vec right_v = pqc_kernel_gf16_gfni_embed(
                pqc_kernel_gf16_vec_load(right + c * lane_count + i));
            for (size_t r = 0; r < rows; ++r) {
                const pqc_kernel_gf16_vec left_v = pqc_kernel_gf16_gfni_embed(
                    pqc_kernel_gf16_vec_load(left + (r * cols + c) * lane_count + i));
                acc[r] = pqc_kernel_gf16_vec_xor(acc[r], pqc_kernel_gf16_gfni_mul_embedded(left_v, right_v));
            }
        }
        for (size_t r = 0; r < rows; ++r) {
            pqc_kernel_gf16_vec_store(out + r * lane_count + i, acc[r]);
        }
    }
    free(acc);
#else
    pqc_kernel_gf16_snova_sts_bitplane(out, left, right, rows, cols, lane_count);
#endif
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        for (size_t i = lane_count & ~(size_t)(GF16_VEC_BYTES - 1); i < lane_count; ++i) {
            uint8_t v = 0;
            for (size_t c = 0; c < cols; ++c) {
                v ^= pqc_kernel_gf16_mul_table[left[(r * cols + c) * lane_count + i] & 0x0f]
                                     [right[c * lane_count + i] & 0x0f];
            }
            dst[i] = v;
        }
    }
}
