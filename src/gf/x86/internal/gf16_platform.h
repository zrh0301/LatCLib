#ifndef GF16_PLATFORM_H
#define GF16_PLATFORM_H

#include <immintrin.h>
#include <stdint.h>

#define GF16_BACKEND_NAME "x86-avx2"
#define GF16_VEC_BYTES 32

typedef __m256i pqc_kernel_gf16_vec;

static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_load(const void *p) { return _mm256_loadu_si256((const __m256i *)p); }
static inline void pqc_kernel_gf16_vec_store(void *p, pqc_kernel_gf16_vec v) { _mm256_storeu_si256((__m256i *)p, v); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_zero(void) { return _mm256_setzero_si256(); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_set1(uint8_t x) { return _mm256_set1_epi8((char)x); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_xor_si256(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_and_si256(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_or_si256(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_cmpeq_epi8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_gt_u8_small(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_cmpgt_epi8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_vec tab, pqc_kernel_gf16_vec idx) { return _mm256_shuffle_epi8(tab, idx); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_srli4(pqc_kernel_gf16_vec a) { return _mm256_and_si256(_mm256_srli_epi16(a, 4), _mm256_set1_epi8(0x0f)); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_slli4(pqc_kernel_gf16_vec a) { return _mm256_slli_epi16(a, 4); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_add_u8(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_add_epi8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_sub_u8(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return _mm256_sub_epi8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec v, uint8_t m) { return _mm256_and_si256(v, _mm256_set1_epi8((char)m)); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_not(pqc_kernel_gf16_vec v) { return _mm256_xor_si256(v, _mm256_set1_epi8((char)0xff)); }

#endif
