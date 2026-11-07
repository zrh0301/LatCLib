/* SPDX-License-Identifier: MIT
 * GF(16), x^4+x+1: linear reduction as GFNI binary affine maps.
 * Matrices reproduce the original maps for every input bit pattern.
 * No secret-indexed memory access and no key-dependent retained state.
 */
#ifndef ROUND3_GF_H
#define ROUND3_GF_H
#if defined(__GFNI__) && defined(__AVX2__)
#include <immintrin.h>
static inline __m256i round3_reduce256(__m256i a) {
    return _mm256_gf2p8affine_epi64_epi8(a,
        _mm256_set1_epi64x(0x113264c800000000ULL), 0);
}
static inline __m128i round3_reduce128(__m128i a) {
    return _mm_gf2p8affine_epi64_epi8(a,
        _mm_set1_epi64x(0x113264c800000000ULL), 0);
}
static inline __m256i round3_compress256(__m256i a) {
    __m256i lo = _mm256_gf2p8affine_epi64_epi8(a,
        _mm256_set1_epi64x(0x8192244800000000ULL), 0);
    __m256i hi = _mm256_gf2p8affine_epi64_epi8(a,
        _mm256_set1_epi64x(0x246cc99200000000ULL), 0);
    return _mm256_and_si256(_mm256_xor_si256(lo, _mm256_srli_epi16(hi, 8)),
        _mm256_set1_epi16(15));
}
static inline __m128i round3_compress128(__m128i a) {
    __m128i lo = _mm_gf2p8affine_epi64_epi8(a,
        _mm_set1_epi64x(0x8192244800000000ULL), 0);
    __m128i hi = _mm_gf2p8affine_epi64_epi8(a,
        _mm_set1_epi64x(0x246cc99200000000ULL), 0);
    return _mm_and_si128(_mm_xor_si128(lo, _mm_srli_epi16(hi, 8)),
        _mm_set1_epi16(15));
}
/* Secret coefficient -> register table, shared by every vector in the row. */
static inline __m256i round3_scalar_table(uint8_t s) {
    const __m256i x = _mm256_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                                      0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    return round3_reduce256(_mm256_gf2p8mul_epi8(_mm256_set1_epi8((char)s),x));
}
static inline void round3_row_scale(uint8_t *row,uint8_t s,int k0,int kend) {
    __m256i t=round3_scalar_table(s);
    for(int k=k0;k<kend;k+=32) {
        __m256i v=_mm256_loadu_si256((const __m256i*)(row+k));
        _mm256_storeu_si256((__m256i*)(row+k),_mm256_shuffle_epi8(t,v));
    }
}
static inline void round3_row_axpy(uint8_t *dst,const uint8_t *src,uint8_t s,int k0,int kend) {
    __m256i t=round3_scalar_table(s);
    for(int k=k0;k<kend;k+=32) {
        __m256i v=_mm256_loadu_si256((const __m256i*)(src+k));
        __m256i d=_mm256_loadu_si256((const __m256i*)(dst+k));
        _mm256_storeu_si256((__m256i*)(dst+k),_mm256_xor_si256(d,_mm256_shuffle_epi8(t,v)));
    }
}
#endif
#endif
