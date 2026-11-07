// SPDX-License-Identifier: Apache-2.0

#ifndef MAYO_AVX2_MATVEC_H
#define MAYO_AVX2_MATVEC_H

#include <immintrin.h>
#include <stdint.h>

#define MAYO_AVX2_MATVEC_MAX_CHUNKS 8

static inline uint8_t gf16_mul_scalar_local(uint8_t a, uint8_t b) {
    uint8_t product = 0;
    a &= 0x0f;
    b &= 0x0f;
    for (int i = 0; i < 4; ++i) {
        const uint8_t bit_mask = (uint8_t)(0u - (b & 1u));
        const uint8_t carry_mask = (uint8_t)(0u - ((a >> 3) & 1u));
        product ^= a & bit_mask;
        a = (uint8_t)(((a << 1) & 0x0f) ^ (0x03 & carry_mask));
        b = (uint8_t)(b >> 1);
    }
    return product & 0x0f;
}

/* GF_cal-inspired multiply-by-x for byte lanes containing GF16 values. */
static inline __m256i mayo_gf16_mul_x_vec(__m256i value) {
    const __m256i low_mask = _mm256_set1_epi8(0x0f);
    const __m256i msb_mask = _mm256_set1_epi8(0x08);
    const __m256i reduce = _mm256_set1_epi8(0x03);
    const __m256i carry = _mm256_and_si256(value, msb_mask);
    const __m256i shifted = _mm256_and_si256(_mm256_slli_epi16(value, 1), low_mask);
    const __m256i carry_mask = _mm256_cmpeq_epi8(carry, msb_mask);
    return _mm256_xor_si256(shifted, _mm256_and_si256(carry_mask, reduce));
}

static inline __m128i mayo_gf16_mul_x_vec128(__m128i value) {
    const __m128i low_mask = _mm_set1_epi8(0x0f);
    const __m128i msb_mask = _mm_set1_epi8(0x08);
    const __m128i reduce = _mm_set1_epi8(0x03);
    const __m128i carry = _mm_and_si128(value, msb_mask);
    const __m128i shifted = _mm_and_si128(_mm_slli_epi16(value, 1), low_mask);
    const __m128i carry_mask = _mm_cmpeq_epi8(carry, msb_mask);
    return _mm_xor_si128(shifted, _mm_and_si128(carry_mask, reduce));
}

static inline __m256i mayo_gf16_mul_precomputed_b(__m256i a,
                                                   __m256i b1,
                                                   __m256i b2,
                                                   __m256i b4,
                                                   __m256i b8) {
    const __m256i low_mask = _mm256_set1_epi8(0x0f);
    const __m256i one = _mm256_set1_epi8(1);
    const __m256i two = _mm256_set1_epi8(2);
    const __m256i four = _mm256_set1_epi8(4);
    const __m256i eight = _mm256_set1_epi8(8);

    a = _mm256_and_si256(a, low_mask);
    const __m256i m1 = _mm256_cmpeq_epi8(_mm256_and_si256(a, one), one);
    const __m256i m2 = _mm256_cmpeq_epi8(_mm256_and_si256(a, two), two);
    const __m256i m4 = _mm256_cmpeq_epi8(_mm256_and_si256(a, four), four);
    const __m256i m8 = _mm256_cmpeq_epi8(_mm256_and_si256(a, eight), eight);

    const __m256i x12 = _mm256_xor_si256(_mm256_and_si256(b1, m1),
                                         _mm256_and_si256(b2, m2));
    const __m256i x48 = _mm256_xor_si256(_mm256_and_si256(b4, m4),
                                         _mm256_and_si256(b8, m8));
    return _mm256_xor_si256(x12, x48);
}

static inline uint8_t mayo_gf16_hxor_256(__m256i value) {
    __m128i folded = _mm_xor_si128(_mm256_castsi256_si128(value),
                                   _mm256_extracti128_si256(value, 1));
    folded = _mm_xor_si128(folded, _mm_srli_si128(folded, 8));
    folded = _mm_xor_si128(folded, _mm_srli_si128(folded, 4));
    folded = _mm_xor_si128(folded, _mm_srli_si128(folded, 2));
    folded = _mm_xor_si128(folded, _mm_srli_si128(folded, 1));
    return (uint8_t)(_mm_cvtsi128_si32(folded) & 0x0f);
}

static inline __m128i mayo_gf16_mul_precomputed_b128(__m128i a,
                                                      __m128i b1,
                                                      __m128i b2,
                                                      __m128i b4,
                                                      __m128i b8) {
    const __m128i low_mask = _mm_set1_epi8(0x0f);
    const __m128i one = _mm_set1_epi8(1);
    const __m128i two = _mm_set1_epi8(2);
    const __m128i four = _mm_set1_epi8(4);
    const __m128i eight = _mm_set1_epi8(8);

    a = _mm_and_si128(a, low_mask);
    const __m128i m1 = _mm_cmpeq_epi8(_mm_and_si128(a, one), one);
    const __m128i m2 = _mm_cmpeq_epi8(_mm_and_si128(a, two), two);
    const __m128i m4 = _mm_cmpeq_epi8(_mm_and_si128(a, four), four);
    const __m128i m8 = _mm_cmpeq_epi8(_mm_and_si128(a, eight), eight);
    const __m128i x12 = _mm_xor_si128(_mm_and_si128(b1, m1),
                                      _mm_and_si128(b2, m2));
    const __m128i x48 = _mm_xor_si128(_mm_and_si128(b4, m4),
                                      _mm_and_si128(b8, m8));
    return _mm_xor_si128(x12, x48);
}

static inline uint8_t mayo_gf16_hxor_128(__m128i value) {
    value = _mm_xor_si128(value, _mm_srli_si128(value, 8));
    value = _mm_xor_si128(value, _mm_srli_si128(value, 4));
    value = _mm_xor_si128(value, _mm_srli_si128(value, 2));
    value = _mm_xor_si128(value, _mm_srli_si128(value, 1));
    return (uint8_t)(_mm_cvtsi128_si32(value) & 0x0f);
}

#ifdef MAYO_GF16_GFNI
static inline __m256i mayo_gf16_embed_gfni(__m256i value) {
    const __m256i embed = _mm256_setr_epi8(
        0x00, 0x01, 0x5c, 0x5d, (char)0xe0, (char)0xe1, (char)0xbc,
        (char)0xbd, 0x50, 0x51, 0x0c, 0x0d, (char)0xb0, (char)0xb1,
        (char)0xec, (char)0xed,
        0x00, 0x01, 0x5c, 0x5d, (char)0xe0, (char)0xe1, (char)0xbc,
        (char)0xbd, 0x50, 0x51, 0x0c, 0x0d, (char)0xb0, (char)0xb1,
        (char)0xec, (char)0xed);
    return _mm256_shuffle_epi8(embed, value);
}

static inline __m128i mayo_gf16_embed_gfni128(__m128i value) {
    const __m128i embed = _mm_setr_epi8(
        0x00, 0x01, 0x5c, 0x5d, (char)0xe0, (char)0xe1, (char)0xbc,
        (char)0xbd, 0x50, 0x51, 0x0c, 0x0d, (char)0xb0, (char)0xb1,
        (char)0xec, (char)0xed);
    return _mm_shuffle_epi8(embed, value);
}

static inline __m256i mayo_gf16_decode_gfni(__m256i value) {
    const __m256i inverse = _mm256_set1_epi64x(0x0104201400000000ULL);
    return _mm256_gf2p8affine_epi64_epi8(value, inverse, 0);
}

static inline __m128i mayo_gf16_decode_gfni128(__m128i value) {
    const __m128i inverse = _mm_set1_epi64x(0x0104201400000000ULL);
    return _mm_gf2p8affine_epi64_epi8(value, inverse, 0);
}

#ifdef MAYO_GF16_GFNI_CROSS_ROW
/*
 * GF_cal-style cross-row ILP helpers.  The right operand is embedded once;
 * four independent row accumulators are then advanced in the same chunk loop.
 * No table lookup or secret-dependent memory access is introduced.
 */
static inline uint8_t mayo_avx2_gfni_dot_skip_last(
    const unsigned char *row_ptr, const unsigned char *r,
    const __m256i *embedded_r, int chunks, __m128i embedded_tail,
    int vector_tail, int tail_offset, int terms) {
    __m256i product = _mm256_setzero_si256();
    for (int chunk = 0; chunk < chunks; ++chunk) {
        const __m256i a = mayo_gf16_embed_gfni(_mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(row_ptr + 32 * chunk)),
            _mm256_set1_epi8(0x0f)));
        product = _mm256_xor_si256(
            product, _mm256_gf2p8mul_epi8(a, embedded_r[chunk]));
    }
    if (vector_tail) {
        __m128i a;
        if (vector_tail == 16) {
            a = _mm_loadu_si128((const __m128i *)(row_ptr + tail_offset));
        } else if (vector_tail == 12) {
            a = _mm_maskload_epi32((const int *)(row_ptr + tail_offset),
                                   _mm_set_epi32(0, -1, -1, -1));
        } else if (vector_tail == 8) {
            a = _mm_maskload_epi32((const int *)(row_ptr + tail_offset),
                                   _mm_set_epi32(0, 0, -1, -1));
        } else {
            a = _mm_maskload_epi32((const int *)(row_ptr + tail_offset),
                                   _mm_set_epi32(0, 0, 0, -1));
        }
        const __m128i tail_product = _mm_gf2p8mul_epi8(
            mayo_gf16_embed_gfni128(_mm_and_si128(a, _mm_set1_epi8(0x0f))),
            embedded_tail);
        product = _mm256_xor_si256(
            product, _mm256_inserti128_si256(_mm256_setzero_si256(),
                                             tail_product, 0));
    }
    uint8_t result = mayo_gf16_hxor_256(mayo_gf16_decode_gfni(product));
    for (int i = tail_offset + vector_tail; i < terms; ++i) {
        result ^= gf16_mul_scalar_local(row_ptr[i], r[i]);
    }
    return result;
}

static inline uint8_t mayo_avx2_gfni_dot_short(
    const unsigned char *row_ptr, const unsigned char *b,
    __m128i embedded_b, __m128i load_mask, int vector_terms, int terms) {
    uint8_t result = 0;
    if (vector_terms) {
        const __m128i a = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)row_ptr, load_mask),
            _mm_set1_epi8(0x0f)));
        result = mayo_gf16_hxor_128(mayo_gf16_decode_gfni128(
            _mm_gf2p8mul_epi8(a, embedded_b)));
    }
    for (int i = vector_terms; i < terms; ++i) {
        result ^= gf16_mul_scalar_local(row_ptr[i], b[i]);
    }
    return result;
}

static inline void mayo_avx2_gfni_matvec_skip_last_cross_row(
    const unsigned char *A, const unsigned char *r, unsigned char *out,
    int rows, int stride) {
    const int terms = stride - 1;
    const int chunks = terms / 32;
    const int tail_offset = chunks * 32;
    const int tail = terms - tail_offset;
    __m256i embedded_r[MAYO_AVX2_MATVEC_MAX_CHUNKS];
    __m128i embedded_tail = _mm_setzero_si128();
    int vector_tail = 0;

    for (int chunk = 0; chunk < chunks; ++chunk) {
        const __m256i rv = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(r + 32 * chunk)),
            _mm256_set1_epi8(0x0f));
        embedded_r[chunk] = mayo_gf16_embed_gfni(rv);
    }
    if (tail >= 16) {
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_loadu_si128((const __m128i *)(r + tail_offset)),
            _mm_set1_epi8(0x0f)));
        vector_tail = 16;
    } else if (tail >= 12) {
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset),
                               _mm_set_epi32(0, -1, -1, -1)),
            _mm_set1_epi8(0x0f)));
        vector_tail = 12;
    } else if (tail >= 8) {
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset),
                               _mm_set_epi32(0, 0, -1, -1)),
            _mm_set1_epi8(0x0f)));
        vector_tail = 8;
    } else if (tail >= 4) {
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset),
                               _mm_set_epi32(0, 0, 0, -1)),
            _mm_set1_epi8(0x0f)));
        vector_tail = 4;
    }

    int row = 0;
    for (; row + 3 < rows; row += 4) {
        const unsigned char *row0 = A + (row + 0) * stride;
        const unsigned char *row1 = A + (row + 1) * stride;
        const unsigned char *row2 = A + (row + 2) * stride;
        const unsigned char *row3 = A + (row + 3) * stride;
        __m256i p0 = _mm256_setzero_si256();
        __m256i p1 = _mm256_setzero_si256();
        __m256i p2 = _mm256_setzero_si256();
        __m256i p3 = _mm256_setzero_si256();
        for (int chunk = 0; chunk < chunks; ++chunk) {
            const __m256i rv = embedded_r[chunk];
            const __m256i a0 = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_loadu_si256((const __m256i *)(row0 + 32 * chunk)),
                _mm256_set1_epi8(0x0f)));
            const __m256i a1 = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_loadu_si256((const __m256i *)(row1 + 32 * chunk)),
                _mm256_set1_epi8(0x0f)));
            const __m256i a2 = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_loadu_si256((const __m256i *)(row2 + 32 * chunk)),
                _mm256_set1_epi8(0x0f)));
            const __m256i a3 = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_loadu_si256((const __m256i *)(row3 + 32 * chunk)),
                _mm256_set1_epi8(0x0f)));
            p0 = _mm256_xor_si256(p0, _mm256_gf2p8mul_epi8(a0, rv));
            p1 = _mm256_xor_si256(p1, _mm256_gf2p8mul_epi8(a1, rv));
            p2 = _mm256_xor_si256(p2, _mm256_gf2p8mul_epi8(a2, rv));
            p3 = _mm256_xor_si256(p3, _mm256_gf2p8mul_epi8(a3, rv));
        }
        if (vector_tail) {
            __m128i a0, a1, a2, a3;
            if (vector_tail == 16) {
                a0 = _mm_loadu_si128((const __m128i *)(row0 + tail_offset));
                a1 = _mm_loadu_si128((const __m128i *)(row1 + tail_offset));
                a2 = _mm_loadu_si128((const __m128i *)(row2 + tail_offset));
                a3 = _mm_loadu_si128((const __m128i *)(row3 + tail_offset));
            } else {
                const __m128i mask = vector_tail == 12
                    ? _mm_set_epi32(0, -1, -1, -1)
                    : (vector_tail == 8 ? _mm_set_epi32(0, 0, -1, -1)
                                        : _mm_set_epi32(0, 0, 0, -1));
                a0 = _mm_maskload_epi32((const int *)(row0 + tail_offset), mask);
                a1 = _mm_maskload_epi32((const int *)(row1 + tail_offset), mask);
                a2 = _mm_maskload_epi32((const int *)(row2 + tail_offset), mask);
                a3 = _mm_maskload_epi32((const int *)(row3 + tail_offset), mask);
            }
            const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
            const __m256i embedded_tail_pair =
                _mm256_set_m128i(embedded_tail, embedded_tail);
            const __m256i e01 = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_set_m128i(a1, a0), nibble_mask));
            const __m256i e23 = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_set_m128i(a3, a2), nibble_mask));
            const __m256i tail01 =
                _mm256_gf2p8mul_epi8(e01, embedded_tail_pair);
            const __m256i tail23 =
                _mm256_gf2p8mul_epi8(e23, embedded_tail_pair);
            p0 = _mm256_xor_si256(p0, _mm256_inserti128_si256(
                _mm256_setzero_si256(), _mm256_castsi256_si128(tail01), 0));
            p1 = _mm256_xor_si256(p1, _mm256_inserti128_si256(
                _mm256_setzero_si256(), _mm256_extracti128_si256(tail01, 1), 0));
            p2 = _mm256_xor_si256(p2, _mm256_inserti128_si256(
                _mm256_setzero_si256(), _mm256_castsi256_si128(tail23), 0));
            p3 = _mm256_xor_si256(p3, _mm256_inserti128_si256(
                _mm256_setzero_si256(), _mm256_extracti128_si256(tail23, 1), 0));
        }
        out[row + 0] = mayo_gf16_hxor_256(mayo_gf16_decode_gfni(p0));
        out[row + 1] = mayo_gf16_hxor_256(mayo_gf16_decode_gfni(p1));
        out[row + 2] = mayo_gf16_hxor_256(mayo_gf16_decode_gfni(p2));
        out[row + 3] = mayo_gf16_hxor_256(mayo_gf16_decode_gfni(p3));
        for (int i = tail_offset + vector_tail; i < terms; ++i) {
            out[row + 0] ^= gf16_mul_scalar_local(row0[i], r[i]);
            out[row + 1] ^= gf16_mul_scalar_local(row1[i], r[i]);
            out[row + 2] ^= gf16_mul_scalar_local(row2[i], r[i]);
            out[row + 3] ^= gf16_mul_scalar_local(row3[i], r[i]);
        }
    }
    for (; row < rows; ++row) {
        out[row] = mayo_avx2_gfni_dot_skip_last(
            A + row * stride, r, embedded_r, chunks, embedded_tail,
            vector_tail, tail_offset, terms);
    }
}

static inline void mayo_avx2_gfni_matvec_short_cross_row(
    const unsigned char *A, const unsigned char *b, unsigned char *out,
    int rows, int terms) {
    __m128i load_mask;
    int vector_terms;
    if (terms >= 16) {
        load_mask = _mm_set1_epi32(-1);
        vector_terms = 16;
    } else if (terms >= 12) {
        load_mask = _mm_set_epi32(0, -1, -1, -1);
        vector_terms = 12;
    } else if (terms >= 8) {
        load_mask = _mm_set_epi32(0, 0, -1, -1);
        vector_terms = 8;
    } else if (terms >= 4) {
        load_mask = _mm_set_epi32(0, 0, 0, -1);
        vector_terms = 4;
    } else {
        load_mask = _mm_setzero_si128();
        vector_terms = 0;
    }
    __m128i embedded_b = _mm_setzero_si128();
    if (vector_terms) {
        embedded_b = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)b, load_mask),
            _mm_set1_epi8(0x0f)));
    }
    const __m256i embedded_b_pair = _mm256_set_m128i(embedded_b, embedded_b);
    const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
    int row = 0;
    for (; row + 1 < rows; row += 2) {
        const unsigned char *row0 = A + row * terms;
        const unsigned char *row1 = row0 + terms;
        uint8_t result0 = 0;
        uint8_t result1 = 0;
        if (vector_terms) {
            const __m128i a0 = _mm_maskload_epi32((const int *)row0,
                                                  load_mask);
            const __m128i a1 = _mm_maskload_epi32((const int *)row1,
                                                  load_mask);
            const __m256i pair = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_set_m128i(a1, a0), nibble_mask));
            const __m256i product = mayo_gf16_decode_gfni(
                _mm256_gf2p8mul_epi8(pair, embedded_b_pair));
            result0 = mayo_gf16_hxor_128(_mm256_castsi256_si128(product));
            result1 = mayo_gf16_hxor_128(
                _mm256_extracti128_si256(product, 1));
        }
        for (int i = vector_terms; i < terms; ++i) {
            result0 ^= gf16_mul_scalar_local(row0[i], b[i]);
            result1 ^= gf16_mul_scalar_local(row1[i], b[i]);
        }
        out[row] = result0;
        out[row + 1] = result1;
    }
    for (; row < rows; ++row) {
        out[row] = mayo_avx2_gfni_dot_short(
            A + row * terms, b, embedded_b, load_mask, vector_terms, terms);
    }
}
#endif

static inline void mayo_avx2_gfni_matvec_skip_last(const unsigned char *A,
                                                   const unsigned char *r,
                                                   unsigned char *out,
                                                   int rows,
                                                   int stride) {
    const int terms = stride - 1;
    const int chunks = terms / 32;
    const int tail_offset = chunks * 32;
    const int tail = terms - tail_offset;
    __m256i embedded_r[MAYO_AVX2_MATVEC_MAX_CHUNKS];
    __m128i embedded_tail = _mm_setzero_si128();
    int vector_tail = 0;

    for (int chunk = 0; chunk < chunks; ++chunk) {
        const __m256i rv = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(r + 32 * chunk)),
            _mm256_set1_epi8(0x0f));
        embedded_r[chunk] = mayo_gf16_embed_gfni(rv);
    }

    if (tail >= 16) {
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_loadu_si128((const __m128i *)(r + tail_offset)),
            _mm_set1_epi8(0x0f)));
        vector_tail = 16;
    } else if (tail >= 12) {
        const __m128i mask3 = _mm_set_epi32(0, -1, -1, -1);
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset), mask3),
            _mm_set1_epi8(0x0f)));
        vector_tail = 12;
    } else if (tail >= 8) {
        const __m128i mask2 = _mm_set_epi32(0, 0, -1, -1);
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset), mask2),
            _mm_set1_epi8(0x0f)));
        vector_tail = 8;
    } else if (tail >= 4) {
        const __m128i mask1 = _mm_set_epi32(0, 0, 0, -1);
        embedded_tail = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset), mask1),
            _mm_set1_epi8(0x0f)));
        vector_tail = 4;
    }

    for (int row = 0; row < rows; ++row) {
        const unsigned char *row_ptr = A + row * stride;
        __m256i product = _mm256_setzero_si256();
        for (int chunk = 0; chunk < chunks; ++chunk) {
            const __m256i a = mayo_gf16_embed_gfni(_mm256_and_si256(
                _mm256_loadu_si256((const __m256i *)(row_ptr + 32 * chunk)),
                _mm256_set1_epi8(0x0f)));
            product = _mm256_xor_si256(
                product, _mm256_gf2p8mul_epi8(a, embedded_r[chunk]));
        }
        if (vector_tail) {
            __m128i a;
            const int *tail_ptr = (const int *)(row_ptr + tail_offset);
            if (vector_tail == 16) {
                a = _mm_loadu_si128((const __m128i *)tail_ptr);
            } else if (vector_tail == 12) {
                a = _mm_maskload_epi32(tail_ptr, _mm_set_epi32(0, -1, -1, -1));
            } else if (vector_tail == 8) {
                a = _mm_maskload_epi32(tail_ptr, _mm_set_epi32(0, 0, -1, -1));
            } else {
                a = _mm_maskload_epi32(tail_ptr, _mm_set_epi32(0, 0, 0, -1));
            }
            const __m128i tail_product = _mm_gf2p8mul_epi8(
                mayo_gf16_embed_gfni128(_mm_and_si128(a, _mm_set1_epi8(0x0f))),
                embedded_tail);
            product = _mm256_xor_si256(
                product, _mm256_inserti128_si256(_mm256_setzero_si256(),
                                                 tail_product, 0));
        }
        uint8_t result = mayo_gf16_hxor_256(mayo_gf16_decode_gfni(product));
        for (int i = tail_offset + vector_tail; i < terms; ++i) {
            result ^= gf16_mul_scalar_local(row_ptr[i], r[i]);
        }
        out[row] = result;
    }
}

static inline void mayo_avx2_gfni_matvec_short(const unsigned char *A,
                                               const unsigned char *b,
                                               unsigned char *out,
                                               int rows,
                                               int terms) {
    __m128i load_mask;
    int vector_terms;
    if (terms >= 16) {
        load_mask = _mm_set1_epi32(-1);
        vector_terms = 16;
    } else if (terms >= 12) {
        load_mask = _mm_set_epi32(0, -1, -1, -1);
        vector_terms = 12;
    } else if (terms >= 8) {
        load_mask = _mm_set_epi32(0, 0, -1, -1);
        vector_terms = 8;
    } else if (terms >= 4) {
        load_mask = _mm_set_epi32(0, 0, 0, -1);
        vector_terms = 4;
    } else {
        load_mask = _mm_setzero_si128();
        vector_terms = 0;
    }

    __m128i embedded_b = _mm_setzero_si128();
    if (vector_terms) {
        embedded_b = mayo_gf16_embed_gfni128(_mm_and_si128(
            _mm_maskload_epi32((const int *)b, load_mask),
            _mm_set1_epi8(0x0f)));
    }
    for (int row = 0; row < rows; ++row) {
        const unsigned char *row_ptr = A + row * terms;
        uint8_t result = 0;
        if (vector_terms) {
            const __m128i a = mayo_gf16_embed_gfni128(_mm_and_si128(
                _mm_maskload_epi32((const int *)row_ptr, load_mask),
                _mm_set1_epi8(0x0f)));
            result = mayo_gf16_hxor_128(mayo_gf16_decode_gfni128(
                _mm_gf2p8mul_epi8(a, embedded_b)));
        }
        for (int i = vector_terms; i < terms; ++i) {
            result ^= gf16_mul_scalar_local(row_ptr[i], b[i]);
        }
        out[row] = result;
    }
}
#endif

/*
 * Compute A * r for sample_solution.  The final matrix column is the
 * right-hand-side slot and is known to be zero at this call site, so only the
 * first (stride - 1) coefficients are consumed.
 */
static inline void mayo_avx2_matvec_skip_last(const unsigned char *A,
                                              const unsigned char *r,
                                              unsigned char *out,
                                              int rows,
                                              int stride) {
    const int terms = stride - 1;
    const int chunks = terms / 32;
    const int tail_offset = chunks * 32;
    const int tail = terms - tail_offset;
    __m256i b1[MAYO_AVX2_MATVEC_MAX_CHUNKS];
    __m256i b2[MAYO_AVX2_MATVEC_MAX_CHUNKS];
    __m256i b4[MAYO_AVX2_MATVEC_MAX_CHUNKS];
    __m256i b8[MAYO_AVX2_MATVEC_MAX_CHUNKS];
    __m128i tail_b1 = _mm_setzero_si128();
    __m128i tail_b2 = _mm_setzero_si128();
    __m128i tail_b4 = _mm_setzero_si128();
    __m128i tail_b8 = _mm_setzero_si128();
    int vector_tail = 0;

    for (int chunk = 0; chunk < chunks; ++chunk) {
        b1[chunk] = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(r + 32 * chunk)),
            _mm256_set1_epi8(0x0f));
        b2[chunk] = mayo_gf16_mul_x_vec(b1[chunk]);
        b4[chunk] = mayo_gf16_mul_x_vec(b2[chunk]);
        b8[chunk] = mayo_gf16_mul_x_vec(b4[chunk]);
    }

    if (tail >= 16) {
        tail_b1 = _mm_and_si128(
            _mm_loadu_si128((const __m128i *)(r + tail_offset)),
            _mm_set1_epi8(0x0f));
        vector_tail = 16;
    } else if (tail >= 12) {
        const __m128i mask3 = _mm_set_epi32(0, -1, -1, -1);
        tail_b1 = _mm_and_si128(
            _mm_maskload_epi32((const int *)(r + tail_offset), mask3),
            _mm_set1_epi8(0x0f));
        vector_tail = 12;
    }
    if (vector_tail) {
        tail_b2 = mayo_gf16_mul_x_vec128(tail_b1);
        tail_b4 = mayo_gf16_mul_x_vec128(tail_b2);
        tail_b8 = mayo_gf16_mul_x_vec128(tail_b4);
    }

    for (int row = 0; row < rows; ++row) {
        const unsigned char *row_ptr = A + row * stride;
        uint8_t result = 0;
        for (int chunk = 0; chunk < chunks; ++chunk) {
            const __m256i a = _mm256_loadu_si256(
                (const __m256i *)(row_ptr + 32 * chunk));
            const __m256i product = mayo_gf16_mul_precomputed_b(
                a, b1[chunk], b2[chunk], b4[chunk], b8[chunk]);
            result ^= mayo_gf16_hxor_256(product);
        }
        if (vector_tail) {
            __m128i a;
            if (vector_tail == 16) {
                a = _mm_loadu_si128((const __m128i *)(row_ptr + tail_offset));
            } else {
                const __m128i mask3 = _mm_set_epi32(0, -1, -1, -1);
                a = _mm_maskload_epi32((const int *)(row_ptr + tail_offset), mask3);
            }
            result ^= mayo_gf16_hxor_128(mayo_gf16_mul_precomputed_b128(
                a, tail_b1, tail_b2, tail_b4, tail_b8));
        }
        for (int i = tail_offset + vector_tail; i < terms; ++i) {
            result ^= gf16_mul_scalar_local(row_ptr[i], r[i]);
        }
        out[row] = result;
    }
}

/* Short dense matrix-vector products used for O*x after solving. */
static inline void mayo_avx2_matvec_short(const unsigned char *A,
                                          const unsigned char *b,
                                          unsigned char *out,
                                          int rows,
                                          int terms) {
    __m128i load_mask;
    int vector_terms;
    if (terms >= 16) {
        load_mask = _mm_set1_epi32(-1);
        vector_terms = 16;
    } else if (terms >= 12) {
        load_mask = _mm_set_epi32(0, -1, -1, -1);
        vector_terms = 12;
    } else if (terms >= 8) {
        load_mask = _mm_set_epi32(0, 0, -1, -1);
        vector_terms = 8;
    } else if (terms >= 4) {
        load_mask = _mm_set_epi32(0, 0, 0, -1);
        vector_terms = 4;
    } else {
        load_mask = _mm_setzero_si128();
        vector_terms = 0;
    }

    __m128i b1 = _mm_setzero_si128();
    __m128i b2 = _mm_setzero_si128();
    __m128i b4 = _mm_setzero_si128();
    __m128i b8 = _mm_setzero_si128();
    if (vector_terms) {
        b1 = _mm_and_si128(_mm_maskload_epi32((const int *)b, load_mask),
                           _mm_set1_epi8(0x0f));
        b2 = mayo_gf16_mul_x_vec128(b1);
        b4 = mayo_gf16_mul_x_vec128(b2);
        b8 = mayo_gf16_mul_x_vec128(b4);
    }

    for (int row = 0; row < rows; ++row) {
        const unsigned char *row_ptr = A + row * terms;
        uint8_t result = 0;
        if (vector_terms) {
            const __m128i a = _mm_maskload_epi32((const int *)row_ptr, load_mask);
            result = mayo_gf16_hxor_128(mayo_gf16_mul_precomputed_b128(
                a, b1, b2, b4, b8));
        }
        for (int i = vector_terms; i < terms; ++i) {
            result ^= gf16_mul_scalar_local(row_ptr[i], b[i]);
        }
        out[row] = result;
    }
}

#endif
