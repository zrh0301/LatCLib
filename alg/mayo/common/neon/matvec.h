/* SPDX-License-Identifier: Apache-2.0 */

#ifndef MAYO_NEON_MATVEC_H
#define MAYO_NEON_MATVEC_H

#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <arithmetic_common.h>
#include <simple_arithmetic.h>

#if defined(__GNUC__) || defined(__clang__)
#define MAYO_NEON_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
#define MAYO_NEON_ALWAYS_INLINE static inline
#endif

/*
 * ARM candidate path for the two matrix-vector products that are otherwise
 * handled by scalar mat_mul().  The input and output format is the same as
 * mat_mul(): one GF(16) value per byte, low nibble significant, row-major A.
 */

static inline uint8x16_t mayo_neon_mul_x(uint8x16_t x) {
    const uint8x16_t low_mask = vdupq_n_u8(0x0f);
    const uint8x16_t high_bit = vdupq_n_u8(0x08);
    const uint8x16_t reduction = vdupq_n_u8(0x03);
    x = vandq_u8(x, low_mask);
    const uint8x16_t carry = vtstq_u8(x, high_bit);
    return vandq_u8(
        veorq_u8(vshlq_n_u8(x, 1), vandq_u8(carry, reduction)),
        low_mask);
}

static inline void mayo_neon_make_coeff_masks(uint8x16_t a,
                                              uint8x16_t *m1,
                                              uint8x16_t *m2,
                                              uint8x16_t *m4,
                                              uint8x16_t *m8) {
    const uint8x16_t low_mask = vdupq_n_u8(0x0f);
    const uint8x16_t one = vdupq_n_u8(0x01);
    const uint8x16_t two = vdupq_n_u8(0x02);
    const uint8x16_t four = vdupq_n_u8(0x04);
    const uint8x16_t eight = vdupq_n_u8(0x08);

    a = vandq_u8(a, low_mask);
    *m1 = vtstq_u8(a, one);
    *m2 = vtstq_u8(a, two);
    *m4 = vtstq_u8(a, four);
    *m8 = vtstq_u8(a, eight);
}

static inline uint8x16_t mayo_neon_mul_vec_precomputed_masks(
    uint8x16_t b, uint8x16_t b2, uint8x16_t b4, uint8x16_t b8,
    uint8x16_t m1, uint8x16_t m2, uint8x16_t m4, uint8x16_t m8) {
    const uint8x16_t low_mask = vdupq_n_u8(0x0f);
    b = vandq_u8(b, low_mask);

    const uint8x16_t t1 = vandq_u8(b, m1);
    const uint8x16_t t2 = vandq_u8(b2, m2);
    const uint8x16_t t4 = vandq_u8(b4, m4);
    const uint8x16_t t8 = vandq_u8(b8, m8);
    return veorq_u8(veorq_u8(t1, t2), veorq_u8(t4, t8));
}

static inline uint8x16_t mayo_neon_mul_vec_precomputed(uint8x16_t a,
                                                       uint8x16_t b,
                                                       uint8x16_t b2,
                                                       uint8x16_t b4,
                                                       uint8x16_t b8) {
    uint8x16_t m1, m2, m4, m8;
    mayo_neon_make_coeff_masks(a, &m1, &m2, &m4, &m8);
    return mayo_neon_mul_vec_precomputed_masks(b, b2, b4, b8,
                                               m1, m2, m4, m8);
}

static inline uint8x16_t mayo_neon_mul_vec(uint8x16_t a, uint8x16_t b) {
    const uint8x16_t b2 = mayo_neon_mul_x(b);
    const uint8x16_t b4 = mayo_neon_mul_x(b2);
    const uint8x16_t b8 = mayo_neon_mul_x(b4);
    return mayo_neon_mul_vec_precomputed(a, b, b2, b4, b8);
}

static inline unsigned char mayo_neon_hxor16(uint8x16_t x) {
    uint8x8_t h = veor_u8(vget_low_u8(x), vget_high_u8(x));
    h = veor_u8(h, vext_u8(h, h, 4));
    h = veor_u8(h, vext_u8(h, h, 2));
    h = veor_u8(h, vext_u8(h, h, 1));
    return (unsigned char)(vget_lane_u8(h, 0) & 0x0f);
}

/* Load a short MAYO row without reading into the next row. */
static inline uint8x16_t mayo_neon_load_short(const unsigned char *p,
                                              int cols) {
    uint8x8_t low = vld1_u8(p);
    uint8x8_t high = vdup_n_u8(0x00);

    if (cols >= 9) {
        high = vld1_lane_u8(p + 8, high, 0);
    }
    if (cols >= 10) {
        high = vld1_lane_u8(p + 9, high, 1);
    }
    if (cols >= 11) {
        high = vld1_lane_u8(p + 10, high, 2);
    }
    if (cols >= 12) {
        high = vld1_lane_u8(p + 11, high, 3);
    }
    return vcombine_u8(low, high);
}

static inline void mayo_neon_matvec_strided(const unsigned char *a,
                                            const unsigned char *b,
                                            unsigned char *c,
                                            int rows, int cols,
                                            int row_stride) {
    const int full_chunks = cols / 16;
    uint8x16_t b2_chunks[full_chunks > 0 ? full_chunks : 1];
    uint8x16_t b4_chunks[full_chunks > 0 ? full_chunks : 1];
    uint8x16_t b8_chunks[full_chunks > 0 ? full_chunks : 1];

    for (int block = 0; block < full_chunks; ++block) {
        const uint8x16_t bv = vld1q_u8(b + block * 16);
        b2_chunks[block] = mayo_neon_mul_x(bv);
        b4_chunks[block] = mayo_neon_mul_x(b2_chunks[block]);
        b8_chunks[block] = mayo_neon_mul_x(b4_chunks[block]);
    }

    /*
     * Keep four independent accumulators live while walking the same b
     * stream.  The previous one-row loop was mathematically identical, but
     * left NEON dependent-XOR latency exposed between consecutive rows.
     * Four rows also matches the natural width of the MAYO row batches while
     * keeping the temporary state in registers.
     */
    int row = 0;
    for (; row + 3 < rows; row += 4) {
        const unsigned char *a0 = a + (size_t)(row + 0) * (size_t)row_stride;
        const unsigned char *a1 = a + (size_t)(row + 1) * (size_t)row_stride;
        const unsigned char *a2 = a + (size_t)(row + 2) * (size_t)row_stride;
        const unsigned char *a3 = a + (size_t)(row + 3) * (size_t)row_stride;
        uint8x16_t acc0 = vdupq_n_u8(0x00);
        uint8x16_t acc1 = vdupq_n_u8(0x00);
        uint8x16_t acc2 = vdupq_n_u8(0x00);
        uint8x16_t acc3 = vdupq_n_u8(0x00);

        for (int block = 0, col = 0; block < full_chunks; ++block, col += 16) {
            const uint8x16_t bv = vld1q_u8(b + col);
            const uint8x16_t b2 = b2_chunks[block];
            const uint8x16_t b4 = b4_chunks[block];
            const uint8x16_t b8 = b8_chunks[block];
            acc0 = veorq_u8(acc0, mayo_neon_mul_vec_precomputed(vld1q_u8(a0 + col), bv, b2, b4, b8));
            acc1 = veorq_u8(acc1, mayo_neon_mul_vec_precomputed(vld1q_u8(a1 + col), bv, b2, b4, b8));
            acc2 = veorq_u8(acc2, mayo_neon_mul_vec_precomputed(vld1q_u8(a2 + col), bv, b2, b4, b8));
            acc3 = veorq_u8(acc3, mayo_neon_mul_vec_precomputed(vld1q_u8(a3 + col), bv, b2, b4, b8));
        }

        unsigned char t0 = mayo_neon_hxor16(acc0);
        unsigned char t1 = mayo_neon_hxor16(acc1);
        unsigned char t2 = mayo_neon_hxor16(acc2);
        unsigned char t3 = mayo_neon_hxor16(acc3);
        for (int col = full_chunks * 16; col < cols; ++col) {
            const unsigned char bv = b[col];
            t0 ^= mul_f(a0[col], bv);
            t1 ^= mul_f(a1[col], bv);
            t2 ^= mul_f(a2[col], bv);
            t3 ^= mul_f(a3[col], bv);
        }
        c[row + 0] = t0;
        c[row + 1] = t1;
        c[row + 2] = t2;
        c[row + 3] = t3;
    }

    for (; row < rows; ++row) {
        const unsigned char *arow = a + (size_t)row * (size_t)row_stride;
        uint8x16_t acc = vdupq_n_u8(0x00);
        int col = 0;
        for (int block = 0; block < full_chunks; ++block, col += 16) {
            const uint8x16_t bv = vld1q_u8(b + col);
            acc = veorq_u8(acc, mayo_neon_mul_vec_precomputed(
                vld1q_u8(arow + col), bv, b2_chunks[block],
                b4_chunks[block], b8_chunks[block]));
        }
        unsigned char tail = mayo_neon_hxor16(acc);
        for (; col < cols; ++col) {
            tail ^= mul_f(arow[col], b[col]);
        }
        c[row] = tail;
    }
}

static inline void mayo_neon_matvec_skip_last(const unsigned char *a,
                                              const unsigned char *b,
                                              unsigned char *c,
                                              int rows, int cols_with_rhs) {
    mayo_neon_matvec_strided(a, b, c, rows, cols_with_rhs - 1,
                             cols_with_rhs);
}

/*
 * Short-width O*x path for MAYO-1/2/3/5.  These parameter sets use eight,
 * ten or twelve columns, so a row must be loaded with a bounded helper
 * instead of a 16-byte load.  Padding lanes stay zero and are harmless to
 * the horizontal XOR reduction.
 */
static inline void mayo_neon_matvec_short_small(const unsigned char *a,
                                                const unsigned char *b,
                                                unsigned char *c,
                                                int rows, int cols) {
    const uint8x16_t bv = mayo_neon_load_short(b, cols);
    const uint8x16_t b2 = mayo_neon_mul_x(bv);
    const uint8x16_t b4 = mayo_neon_mul_x(b2);
    const uint8x16_t b8 = mayo_neon_mul_x(b4);
    int row = 0;

    for (; row + 3 < rows; row += 4) {
        const unsigned char *a0 = a + (size_t)(row + 0) * (size_t)cols;
        const unsigned char *a1 = a + (size_t)(row + 1) * (size_t)cols;
        const unsigned char *a2 = a + (size_t)(row + 2) * (size_t)cols;
        const unsigned char *a3 = a + (size_t)(row + 3) * (size_t)cols;
        const uint8x16_t x0 = mayo_neon_load_short(a0, cols);
        const uint8x16_t x1 = mayo_neon_load_short(a1, cols);
        const uint8x16_t x2 = mayo_neon_load_short(a2, cols);
        const uint8x16_t x3 = mayo_neon_load_short(a3, cols);
        const uint8x16_t y0 = mayo_neon_mul_vec_precomputed(x0, bv, b2, b4, b8);
        const uint8x16_t y1 = mayo_neon_mul_vec_precomputed(x1, bv, b2, b4, b8);
        const uint8x16_t y2 = mayo_neon_mul_vec_precomputed(x2, bv, b2, b4, b8);
        const uint8x16_t y3 = mayo_neon_mul_vec_precomputed(x3, bv, b2, b4, b8);

        c[row + 0] = mayo_neon_hxor16(y0);
        c[row + 1] = mayo_neon_hxor16(y1);
        c[row + 2] = mayo_neon_hxor16(y2);
        c[row + 3] = mayo_neon_hxor16(y3);
    }

    for (; row < rows; ++row) {
        const unsigned char *arow = a + (size_t)row * (size_t)cols;
        const uint8x16_t x = mayo_neon_load_short(arow, cols);
        c[row] = mayo_neon_hxor16(
            mayo_neon_mul_vec_precomputed(x, bv, b2, b4, b8));
    }
}

MAYO_NEON_ALWAYS_INLINE void mayo_neon_prepare_short_masks(const unsigned char *a,
                                                           int rows, int cols,
                                                           uint8x16_t *m1,
                                                           uint8x16_t *m2,
                                                           uint8x16_t *m4,
                                                           uint8x16_t *m8) {
    for (int row = 0; row < rows; ++row) {
        const unsigned char *arow = a + (size_t)row * (size_t)cols;
        const uint8x16_t av = mayo_neon_load_short(arow, cols);
        mayo_neon_make_coeff_masks(av, &m1[row], &m2[row],
                                   &m4[row], &m8[row]);
    }
}

MAYO_NEON_ALWAYS_INLINE void mayo_neon_matvec_short_prepared(
    const unsigned char *b, unsigned char *c, int rows, int cols,
    const uint8x16_t *m1, const uint8x16_t *m2,
    const uint8x16_t *m4, const uint8x16_t *m8) {
    const uint8x16_t bv = mayo_neon_load_short(b, cols);
    const uint8x16_t b2 = mayo_neon_mul_x(bv);
    const uint8x16_t b4 = mayo_neon_mul_x(b2);
    const uint8x16_t b8 = mayo_neon_mul_x(b4);
    int row = 0;

    for (; row + 3 < rows; row += 4) {
        c[row + 0] = mayo_neon_hxor16(
            mayo_neon_mul_vec_precomputed_masks(bv, b2, b4, b8,
                                                m1[row + 0], m2[row + 0],
                                                m4[row + 0], m8[row + 0]));
        c[row + 1] = mayo_neon_hxor16(
            mayo_neon_mul_vec_precomputed_masks(bv, b2, b4, b8,
                                                m1[row + 1], m2[row + 1],
                                                m4[row + 1], m8[row + 1]));
        c[row + 2] = mayo_neon_hxor16(
            mayo_neon_mul_vec_precomputed_masks(bv, b2, b4, b8,
                                                m1[row + 2], m2[row + 2],
                                                m4[row + 2], m8[row + 2]));
        c[row + 3] = mayo_neon_hxor16(
            mayo_neon_mul_vec_precomputed_masks(bv, b2, b4, b8,
                                                m1[row + 3], m2[row + 3],
                                                m4[row + 3], m8[row + 3]));
    }

    for (; row < rows; ++row) {
        c[row] = mayo_neon_hxor16(
            mayo_neon_mul_vec_precomputed_masks(bv, b2, b4, b8,
                                                m1[row], m2[row],
                                                m4[row], m8[row]));
    }
}

static inline void mayo_neon_matvec_short(const unsigned char *a,
                                          const unsigned char *b,
                                          unsigned char *c,
                                          int rows, int cols) {
    if (cols >= 8 && cols <= 12) {
        mayo_neon_matvec_short_small(a, b, c, rows, cols);
        return;
    }
    mayo_neon_matvec_strided(a, b, c, rows, cols, cols);
}

#if defined(MAYO_NEON_OPT)
static inline int mayo_neon_round_rows16(int rows) {
    return (rows + 15) & ~15;
}

MAYO_NEON_ALWAYS_INLINE void mayo_neon_prepare_short_columns(
    const unsigned char *a, unsigned char *columns, int rows, int cols) {
    const int padded_rows = mayo_neon_round_rows16(rows);
    memset(columns, 0, (size_t)padded_rows * (size_t)cols);
    for (int col = 0; col < cols; ++col) {
        unsigned char *dst = columns + (size_t)col * (size_t)padded_rows;
        for (int row = 0; row < rows; ++row) {
            dst[row] = a[(size_t)row * (size_t)cols + (size_t)col];
        }
    }
}

MAYO_NEON_ALWAYS_INLINE void mayo_neon_matvec_short_columns(
    const unsigned char *columns, const unsigned char *b, unsigned char *c,
    int rows, int cols) {
    const int padded_rows = mayo_neon_round_rows16(rows);
    const int full_rows = rows & ~15;
    int row = 0;

    for (; row < full_rows; row += 16) {
        uint8x16_t acc = vdupq_n_u8(0);
        for (int col = 0; col < cols; ++col) {
            const uint8x16_t table = gf16v_get_multab(b[col]);
            const uint8x16_t coeff = vld1q_u8(
                columns + (size_t)col * (size_t)padded_rows + (size_t)row);
            acc = veorq_u8(acc, vqtbl1q_u8(table, coeff));
        }
        vst1q_u8(c + row, acc);
    }

    if (row < rows) {
        uint8x16_t acc = vdupq_n_u8(0);
        unsigned char tail[16];
        for (int col = 0; col < cols; ++col) {
            const uint8x16_t table = gf16v_get_multab(b[col]);
            const uint8x16_t coeff = vld1q_u8(
                columns + (size_t)col * (size_t)padded_rows + (size_t)row);
            acc = veorq_u8(acc, vqtbl1q_u8(table, coeff));
        }
        vst1q_u8(tail, acc);
        memcpy(c + row, tail, (size_t)(rows - row));
    }
}
#endif

#endif

#undef MAYO_NEON_ALWAYS_INLINE
