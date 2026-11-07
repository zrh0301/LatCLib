#include "gfext.h"

uint8_t pqc_kernel_gf2_parity_u64(uint64_t x) {
    x ^= x >> 32;
    x ^= x >> 16;
    x ^= x >> 8;
    x ^= x >> 4;
    x ^= x >> 2;
    x ^= x >> 1;
    return (uint8_t)(x & 1u);
}

uint8_t pqc_kernel_gf2_dot_u64(uint64_t a, uint64_t b) {
    return pqc_kernel_gf2_parity_u64(a & b);
}

void pqc_kernel_gf2_xor_words(uint64_t *dst, const uint64_t *src, size_t words) {
    for (size_t i = 0; i < words; ++i) {
        dst[i] ^= src[i];
    }
}

uint8_t pqc_kernel_gf2_dot_words(const uint64_t *a, const uint64_t *b, size_t words) {
    uint8_t result = 0;
    for (size_t i = 0; i < words; ++i) {
        result ^= pqc_kernel_gf2_dot_u64(a[i], b[i]);
    }
    return result;
}

void pqc_kernel_gf2_mat_vec_u64(uint64_t *out, const uint64_t *matrix_rows,
                     const uint64_t *vector, size_t rows, size_t words_per_row) {
    const size_t out_words = (rows + 63u) / 64u;
    for (size_t i = 0; i < out_words; ++i) {
        out[i] = 0;
    }
    for (size_t row = 0; row < rows; ++row) {
        const uint8_t bit = pqc_kernel_gf2_dot_words(matrix_rows + row * words_per_row,
                                          vector, words_per_row);
        out[row >> 6] |= (uint64_t)bit << (row & 63u);
    }
}
