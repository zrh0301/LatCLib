#include "gf16.h"

#include <arm_neon.h>
#include <stdint.h>
#include <string.h>

/* See the x86 implementation for the packed layout. */

static inline uint64_t load_u64(const uint8_t *p) {
    uint64_t value;
    memcpy(&value, p, sizeof(value));
    return value;
}

static inline void store_u64(uint8_t *p, uint64_t value) {
    memcpy(p, &value, sizeof(value));
}

static inline uint64_t clmul_pair(uint64_t a, uint64_t b) {
    const poly128_t product = vmull_p64((poly64_t)a, (poly64_t)b);
    const uint64x2_t words = vreinterpretq_u64_p128(product);
    return vgetq_lane_u64(words, 0);
}

static inline uint32_t reduce7(uint32_t p) {
    const uint32_t high = (p >> 4) & 0x07u;
    return (p ^ high ^ (high << 1)) & 0x0fu;
}

static inline uint32_t reduce_pair_product(uint64_t product) {
    const uint32_t p0 = (uint32_t)(product & UINT64_C(0x7f));
    const uint32_t p1 = (uint32_t)((product >> 16) & UINT64_C(0x7f));
    return reduce7(p0) | (reduce7(p1) << 8);
}

size_t pqc_kernel_gf16_clmul_packed_bytes(size_t n) {
    return ((n + 1u) / 2u) * 8u;
}

void pqc_kernel_gf16_clmul_pack_array(uint8_t *packed, const uint8_t *raw, size_t n) {
    size_t i = 0;
    for (; i + 8u <= n; i += 8u) {
        const uint64_t values = load_u64(raw + i);
        for (size_t pair = 0; pair < 4u; ++pair) {
            const uint64_t value =
                (values >> (pair * 16u)) & UINT64_C(0x0f0f);
            store_u64(packed + (i / 2u + pair) * 8u, value);
        }
    }
    for (; i + 1u < n; i += 2u) {
        const uint64_t value = (uint64_t)(raw[i] & 0x0f) |
                               ((uint64_t)(raw[i + 1u] & 0x0f) << 8);
        store_u64(packed + (i / 2u) * 8u, value);
    }
    if (i < n) {
        store_u64(packed + (i / 2u) * 8u, (uint64_t)(raw[i] & 0x0f));
    }
}

void pqc_kernel_gf16_clmul_unpack_array(uint8_t *raw, const uint8_t *packed, size_t n) {
    size_t i = 0;
    size_t group = 0;
    for (; i + 8u <= n; i += 8u, group += 4u) {
        for (size_t pair = 0; pair < 4u; ++pair) {
            const uint64_t value = load_u64(packed + (group + pair) * 8u);
            raw[i + 2u * pair] = (uint8_t)(value & 0x0f);
            raw[i + 2u * pair + 1u] = (uint8_t)((value >> 8) & 0x0f);
        }
    }
    for (; i + 1u < n; i += 2u, ++group) {
        const uint64_t value = load_u64(packed + group * 8u);
        raw[i] = (uint8_t)(value & 0x0f);
        raw[i + 1u] = (uint8_t)((value >> 8) & 0x0f);
    }
    if (i < n) raw[i] = (uint8_t)(load_u64(packed + group * 8u) & 0x0f);
}

void pqc_kernel_gf16_clmul_product_array(uint8_t *product, const uint8_t *a_packed,
                              const uint8_t *b_packed, size_t n) {
    const size_t groups = (n + 1u) / 2u;
    size_t g = 0;
    for (; g + 4u <= groups; g += 4u) {
        const uint64_t a0 = load_u64(a_packed + (g + 0u) * 8u);
        const uint64_t a1 = load_u64(a_packed + (g + 1u) * 8u);
        const uint64_t a2 = load_u64(a_packed + (g + 2u) * 8u);
        const uint64_t a3 = load_u64(a_packed + (g + 3u) * 8u);
        const uint64_t b0 = load_u64(b_packed + (g + 0u) * 8u);
        const uint64_t b1 = load_u64(b_packed + (g + 1u) * 8u);
        const uint64_t b2 = load_u64(b_packed + (g + 2u) * 8u);
        const uint64_t b3 = load_u64(b_packed + (g + 3u) * 8u);
        store_u64(product + (g + 0u) * 8u, clmul_pair(a0, b0));
        store_u64(product + (g + 1u) * 8u, clmul_pair(a1, b1));
        store_u64(product + (g + 2u) * 8u, clmul_pair(a2, b2));
        store_u64(product + (g + 3u) * 8u, clmul_pair(a3, b3));
    }
    for (; g < groups; ++g) {
        store_u64(product + g * 8u,
                  clmul_pair(load_u64(a_packed + g * 8u),
                             load_u64(b_packed + g * 8u)));
    }
}

void pqc_kernel_gf16_clmul_reduce_array(uint8_t *out_packed, const uint8_t *product,
                             size_t n) {
    const size_t groups = (n + 1u) / 2u;
    size_t g = 0;
    for (; g + 4u <= groups; g += 4u) {
        const uint64_t p0 = load_u64(product + (g + 0u) * 8u);
        const uint64_t p1 = load_u64(product + (g + 1u) * 8u);
        const uint64_t p2 = load_u64(product + (g + 2u) * 8u);
        const uint64_t p3 = load_u64(product + (g + 3u) * 8u);
        store_u64(out_packed + (g + 0u) * 8u, reduce_pair_product(p0));
        store_u64(out_packed + (g + 1u) * 8u, reduce_pair_product(p1));
        store_u64(out_packed + (g + 2u) * 8u, reduce_pair_product(p2));
        store_u64(out_packed + (g + 3u) * 8u, reduce_pair_product(p3));
    }
    for (; g < groups; ++g) {
        store_u64(out_packed + g * 8u,
                  reduce_pair_product(load_u64(product + g * 8u)));
    }
}

void pqc_kernel_gf16_clmul_mul_packed_array(uint8_t *out_packed, const uint8_t *a_packed,
                                 const uint8_t *b_packed, size_t n) {
    const size_t groups = (n + 1u) / 2u;
    size_t g = 0;
    for (; g + 4u <= groups; g += 4u) {
        const uint64_t a0 = load_u64(a_packed + (g + 0u) * 8u);
        const uint64_t a1 = load_u64(a_packed + (g + 1u) * 8u);
        const uint64_t a2 = load_u64(a_packed + (g + 2u) * 8u);
        const uint64_t a3 = load_u64(a_packed + (g + 3u) * 8u);
        const uint64_t b0 = load_u64(b_packed + (g + 0u) * 8u);
        const uint64_t b1 = load_u64(b_packed + (g + 1u) * 8u);
        const uint64_t b2 = load_u64(b_packed + (g + 2u) * 8u);
        const uint64_t b3 = load_u64(b_packed + (g + 3u) * 8u);
        store_u64(out_packed + (g + 0u) * 8u, reduce_pair_product(clmul_pair(a0, b0)));
        store_u64(out_packed + (g + 1u) * 8u, reduce_pair_product(clmul_pair(a1, b1)));
        store_u64(out_packed + (g + 2u) * 8u, reduce_pair_product(clmul_pair(a2, b2)));
        store_u64(out_packed + (g + 3u) * 8u, reduce_pair_product(clmul_pair(a3, b3)));
    }
    for (; g < groups; ++g) {
        store_u64(out_packed + g * 8u,
                  reduce_pair_product(clmul_pair(load_u64(a_packed + g * 8u),
                                                 load_u64(b_packed + g * 8u))));
    }
}

void pqc_kernel_gf16_clmul_ps_packed(uint8_t *out_packed, const uint64_t *coeff_words,
                          const uint8_t *sig_packed, size_t rows, size_t cols,
                          size_t lanes) {
    const size_t groups = (lanes + 1u) / 2u;
    memset(out_packed, 0, rows * groups * 8u);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out_packed + r * groups * 8u;
        for (size_t c = 0; c < cols; ++c) {
            const uint64_t coeff = coeff_words[r * cols + c];
            const uint8_t *src = sig_packed + c * groups * 8u;
            for (size_t g = 0; g < groups; ++g) {
                const uint64_t value = load_u64(dst + g * 8u) ^
                    reduce_pair_product(clmul_pair(coeff, load_u64(src + g * 8u)));
                store_u64(dst + g * 8u, value);
            }
        }
    }
}

void pqc_kernel_gf16_clmul_sts_packed(uint8_t *out_packed, const uint8_t *left_packed,
                           const uint8_t *right_packed, size_t rows,
                           size_t cols, size_t lanes) {
    const size_t groups = (lanes + 1u) / 2u;
    memset(out_packed, 0, rows * groups * 8u);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out_packed + r * groups * 8u;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *left = left_packed + (r * cols + c) * groups * 8u;
            const uint8_t *right = right_packed + c * groups * 8u;
            for (size_t g = 0; g < groups; ++g) {
                const uint64_t value = load_u64(dst + g * 8u) ^
                    reduce_pair_product(clmul_pair(load_u64(left + g * 8u),
                                                   load_u64(right + g * 8u)));
                store_u64(dst + g * 8u, value);
            }
        }
    }
}
