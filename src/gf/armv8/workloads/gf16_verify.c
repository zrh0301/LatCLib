#include "gf16_verify.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct pqc_kernel_gf16_verify_case {
    size_t l;
    size_t pairs;
    size_t rows;
    size_t cols;
    size_t lanes;
    size_t stream_out_bytes;
    size_t stream_in_bytes;
    size_t stream_left_bytes;
    size_t coeff_ps_bytes;
    size_t coeff_sec_bytes;

    uint8_t *coeff_ps;
    uint8_t *coeff_sec;
    uint8_t *state;
    uint8_t *left;
    uint8_t *right;

    uint8_t *ps;
    uint8_t *variable;
    uint8_t *secondary;
    uint8_t *output;
    uint8_t *reference;

    uint8_t *pair_coeff_ps;
    uint8_t *pair_coeff_sec;
    uint8_t *packed_state;
    uint8_t *packed_left;
    uint8_t *packed_right;
    uint8_t *packed_ps;
    uint8_t *packed_variable;
    uint8_t *packed_secondary;
    uint8_t *packed_output;
};

static uint64_t verify_now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

static void *verify_xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) {
        fprintf(stderr, "gf16 verify allocation failed: %zu bytes\n", n);
        abort();
    }
    return p;
}

static uint8_t verify_next(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return (uint8_t)((*state >> 28) & 0x0f);
}

static void verify_fill(uint8_t *p, size_t n, uint32_t *state) {
    for (size_t i = 0; i < n; ++i) {
        p[i] = verify_next(state);
    }
}

static void verify_pack_pair(uint8_t *dst, const uint8_t *lo,
                             const uint8_t *hi, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec l = pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(lo + i), 0x0f);
        pqc_kernel_gf16_vec h = pqc_kernel_gf16_vec_slli4(pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(hi + i), 0x0f));
        pqc_kernel_gf16_vec_store(dst + i, pqc_kernel_gf16_vec_or(l, h));
    }
    for (; i < n; ++i) {
        dst[i] = (uint8_t)((lo[i] & 0x0f) | ((hi[i] & 0x0f) << 4));
    }
}

static void verify_unpack_pair(uint8_t *lo, uint8_t *hi,
                               const uint8_t *src, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec v = pqc_kernel_gf16_vec_load(src + i);
        pqc_kernel_gf16_vec_store(lo + i, pqc_kernel_gf16_vec_mask_u8(v, 0x0f));
        pqc_kernel_gf16_vec_store(hi + i, pqc_kernel_gf16_vec_srli4(v));
    }
    for (; i < n; ++i) {
        lo[i] = src[i] & 0x0f;
        hi[i] = (src[i] >> 4) & 0x0f;
    }
}

static void verify_pack_inputs(struct pqc_kernel_gf16_verify_case *c) {
    for (size_t p = 0; p < c->pairs; ++p) {
        const size_t lo = 2u * p;
        const size_t hi = lo + 1u;
        verify_pack_pair(c->packed_state + p * c->stream_in_bytes,
                         c->state + lo * c->stream_in_bytes,
                         c->state + hi * c->stream_in_bytes,
                         c->stream_in_bytes);
        verify_pack_pair(c->packed_right + p * c->stream_in_bytes,
                         c->right + lo * c->stream_in_bytes,
                         c->right + hi * c->stream_in_bytes,
                         c->stream_in_bytes);
        verify_pack_pair(c->packed_left + p * c->stream_left_bytes,
                         c->left + lo * c->stream_left_bytes,
                         c->left + hi * c->stream_left_bytes,
                         c->stream_left_bytes);
    }
}

static void verify_unpack_outputs(struct pqc_kernel_gf16_verify_case *c) {
    for (size_t p = 0; p < c->pairs; ++p) {
        const size_t lo = 2u * p;
        const size_t hi = lo + 1u;
        verify_unpack_pair(c->output + lo * c->stream_out_bytes,
                           c->output + hi * c->stream_out_bytes,
                           c->packed_output + p * c->stream_out_bytes,
                           c->stream_out_bytes);
    }
}

static void verify_build_pair_coefficients(struct pqc_kernel_gf16_verify_case *c) {
    for (size_t p = 0; p < c->pairs; ++p) {
        const size_t lo = 2u * p;
        const size_t hi = lo + 1u;
        uint8_t *ps = c->pair_coeff_ps + p * 2u * c->coeff_ps_bytes;
        uint8_t *sec = c->pair_coeff_sec + p * 2u * c->coeff_sec_bytes;
        for (size_t i = 0; i < c->coeff_ps_bytes; ++i) {
            ps[2u * i] = c->coeff_ps[lo * c->coeff_ps_bytes + i];
            ps[2u * i + 1u] = c->coeff_ps[hi * c->coeff_ps_bytes + i];
        }
        for (size_t i = 0; i < c->coeff_sec_bytes; ++i) {
            sec[2u * i] = c->coeff_sec[lo * c->coeff_sec_bytes + i];
            sec[2u * i + 1u] = c->coeff_sec[hi * c->coeff_sec_bytes + i];
        }
    }
}

static void verify_xor3(uint8_t *dst, const uint8_t *a,
                        const uint8_t *b, const uint8_t *d, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        dst[i] = (uint8_t)(a[i] ^ b[i] ^ d[i]);
    }
}

static void verify_xor3_packed(uint8_t *dst, const uint8_t *a,
                               const uint8_t *b, const uint8_t *d, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i));
        r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_load(d + i));
        pqc_kernel_gf16_vec_store(dst + i, r);
    }
    for (; i < n; ++i) {
        dst[i] = (uint8_t)(a[i] ^ b[i] ^ d[i]);
    }
}

static void verify_reference_stream(struct pqc_kernel_gf16_verify_case *c, size_t s) {
    const uint8_t *coeff = c->coeff_ps + s * c->coeff_ps_bytes;
    const uint8_t *sec_coeff = c->coeff_sec + s * c->coeff_sec_bytes;
    const uint8_t *state = c->state + s * c->stream_in_bytes;
    const uint8_t *left = c->left + s * c->stream_left_bytes;
    const uint8_t *right = c->right + s * c->stream_in_bytes;
    uint8_t *ps = c->ps + s * c->stream_out_bytes;
    uint8_t *var = c->variable + s * c->stream_out_bytes;
    uint8_t *sec = c->secondary + s * c->stream_out_bytes;
    uint8_t *out = c->reference + s * c->stream_out_bytes;

    pqc_kernel_gf16_snova_ps_table(ps, coeff, state, c->rows, c->cols, c->lanes);
    pqc_kernel_gf16_snova_sts_table(var, left, right, c->rows, c->cols, c->lanes);
    pqc_kernel_gf16_snova_ps_table(sec, sec_coeff, var, c->rows, c->rows, c->lanes);
    verify_xor3(out, ps, var, sec, c->stream_out_bytes);
}

#if 0
static void verify_run_raw(struct pqc_kernel_gf16_verify_case *c, pqc_kernel_gf16_verify_method method) {
    for (size_t s = 0; s < c->l; ++s) {
        const uint8_t *coeff = c->coeff_ps + s * c->coeff_ps_bytes;
        const uint8_t *sec_coeff = c->coeff_sec + s * c->coeff_sec_bytes;
        const uint8_t *state = c->state + s * c->stream_in_bytes;
        const uint8_t *left = c->left + s * c->stream_left_bytes;
        const uint8_t *right = c->right + s * c->stream_in_bytes;
        uint8_t *ps = c->ps + s * c->stream_out_bytes;
        uint8_t *var = c->variable + s * c->stream_out_bytes;
        uint8_t *sec = c->secondary + s * c->stream_out_bytes;
        uint8_t *out = c->output + s * c->stream_out_bytes;

        if (method == GF16_VERIFY_TABLE) {
            pqc_kernel_gf16_snova_ps_table(ps, coeff, state, c->rows, c->cols, c->lanes);
            pqc_kernel_gf16_snova_sts_table(var, left, right, c->rows, c->cols, c->lanes);
            pqc_kernel_gf16_snova_ps_table(sec, sec_coeff, var, c->rows, c->rows, c->lanes);
        } else if (method == GF16_VERIFY_OFFICIAL_EXPANDED) {
            pqc_kernel_gf16_snova_ps_official(ps, coeff, state, c->rows, c->cols, c->lanes);
            pqc_kernel_gf16_snova_sts_official(var, left, right, c->rows, c->cols, c->lanes);
            pqc_kernel_gf16_snova_ps_official(sec, sec_coeff, var, c->rows, c->rows, c->lanes);
        } else {
            pqc_kernel_gf16_snova_ps_mtk2_pair(ps, coeff, state, c->rows, c->cols, c->lanes);
            pqc_kernel_gf16_snova_sts_bitplane_ilp(var, left, right, c->rows, c->cols, c->lanes);
            pqc_kernel_gf16_snova_ps_shuffle(sec, sec_coeff, var, c->rows, c->rows, c->lanes);
        }
        verify_xor3(out, ps, var, sec, c->stream_out_bytes);
    }
}
#endif

static void verify_sts_packed_pair_split(uint8_t *out_packed,
                                         const uint8_t *left_packed,
                                         const uint8_t *right_packed,
                                         size_t rows, size_t cols,
                                         size_t lanes) {
    memset(out_packed, 0, rows * lanes);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out_packed + r * lanes;
        const uint8_t *lrow = left_packed + r * cols * lanes;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc_lo = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec acc_hi = pqc_kernel_gf16_vec_zero();
            for (size_t col = 0; col < cols; ++col) {
                pqc_kernel_gf16_vec lo;
                pqc_kernel_gf16_vec hi;
                pqc_kernel_gf16_vec_mul_packed_ab_split_dual(
                    pqc_kernel_gf16_vec_load(lrow + col * lanes + i),
                    pqc_kernel_gf16_vec_load(right_packed + col * lanes + i),
                    &lo, &hi);
                acc_lo = pqc_kernel_gf16_vec_xor(acc_lo, lo);
                acc_hi = pqc_kernel_gf16_vec_xor(acc_hi, hi);
            }
            pqc_kernel_gf16_vec packed = pqc_kernel_gf16_vec_or(
                pqc_kernel_gf16_vec_mask_u8(acc_lo, 0x0f),
                pqc_kernel_gf16_vec_slli4(pqc_kernel_gf16_vec_mask_u8(acc_hi, 0x0f)));
            pqc_kernel_gf16_vec_store(dst + i, packed);
        }
        for (; i < lanes; ++i) {
            uint8_t lo = 0;
            uint8_t hi = 0;
            for (size_t col = 0; col < cols; ++col) {
                const uint8_t a = lrow[col * lanes + i];
                const uint8_t b = right_packed[col * lanes + i];
                lo ^= pqc_kernel_gf16_mul_table[a & 0x0f][b & 0x0f];
                hi ^= pqc_kernel_gf16_mul_table[(a >> 4) & 0x0f][(b >> 4) & 0x0f];
            }
            dst[i] = (uint8_t)((lo & 0x0f) | ((hi & 0x0f) << 4));
        }
    }
}

#if 0
static void verify_run_packed(struct pqc_kernel_gf16_verify_case *c, int include_conversion) {
    if (include_conversion) {
        verify_pack_inputs(c);
    }
    for (size_t p = 0; p < c->pairs; ++p) {
        const uint8_t *ps_coeff = c->pair_coeff_ps + p * 2u * c->coeff_ps_bytes;
        const uint8_t *sec_coeff = c->pair_coeff_sec + p * 2u * c->coeff_sec_bytes;
        uint8_t *ps = c->packed_ps + p * c->stream_out_bytes;
        uint8_t *var = c->packed_variable + p * c->stream_out_bytes;
        uint8_t *sec = c->packed_secondary + p * c->stream_out_bytes;
        const uint8_t *state = c->packed_state + p * c->stream_in_bytes;
        const uint8_t *left = c->packed_left + p * c->stream_left_bytes;
        const uint8_t *right = c->packed_right + p * c->stream_in_bytes;
        pqc_kernel_gf16_snova_packed_pair_contract(ps, ps_coeff, state,
                                         c->rows, c->cols, c->lanes);
        verify_sts_packed_pair_split(var, left, right,
                                     c->rows, c->cols, c->lanes);
        pqc_kernel_gf16_snova_packed_pair_contract(sec, sec_coeff, var,
                                         c->rows, c->rows, c->lanes);
        verify_xor3_packed(c->packed_output, ps, var, sec,
                           c->stream_out_bytes);
    }
    if (include_conversion) {
        verify_unpack_outputs(c);
    } else {
        /* Compute-only still exposes a raw output for correctness checks,
         * but conversion is intentionally outside the timed region. */
        verify_unpack_outputs(c);
    }
}
#endif

static void zero_metrics(pqc_kernel_gf16_verify_metrics *m) {
    if (m) {
        memset(m, 0, sizeof(*m));
    }
}

static void measure_raw_stage(struct pqc_kernel_gf16_verify_case *c,
                             pqc_kernel_gf16_verify_method method,
                             pqc_kernel_gf16_verify_metrics *m) {
    uint64_t t = 0;
    if (m) t = verify_now_ns();
    for (size_t s = 0; s < c->l; ++s) {
        const uint8_t *coeff = c->coeff_ps + s * c->coeff_ps_bytes;
        const uint8_t *state = c->state + s * c->stream_in_bytes;
        uint8_t *ps = c->ps + s * c->stream_out_bytes;
        if (method == GF16_VERIFY_TABLE) {
            pqc_kernel_gf16_snova_ps_table(ps, coeff, state, c->rows, c->cols, c->lanes);
        } else if (method == GF16_VERIFY_OFFICIAL_EXPANDED) {
            pqc_kernel_gf16_snova_ps_official(ps, coeff, state, c->rows, c->cols, c->lanes);
        } else {
            pqc_kernel_gf16_snova_ps_shuffle(ps, coeff, state, c->rows, c->cols, c->lanes);
        }
    }
    if (m) m->ps_ns += verify_now_ns() - t;

    if (m) t = verify_now_ns();
    for (size_t s = 0; s < c->l; ++s) {
        const uint8_t *left = c->left + s * c->stream_left_bytes;
        const uint8_t *right = c->right + s * c->stream_in_bytes;
        uint8_t *var = c->variable + s * c->stream_out_bytes;
        if (method == GF16_VERIFY_TABLE) {
            pqc_kernel_gf16_snova_sts_table(var, left, right, c->rows, c->cols, c->lanes);
        } else if (method == GF16_VERIFY_OFFICIAL_EXPANDED) {
            pqc_kernel_gf16_snova_sts_official(var, left, right, c->rows, c->cols, c->lanes);
        } else {
            pqc_kernel_gf16_snova_sts_same_b4(var, left, right, c->rows, c->cols, c->lanes);
        }
    }
    if (m) m->variable_ns += verify_now_ns() - t;

    if (m) t = verify_now_ns();
    for (size_t s = 0; s < c->l; ++s) {
        const uint8_t *coeff = c->coeff_sec + s * c->coeff_sec_bytes;
        const uint8_t *var = c->variable + s * c->stream_out_bytes;
        uint8_t *sec = c->secondary + s * c->stream_out_bytes;
        if (method == GF16_VERIFY_TABLE) {
            pqc_kernel_gf16_snova_ps_table(sec, coeff, var, c->rows, c->rows, c->lanes);
        } else if (method == GF16_VERIFY_OFFICIAL_EXPANDED) {
            pqc_kernel_gf16_snova_ps_official(sec, coeff, var, c->rows, c->rows, c->lanes);
        } else {
            pqc_kernel_gf16_snova_ps_mtk2_pair(sec, coeff, var, c->rows, c->rows, c->lanes);
        }
    }
    if (m) m->secondary_ns += verify_now_ns() - t;

    for (size_t s = 0; s < c->l; ++s) {
        verify_xor3(c->output + s * c->stream_out_bytes,
                    c->ps + s * c->stream_out_bytes,
                    c->variable + s * c->stream_out_bytes,
                    c->secondary + s * c->stream_out_bytes,
                    c->stream_out_bytes);
    }
}

pqc_kernel_gf16_verify_case *pqc_kernel_gf16_verify_create(size_t l, size_t rows, size_t cols,
                                     size_t lanes, uint32_t seed) {
    if ((l != 2u && l != 4u) || rows == 0 || cols == 0 || lanes == 0) {
        return NULL;
    }
    pqc_kernel_gf16_tables_init();
    pqc_kernel_gf16_verify_case *c = (pqc_kernel_gf16_verify_case *)calloc(1, sizeof(*c));
    if (!c) return NULL;
    c->l = l;
    c->pairs = l / 2u;
    c->rows = rows;
    c->cols = cols;
    c->lanes = lanes;
    c->stream_out_bytes = rows * lanes;
    c->stream_in_bytes = cols * lanes;
    c->stream_left_bytes = rows * cols * lanes;
    c->coeff_ps_bytes = rows * cols;
    c->coeff_sec_bytes = rows * rows;

    const size_t streams = l;
    const size_t pairs = c->pairs;
    c->coeff_ps = (uint8_t *)verify_xmalloc(streams * c->coeff_ps_bytes);
    c->coeff_sec = (uint8_t *)verify_xmalloc(streams * c->coeff_sec_bytes);
    c->state = (uint8_t *)verify_xmalloc(streams * c->stream_in_bytes);
    c->left = (uint8_t *)verify_xmalloc(streams * c->stream_left_bytes);
    c->right = (uint8_t *)verify_xmalloc(streams * c->stream_in_bytes);
    c->ps = (uint8_t *)verify_xmalloc(streams * c->stream_out_bytes);
    c->variable = (uint8_t *)verify_xmalloc(streams * c->stream_out_bytes);
    c->secondary = (uint8_t *)verify_xmalloc(streams * c->stream_out_bytes);
    c->output = (uint8_t *)verify_xmalloc(streams * c->stream_out_bytes);
    c->reference = (uint8_t *)verify_xmalloc(streams * c->stream_out_bytes);

    c->pair_coeff_ps = (uint8_t *)verify_xmalloc(pairs * 2u * c->coeff_ps_bytes);
    c->pair_coeff_sec = (uint8_t *)verify_xmalloc(pairs * 2u * c->coeff_sec_bytes);
    c->packed_state = (uint8_t *)verify_xmalloc(pairs * c->stream_in_bytes);
    c->packed_left = (uint8_t *)verify_xmalloc(pairs * c->stream_left_bytes);
    c->packed_right = (uint8_t *)verify_xmalloc(pairs * c->stream_in_bytes);
    c->packed_ps = (uint8_t *)verify_xmalloc(pairs * c->stream_out_bytes);
    c->packed_variable = (uint8_t *)verify_xmalloc(pairs * c->stream_out_bytes);
    c->packed_secondary = (uint8_t *)verify_xmalloc(pairs * c->stream_out_bytes);
    c->packed_output = (uint8_t *)verify_xmalloc(pairs * c->stream_out_bytes);

    verify_fill(c->coeff_ps, streams * c->coeff_ps_bytes, &seed);
    verify_fill(c->coeff_sec, streams * c->coeff_sec_bytes, &seed);
    verify_fill(c->state, streams * c->stream_in_bytes, &seed);
    verify_fill(c->left, streams * c->stream_left_bytes, &seed);
    verify_fill(c->right, streams * c->stream_in_bytes, &seed);
    verify_build_pair_coefficients(c);
    verify_pack_inputs(c);

    for (size_t s = 0; s < l; ++s) {
        verify_reference_stream(c, s);
    }
    return c;
}

void pqc_kernel_gf16_verify_destroy(pqc_kernel_gf16_verify_case *c) {
    if (!c) return;
    free(c->coeff_ps);
    free(c->coeff_sec);
    free(c->state);
    free(c->left);
    free(c->right);
    free(c->ps);
    free(c->variable);
    free(c->secondary);
    free(c->output);
    free(c->reference);
    free(c->pair_coeff_ps);
    free(c->pair_coeff_sec);
    free(c->packed_state);
    free(c->packed_left);
    free(c->packed_right);
    free(c->packed_ps);
    free(c->packed_variable);
    free(c->packed_secondary);
    free(c->packed_output);
    free(c);
}

int pqc_kernel_gf16_verify_run(pqc_kernel_gf16_verify_case *c, pqc_kernel_gf16_verify_method method,
                    pqc_kernel_gf16_verify_metrics *metrics) {
    if (!c || method < GF16_VERIFY_TABLE ||
        method > GF16_VERIFY_PACKED_SPLIT_COMPUTE) {
        return -1;
    }
    zero_metrics(metrics);
    const uint64_t total_start = metrics ? verify_now_ns() : 0;

    if (method == GF16_VERIFY_PACKED_SPLIT_FULL ||
        method == GF16_VERIFY_PACKED_SPLIT_COMPUTE) {
        const int full = method == GF16_VERIFY_PACKED_SPLIT_FULL;
        uint64_t t = 0;
        if (full) {
            if (metrics) t = verify_now_ns();
            verify_pack_inputs(c);
            if (metrics) metrics->pack_ns = verify_now_ns() - t;
        }

        if (metrics) t = verify_now_ns();
        for (size_t p = 0; p < c->pairs; ++p) {
            const uint8_t *coeff = c->pair_coeff_ps + p * 2u * c->coeff_ps_bytes;
            pqc_kernel_gf16_snova_packed_pair_contract(c->packed_ps + p * c->stream_out_bytes,
                                             coeff,
                                             c->packed_state + p * c->stream_in_bytes,
                                             c->rows, c->cols, c->lanes);
        }
        if (metrics) metrics->ps_ns = verify_now_ns() - t;

        if (metrics) t = verify_now_ns();
        for (size_t p = 0; p < c->pairs; ++p) {
            verify_sts_packed_pair_split(
                c->packed_variable + p * c->stream_out_bytes,
                c->packed_left + p * c->stream_left_bytes,
                c->packed_right + p * c->stream_in_bytes,
                c->rows, c->cols, c->lanes);
        }
        if (metrics) metrics->variable_ns = verify_now_ns() - t;

        if (metrics) t = verify_now_ns();
        for (size_t p = 0; p < c->pairs; ++p) {
            const uint8_t *coeff = c->pair_coeff_sec + p * 2u * c->coeff_sec_bytes;
            pqc_kernel_gf16_snova_packed_pair_contract(c->packed_secondary + p * c->stream_out_bytes,
                                             coeff,
                                             c->packed_variable + p * c->stream_out_bytes,
                                             c->rows, c->rows, c->lanes);
            verify_xor3_packed(c->packed_output + p * c->stream_out_bytes,
                               c->packed_ps + p * c->stream_out_bytes,
                               c->packed_variable + p * c->stream_out_bytes,
                               c->packed_secondary + p * c->stream_out_bytes,
                               c->stream_out_bytes);
        }
        if (metrics) metrics->secondary_ns = verify_now_ns() - t;

        if (full) {
            if (metrics) t = verify_now_ns();
            verify_unpack_outputs(c);
            if (metrics) metrics->unpack_ns = verify_now_ns() - t;
        }
    } else {
        measure_raw_stage(c, method, metrics);
    }

    if (metrics) {
        metrics->total_ns = verify_now_ns() - total_start;
    }
    return 0;
}

int pqc_kernel_gf16_verify_check(pqc_kernel_gf16_verify_case *c, pqc_kernel_gf16_verify_method method) {
    if (!c) return 0;
    if (method == GF16_VERIFY_PACKED_SPLIT_COMPUTE) {
        verify_unpack_outputs(c);
    }
    if (!c->output || !c->reference) return 0;
    return memcmp(c->output, c->reference,
                  c->l * c->stream_out_bytes) == 0;
}

const char *pqc_kernel_gf16_verify_method_name(pqc_kernel_gf16_verify_method method) {
    switch (method) {
    case GF16_VERIFY_TABLE: return "scalar-table";
    case GF16_VERIFY_OFFICIAL_EXPANDED: return "official-expanded";
    case GF16_VERIFY_SHUFFLE_BITPLANE: return "mtk2+same-b4";
    case GF16_VERIFY_PACKED_SPLIT_FULL: return "packed-a/b-split-full";
    case GF16_VERIFY_PACKED_SPLIT_COMPUTE: return "packed-a/b-split-compute";
    default: return "unknown";
    }
}

const uint8_t *pqc_kernel_gf16_verify_output(const pqc_kernel_gf16_verify_case *c) {
    return c ? c->output : NULL;
}

size_t pqc_kernel_gf16_verify_output_bytes(const pqc_kernel_gf16_verify_case *c) {
    return c ? c->l * c->stream_out_bytes : 0;
}

size_t pqc_kernel_gf16_verify_l(const pqc_kernel_gf16_verify_case *c) {
    return c ? c->l : 0;
}
