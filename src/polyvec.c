#include "polyvec.h"
#include <string.h>

lssi_t poly_reduce_mod_one(lssi_t a, lssi_t q) {
    return fp32_reduce_raw(a, q);
}

void poly_reduce_mod(lssi_t *a, int n, lssi_t q) {
    for (int i = 0; i < n; ++i) a[i] = fp32_reduce_raw(a[i], q);
}

void poly_tomont(lssi_t *a, int n, lp_ctx *ctx) {
    for(int i = 0; i < n; i++) {
        a[i] = fp32_to_mont(a[i], &ctx->ring.field);
    }
}

void polyvec_tomont(poly *v, int len, lp_ctx *ctx) {
    for(int i = 0; i < len; i++) {
        poly_tomont(v[i], ctx->ring.n, ctx);
    }
}

void polyvec_ntt(poly *v, int len, lp_ctx *ctx) {
    for(int i = 0; i < len; i++) {
        poly_ntt(v[i], ctx);
    }
}

void polyvec_invntt(poly *v, int len, lp_ctx *ctx) {
    for(int i = 0; i < len; i++) {
        poly_invntt(v[i], ctx);
    }
}

void polyvec_invntt_tomont(poly *v, int len, lp_ctx *ctx) {
    polyvec_invntt(v, len, ctx);
}

void polyvec_add(poly *w, const poly *u, const poly *v, int len) {
    for(int i = 0; i < len; i++) {
        POLY_ADD(w[i], u[i], v[i], MAX_POLY_COEFS);
    }
}

void polyvec_sub(poly *w, const poly *u, const poly *v, int len) {
    for(int i = 0; i < len; i++) {
        POLY_SUB(w[i], u[i], v[i], MAX_POLY_COEFS);
    }
}

void polyvec_copy(poly *dst, const poly *src, int len) {
    for(int i = 0; i < len; i++) {
        memcpy(dst[i], src[i], MAX_POLY_COEFS * sizeof(lssi_t));
    }
}

void polyvec_reduce_mod(poly *v, int len, lssi_t q) {
    for(int i = 0; i < len; i++) {
        poly_reduce_mod(v[i], MAX_POLY_COEFS, q);
    }
}

void polyvec_caddq_mod(poly *v, int len, lssi_t q) {
    polyvec_reduce_mod(v, len, q);
}

void polyvec_shiftl(poly *v, int len, int shift) {
    for(int i = 0; i < len; i++) {
        POLY_LSH_LOW(v[i], v[i], shift, MAX_POLY_COEFS);
    }
}

void polyvec_pointwise_acc(poly r, const poly *u, const poly *v, int len,
                           lp_ctx *ctx) {
    poly tmp;
    int64_t sum[MAX_POLY_COEFS] = {0};

    for (int i = 0; i < len; ++i) {
        if (ctx->block_degree == 1) {
            poly_pointwise(tmp, u[i], v[i], ctx);
        } else {
            poly_base_pointwise(tmp, u[i], v[i], ctx);
        }
        for (int j = 0; j < ctx->ring.n; ++j) sum[j] += tmp[j];
    }
    for (int j = 0; j < ctx->ring.n; ++j) {
        r[j] = fp32_reduce(sum[j], &ctx->ring.field);
    }
}

void polyvec_pointwise_acc_montgomery(poly *w, const poly *u, const poly *v,
                                      int len, lp_ctx *ctx) {
    polyvec_pointwise_acc(*w, u, v, len, ctx);
}

void polyvec_matrix_pointwise(poly *t, const poly *mat, int rows, int cols,
                              const poly *v, lp_ctx *ctx) {
    for(int i = 0; i < rows; i++) {
        polyvec_pointwise_acc(t[i], &mat[i * cols], v, cols, ctx);
    }
}

void polyvec_matrix_pointwise_montgomery(poly *t, const poly *mat, int k,
                                         const poly *v, int l, lp_ctx *ctx) {
    polyvec_matrix_pointwise(t, mat, k, l, v, ctx);
}

void polyvec_matrix_transpose_pointwise(poly *t, const poly *mat, int rows,
                                        int cols, const poly *v, lp_ctx *ctx) {
    poly tmp;
    for (int j = 0; j < cols; ++j) {
        int64_t sum[MAX_POLY_COEFS] = {0};
        for (int i = 0; i < rows; ++i) {
            if (ctx->block_degree == 1) {
                poly_pointwise(tmp, mat[i * cols + j], v[i], ctx);
            } else {
                poly_base_pointwise(tmp, mat[i * cols + j], v[i], ctx);
            }
            for (int k = 0; k < ctx->ring.n; ++k) sum[k] += tmp[k];
        }
        for (int k = 0; k < ctx->ring.n; ++k) {
            t[j][k] = fp32_reduce(sum[k], &ctx->ring.field);
        }
    }
}

void polyvec_matrix_expand_with_sampler(poly *mat, int rows, int cols,
                                        const uint8_t rho[32], lp_ctx *ctx,
                                        polyvec_matrix_sample_fn sampler,
                                        void *arg) {
    uint8_t seed[34];

    memcpy(seed, rho, 32);
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            seed[32] = (uint8_t)j;
            seed[33] = (uint8_t)i;
            sampler(mat[i * cols + j], seed, ctx, arg);
        }
    }
}
