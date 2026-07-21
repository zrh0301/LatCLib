#include "polyvec.h"
#include <string.h>

#define DILI_Q_VALUE 8380417

static inline lssi_t dili_reduce32(lssi_t a) {
    lssi_t t = (lssi_t)(((lsdi_t)a + (1LL << 22)) >> 23);
    return (lssi_t)((lsdi_t)a - (lsdi_t)t * DILI_Q_VALUE);
}

static inline lssi_t dili_freeze(lssi_t a) {
    a = dili_reduce32(a);
    a += (a >> 31) & DILI_Q_VALUE;
    return a;
}

lssi_t poly_reduce_mod_one(lssi_t a, lssi_t q) {
    if(q == DILI_Q_VALUE) {
        return dili_freeze(a);
    }
    a %= q;
    if(a < 0) {
        a += q;
    }
    return a;
}

void poly_reduce_mod(lssi_t *a, int n, lssi_t q) {
    if(q == DILI_Q_VALUE) {
        for(int i = 0; i < n; i++) {
            a[i] = dili_freeze(a[i]);
        }
        return;
    }
    for(int i = 0; i < n; i++) {
        a[i] = poly_reduce_mod_one(a[i], q);
    }
}

void poly_tomont(lssi_t *a, int n, lp_ctx *ctx) {
    for(int i = 0; i < n; i++) {
        a[i] = (lssi_t)TO_MONT(poly_reduce_mod_one(a[i], ctx->q), ctx->q);
    }
}

void polyvec_tomont(poly *v, int len, lp_ctx *ctx) {
    for(int i = 0; i < len; i++) {
        poly_tomont(v[i], ctx->n, ctx);
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
                           lp_ctx *ctx, polyvec_mul_mode_t mode) {
    poly tmp;

    if(len <= 0) {
        POLY_ZERO(r, ctx->n);
        return;
    }

    if(mode == POLYVEC_MUL_FULL) {
        for(int j = 0; j < ctx->n; j++) {
            r[j] = MONTGOMERY_REDUCE((lsdi_t)u[0][j] * v[0][j],
                                     ctx->q, ctx->qc);
        }
        for(int i = 1; i < len; i++) {
            for(int j = 0; j < ctx->n; j++) {
                r[j] += MONTGOMERY_REDUCE((lsdi_t)u[i][j] * v[i][j],
                                          ctx->q, ctx->qc);
            }
        }
        return;
    }

    POLY_ZERO(r, ctx->n);
    for(int i = 0; i < len; i++) {
        poly_base_pointwise(tmp, u[i], v[i], ctx);
        POLY_ADD(r, r, tmp, ctx->n);
    }
    poly_reduce_mod(r, ctx->n, ctx->q);
}

void polyvec_pointwise_acc_montgomery(poly *w, const poly *u, const poly *v,
                                      int len, lp_ctx *ctx) {
    polyvec_pointwise_acc(*w, u, v, len, ctx, POLYVEC_MUL_FULL);
}

void polyvec_matrix_pointwise(poly *t, const poly *mat, int rows, int cols,
                              const poly *v, lp_ctx *ctx,
                              polyvec_mul_mode_t mode) {
    for(int i = 0; i < rows; i++) {
        polyvec_pointwise_acc(t[i], &mat[i * cols], v, cols, ctx, mode);
    }
}

void polyvec_matrix_pointwise_montgomery(poly *t, const poly *mat, int k,
                                         const poly *v, int l, lp_ctx *ctx) {
    polyvec_matrix_pointwise(t, mat, k, l, v, ctx, POLYVEC_MUL_FULL);
}

void polyvec_matrix_transpose_pointwise(poly *t, const poly *mat, int rows,
                                        int cols, const poly *v, lp_ctx *ctx,
                                        polyvec_mul_mode_t mode) {
    poly tmp;

    for(int j = 0; j < cols; j++) {
        POLY_ZERO(t[j], ctx->n);
        for(int i = 0; i < rows; i++) {
            if(mode == POLYVEC_MUL_FULL) {
                poly_pointwise(tmp, mat[i * cols + j], v[i], ctx);
            } else {
                poly_base_pointwise(tmp, mat[i * cols + j], v[i], ctx);
            }
            POLY_ADD(t[j], t[j], tmp, ctx->n);
        }
        if(mode == POLYVEC_MUL_BASE) {
            poly_reduce_mod(t[j], ctx->n, ctx->q);
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
