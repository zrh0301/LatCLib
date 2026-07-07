#include "polyvec.h"
#include "poly.h"
#include "keccak.h"
#include <string.h>

// Reduce helper macro
#define SGN_MASK(x) ((int32_t)(x)>>31)

#define REDUCE32(x) ({\
    int32_t _a = (x); \
    int32_t _t = (_a + (1<<22)) >> 23; \
    _t = _a - _t*DILI_Q; \
    _t += (DILI_Q & SGN_MASK(_t)); \
    _t; \
})

void polyvec_matrix_expand(poly *mat, int k, int l, const uint8_t rho[32], lp_ctx *ctx) {
    unsigned int i, j;
    uint8_t seed[34];
    uint8_t buf[840];   // 5 x SHAKE128 rate(168); 840=3*280, 无 3 字节端数
    KECCAK1600_CTX kc;
    int got;

    // 生成矩阵 A
    for(i = 0; i < k; i++) {
        for(j = 0; j < l; j++) {
            memcpy(seed, rho, 32);
            seed[32] = (uint8_t)j;
            seed[33] = (uint8_t)i;

            shake128_init(&kc);
            sha3_update(&kc, seed, 34);
            shake_xof_final(&kc);

            // 拒绝采样：持续 squeeze 直到填满全部系数，避免残留未初始化值
            got = 0;
            while(got < ctx->n) {
                shake_squeeze(&kc, buf, sizeof(buf));
                got += sample_uniform(mat[i * l + j] + got, ctx->n - got,
                                      buf, sizeof(buf), ctx);
            }
        }
    }

    for(i = 0; i < k * l; i++) {
        poly_ntt(mat[i], ctx);
    }
}

void polyvec_matrix_pointwise_montgomery(poly *t, const poly *mat, int k, const poly *v, int l, lp_ctx *ctx) {
    unsigned int i;
    
    for(i = 0; i < k; i++)
        polyvec_pointwise_acc_montgomery(&t[i], &mat[i * l], v, l, ctx);
}

void polyvec_pointwise_acc_montgomery(poly *w, const poly *u, const poly *v, int len, lp_ctx *ctx) {
    unsigned int i, j;
    int32_t q = ctx->q;
    int32_t qc = ctx->qc;
    
    for(j = 0; j < DILI_N; j++) {
        (*w)[j] = MONTGOMERY_REDUCE((lsdi_t)u[0][j] * v[0][j], q, qc);
    }

    for(i = 1; i < len; i++) {
        for(j = 0; j < DILI_N; j++) {
            (*w)[j] += MONTGOMERY_REDUCE((lsdi_t)u[i][j] * v[i][j], q, qc);
        }
    }
}

void polyvec_ntt(poly *v, int len, lp_ctx *ctx) {
    ntt_vector(v, len, ctx);
}

void polyvec_invntt_tomont(poly *v, int len, lp_ctx *ctx) {
    unsigned int i;
    
    for(i = 0; i < len; i++)
        poly_invntt(v[i], ctx);
}

void polyvec_add(poly *w, const poly *u, const poly *v, int len) {
    unsigned int i;
    
    for(i = 0; i < len; i++)
        POLY_ADD(w[i], u[i], v[i], DILI_N);
}

void polyvec_sub(poly *w, const poly *u, const poly *v, int len) {
    unsigned int i;
    
    for(i = 0; i < len; i++)
        POLY_SUB(w[i], u[i], v[i], DILI_N);
}

void polyvec_caddq(poly *v, int len) {
    unsigned int i, j;

    for(i = 0; i < len; i++) {
        for(j = 0; j < DILI_N; j++) {
            v[i][j] = REDUCE32(v[i][j]);
        }
    }
}

void polyvec_reduce(poly *v, int len) {
    unsigned int i, j;
    
    for(i = 0; i < len; i++) {
        for(j = 0; j < DILI_N; j++) {
            v[i][j] = REDUCE32(v[i][j]);
        }
    }
}

void polyvec_shiftl(poly *v, int len) {
    unsigned int i, j;
    
    for(i = 0; i < len; i++) {
        for(j = 0; j < DILI_N; j++) {
            v[i][j] <<= DILI_D;
        }
    }
}

void polyvec_copy(poly *dst, const poly *src, int len) {
    unsigned int i;
    
    for(i = 0; i < len; i++)
        memcpy(dst[i], src[i], DILI_N * sizeof(lssi_t));
}
