#include "polyvec.h"
#include "poly.h"
#include "keccak.h"
#include <string.h>

// Reduce helper macro
#define SGN_MASK(x) ((int32_t)(x)>>31)

#define REDUCE32(x) ({\
    int32_t _r = (x)%DILI_Q; \
    _r += (DILI_Q & SGN_MASK(_r)); \
    _r; \
})

void polyvec_matrix_expand(poly *mat, int k, int l, const uint8_t rho[32], lp_ctx *ctx) {
    unsigned int i, j;
    uint8_t seed[34];
    uint8_t buf[1024];
    
    for(i = 0; i < k; i++) {
        for(j = 0; j < l; j++) {
            memcpy(seed, rho, 32);
            seed[32] = (uint8_t)j;
            seed[33] = (uint8_t)i;
            shake128(buf, sizeof(buf), seed, 34);
            sample_uniform(mat[i * l + j], ctx->n, buf, sizeof(buf), ctx);
        }
    }
}

void polyvec_matrix_pointwise_montgomery(poly *t, const poly *mat, int k, const poly *v, int l, lp_ctx *ctx) {
    unsigned int i;
    
    for(i = 0; i < k; i++)
        polyvec_pointwise_acc_montgomery(&t[i], &mat[i * l], v, l, ctx);
}

void polyvec_pointwise_acc_montgomery(poly *w, const poly *u, const poly *v, int len, lp_ctx *ctx) {
    unsigned int i, j;
    poly t;
    
    poly_pointwise(*w, u[0], v[0], ctx);
    for(i = 1; i < len; i++) {
        poly_pointwise(t, u[i], v[i], ctx);
        for(j = 0; j < DILI_N; j++) {
            (*w)[j] += t[j];
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
            if(v[i][j] < 0)
                v[i][j] += DILI_Q;
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
