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

void polyvec_matrix_expand(polyvecl mat[DILI_K], const uint8_t rho[32], lp_ctx *ctx) {
    unsigned int i, j;
    uint8_t seed[34];
    uint8_t buf[1024];
    
    for(i = 0; i < DILI_K; i++) {
        for(j = 0; j < DILI_L; j++) {
            memcpy(seed, rho, 32);
            seed[32] = (uint8_t)j;
            seed[33] = (uint8_t)i;
            shake128(buf, sizeof(buf), seed, 34);
            sample_uniform(mat[i].vec[j], ctx->n, buf, sizeof(buf), ctx);
        }
    }
}

void polyvec_matrix_pointwise_montgomery(polyveck *t, const polyvecl mat[DILI_K], const polyvecl *v, lp_ctx *ctx) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        polyvecl_pointwise_acc_montgomery(&t->vec[i], &mat[i], v, ctx);
}

void polyvecl_pointwise_acc_montgomery(poly *w, const polyvecl *u, const polyvecl *v, lp_ctx *ctx) {
    unsigned int i, j;
    poly t;
    
    poly_pointwise(*w, u->vec[0], v->vec[0], ctx);
    for(i = 1; i < DILI_L; i++) {
        poly_pointwise(t, u->vec[i], v->vec[i], ctx);
        for(j = 0; j < DILI_N; j++) {
            (*w)[j] += t[j];
        }
    }
}

void polyvecl_ntt(polyvecl *v, lp_ctx *ctx) {
    ntt_vector(v->vec, DILI_L, ctx);
}

void polyvecl_invntt_tomont(polyvecl *v, lp_ctx *ctx) {
    unsigned int i;
    
    for(i = 0; i < DILI_L; i++)
        poly_invntt(v->vec[i], ctx);
}

void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v) {
    unsigned int i;
    
    for(i = 0; i < DILI_L; i++)
        POLY_ADD(w->vec[i], u->vec[i], v->vec[i], DILI_N);
}

void polyveck_ntt(polyveck *v, lp_ctx *ctx) {
    ntt_vector(v->vec, DILI_K, ctx);
}

void polyveck_invntt_tomont(polyveck *v, lp_ctx *ctx) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        poly_invntt(v->vec[i], ctx);
}

void polyveck_caddq(polyveck *v) {
    unsigned int i, j;
    
    for(i = 0; i < DILI_K; i++) {
        for(j = 0; j < DILI_N; j++) {
            v->vec[i][j] = REDUCE32(v->vec[i][j]);
            if(v->vec[i][j] < 0)
                v->vec[i][j] += DILI_Q;
        }
    }
}

void polyveck_reduce(polyveck *v) {
    unsigned int i, j;
    
    for(i = 0; i < DILI_K; i++) {
        for(j = 0; j < DILI_N; j++) {
            v->vec[i][j] = REDUCE32(v->vec[i][j]);
        }
    }
}

void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        POLY_ADD(w->vec[i], u->vec[i], v->vec[i], DILI_N);
}

void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        POLY_SUB(w->vec[i], u->vec[i], v->vec[i], DILI_N);
}

void polyveck_shiftl(polyveck *v) {
    unsigned int i, j;
    
    for(i = 0; i < DILI_K; i++) {
        for(j = 0; j < DILI_N; j++) {
            v->vec[i][j] <<= DILI_D;
        }
    }
}

void polyvecl_copy(polyvecl *dst, const polyvecl *src) {
    unsigned int i;
    
    for(i = 0; i < DILI_L; i++)
        memcpy(dst->vec[i], src->vec[i], DILI_N * sizeof(lssi_t));
}

void polyvecl_copy_from_array(polyvecl *dst, const poly *src) {
    unsigned int i;
    
    for(i = 0; i < DILI_L; i++)
        memcpy(dst->vec[i], src[i], DILI_N * sizeof(lssi_t));
}

void polyvecl_copy_to_array(poly *dst, const polyvecl *src) {
    unsigned int i;
    
    for(i = 0; i < DILI_L; i++)
        memcpy(dst[i], src->vec[i], DILI_N * sizeof(lssi_t));
}

void polyveck_copy(polyveck *dst, const polyveck *src) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        memcpy(dst->vec[i], src->vec[i], DILI_N * sizeof(lssi_t));
}

void polyveck_copy_from_array(polyveck *dst, const poly *src) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        memcpy(dst->vec[i], src[i], DILI_N * sizeof(lssi_t));
}

void polyveck_copy_to_array(poly *dst, const polyveck *src) {
    unsigned int i;
    
    for(i = 0; i < DILI_K; i++)
        memcpy(dst[i], src->vec[i], DILI_N * sizeof(lssi_t));
}
