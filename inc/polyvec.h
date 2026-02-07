#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "poly.h"

// Forward declarations for Dilithium parameters
#ifndef DILI_K
#define DILI_K 4
#endif
#ifndef DILI_L
#define DILI_L 4
#endif
#ifndef DILI_N
#define DILI_N 256
#endif
#ifndef DILI_D
#define DILI_D 13
#endif
#ifndef DILI_Q
#define DILI_Q 8380417
#endif

typedef struct {
    poly vec[DILI_L];
} polyvecl;

typedef struct {
    poly vec[DILI_K];
} polyveck;

void polyvec_matrix_expand(polyvecl mat[DILI_K], const uint8_t rho[32], lp_ctx *ctx);
void polyvec_matrix_pointwise_montgomery(polyveck *t, const polyvecl mat[DILI_K], const polyvecl *v, lp_ctx *ctx);

void polyvecl_ntt(polyvecl *v, lp_ctx *ctx);
void polyvecl_invntt_tomont(polyvecl *v, lp_ctx *ctx);
void polyvecl_pointwise_acc_montgomery(poly *w, const polyvecl *u, const polyvecl *v, lp_ctx *ctx);
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v);
void polyvecl_copy(polyvecl *dst, const polyvecl *src);
void polyvecl_copy_from_array(polyvecl *dst, const poly *src);
void polyvecl_copy_to_array(poly *dst, const polyvecl *src);

void polyveck_ntt(polyveck *v, lp_ctx *ctx);
void polyveck_invntt_tomont(polyveck *v, lp_ctx *ctx);
void polyveck_caddq(polyveck *v);
void polyveck_reduce(polyveck *v);
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v);
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v);
void polyveck_shiftl(polyveck *v);
void polyveck_copy(polyveck *dst, const polyveck *src);
void polyveck_copy_from_array(polyveck *dst, const poly *src);
void polyveck_copy_to_array(poly *dst, const polyveck *src);

#endif
