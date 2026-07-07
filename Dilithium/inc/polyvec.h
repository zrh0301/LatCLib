#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "poly.h"
#include "dili.h"


// Matrix operations
void polyvec_matrix_expand(poly *mat, int k, int l, const uint8_t rho[32], lp_ctx *ctx);
void polyvec_matrix_pointwise_montgomery(poly *t, const poly *mat, int k, const poly *v, int l, lp_ctx *ctx);

// Vector operations (统一接口，通过size参数区分)
void polyvec_ntt(poly *v, int len, lp_ctx *ctx);
void polyvec_invntt_tomont(poly *v, int len, lp_ctx *ctx);
void polyvec_pointwise_acc_montgomery(poly *w, const poly *u, const poly *v, int len, lp_ctx *ctx);
void polyvec_add(poly *w, const poly *u, const poly *v, int len);
void polyvec_sub(poly *w, const poly *u, const poly *v, int len);
void polyvec_copy(poly *dst, const poly *src, int len);
void polyvec_caddq(poly *v, int len);
void polyvec_reduce(poly *v, int len);
void polyvec_shiftl(poly *v, int len);

#endif