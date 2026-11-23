#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "poly.h"

typedef void (*polyvec_matrix_sample_fn)(poly out, const uint8_t seed[34],
                                         lp_ctx *ctx, void *arg);

lssi_t poly_reduce_mod_one(lssi_t a, lssi_t q);
void poly_reduce_mod(lssi_t *a, int n, lssi_t q);
void poly_tomont(lssi_t *a, int n, lp_ctx *ctx);

void polyvec_tomont(poly *v, int len, lp_ctx *ctx);
void polyvec_ntt(poly *v, int len, lp_ctx *ctx);
void polyvec_invntt(poly *v, int len, lp_ctx *ctx);
void polyvec_invntt_tomont(poly *v, int len, lp_ctx *ctx);

void polyvec_add(poly *w, const poly *u, const poly *v, int len);
void polyvec_sub(poly *w, const poly *u, const poly *v, int len);
void polyvec_copy(poly *dst, const poly *src, int len);
void polyvec_reduce_mod(poly *v, int len, lssi_t q);
void polyvec_caddq_mod(poly *v, int len, lssi_t q);
void polyvec_shiftl(poly *v, int len, int shift);

/* Fixed 256-coefficient vectors. Transform products follow the plan and
 * poly.h input bounds; accumulation uses int64 and returns canonical [0,q).
 * Dimensions are nonnegative; matrix outputs must not overlap inputs. */
void polyvec_pointwise_acc(poly r, const poly *u, const poly *v, int len,
                           lp_ctx *ctx);
void polyvec_pointwise_acc_montgomery(poly *w, const poly *u, const poly *v,
                                      int len, lp_ctx *ctx);

void polyvec_matrix_pointwise(poly *t, const poly *mat, int rows, int cols,
                              const poly *v, lp_ctx *ctx);
void polyvec_matrix_pointwise_montgomery(poly *t, const poly *mat, int k,
                                         const poly *v, int l, lp_ctx *ctx);
void polyvec_matrix_transpose_pointwise(poly *t, const poly *mat, int rows,
                                        int cols, const poly *v, lp_ctx *ctx);

void polyvec_matrix_expand_with_sampler(poly *mat, int rows, int cols,
                                        const uint8_t rho[32], lp_ctx *ctx,
                                        polyvec_matrix_sample_fn sampler,
                                        void *arg);

#endif
