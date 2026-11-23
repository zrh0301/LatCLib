#ifndef PQC_NTT_H
#define PQC_NTT_H

#include <stddef.h>
#include "ring.h"

/* This first backend is explicitly limited to n=256. The plan owns its ring
 * by value; zetas points to immutable, process-lifetime tables. Copying a plan
 * is safe. Do not change fields (including layers) after initialization. */
typedef struct {
    ring_ctx ring;
    int layers;
    int block_degree;        /* 1: scalar product, 2: quadratic base product. */
    int32_t inverse_scale;  /* R^2 / 2^layers modulo q, R=2^32. */
    const int32_t *zetas;   /* centered Montgomery roots, bit-reversed order */
    size_t zetas_count;
} ntt_plan;

/* Supported public configurations: (q,n,layers) = (8380417,256,8),
 * (8380417,256,7), (3329,256,7). No implicit fallback or arbitrary roots.
 * Return 0 on success, -1 for unsupported/invalid input; leave out unchanged.
 * Forward output uses bit-reversed blocks of block_degree coefficients.
 * Forward preserves Montgomery exponent; mul_ntt subtracts one exponent;
 * inverse adds one exponent. Thus invntt(ntt(a)) == R*a (mod q), while
 * invntt(mul_ntt(ntt(a),ntt(b))) == a*b (mod X^256+1,q).
 * Coefficient ranges and aliasing are documented with the poly functions. */
int ntt_plan_init_256(ntt_plan *out, const ring_ctx *ring, int layers);

#endif
