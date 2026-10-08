#ifndef PQC_POLY_SAMPLE_H
#define PQC_POLY_SAMPLE_H

#include "poly.h"
#include "rand.h"

/* Compatibility sampling utilities, outside the deterministic algebra layer.
 * The injected variants return 0/-1 and leave output unchanged on RNG error.
 * NULL rng requests the OS source. Legacy void wrappers abort on failure. */
int poly_rand_with_rng(lssi_t *a, int coefs, pqc_rng_func rng, void *rng_ctx);
void poly_rand(lssi_t *a, int coefs);
int poly_uniform_CBD_with_rng(poly *a, lp_ctx *ctx,
                              pqc_rng_func rng, void *rng_ctx);
void poly_uniform_CBD(poly *a, lp_ctx *ctx);

/* Legacy eta=2 CBD decoder. Consumes at most buflen bytes, emits at most len
 * coefficients and returns the number emitted. ctx is unused. */
unsigned int sample_CBD(int32_t *a, unsigned int len, const uint8_t *buf,
                        unsigned int buflen, lp_ctx *ctx);

#endif
