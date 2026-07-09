#ifndef DILI_SAMPLE_H
#define DILI_SAMPLE_H

#include "poly.h"

void poly_uniform(poly *a, lp_ctx *ctx);
void poly_uniform_eta(poly *a, lp_ctx *ctx, int ETA);
unsigned int sample_uniform(int32_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen, lp_ctx *ctx);
unsigned int sample_eta(int32_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen, lp_ctx *ctx, int ETA);
unsigned int sample_gamma1(int32_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen, lp_ctx *ctx);
unsigned int sample_gamma2(int32_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen, lp_ctx *ctx);

#endif // DILI_SAMPLE_H
