#ifndef KYBER_SAMPLE_H
#define KYBER_SAMPLE_H

#include <stddef.h>
#include <stdint.h>
#include "kyber.h"

void kyber_g(uint8_t out[64], const uint8_t *in, size_t inlen);
void kyber_prf(uint8_t *out, size_t outlen, const uint8_t key[KYBER_SYMBYTES],
               uint8_t nonce);

void kyber_sample_ntt(poly a, const uint8_t seed[34], lp_ctx *ctx);
void kyber_sample_cbd(poly a, const uint8_t *buf, int eta);
void kyber_sample_noise(poly a, const uint8_t sigma[KYBER_SYMBYTES],
                        uint8_t nonce, int eta);

#endif
