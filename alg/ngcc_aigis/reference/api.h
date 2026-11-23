#ifndef PQC_AIGIS_REFERENCE_API_H
#define PQC_AIGIS_REFERENCE_API_H

#include <stdint.h>
#include "params.h"

#define KEM_SECRETKEYBYTES SK_BYTES
#define KEM_PUBLICKEYBYTES PK_BYTES
#define KEM_BYTES SEED_BYTES
#define KEM_CIPHERTEXTBYTES CT_BYTES

#define CRYPTO_SECRETKEYBYTES SK_BYTES
#define CRYPTO_PUBLICKEYBYTES PK_BYTES
#define CRYPTO_BYTES SEED_BYTES
#define CRYPTO_CIPHERTEXTBYTES CT_BYTES

int mkem_keygen(uint8_t *pk, uint8_t *sk);
int mkem_enc(uint8_t *pk, uint8_t *ss, uint8_t *ct);
int mdkem_enc(uint8_t *pk, uint8_t *rnd, uint8_t *ss, uint8_t *ct);
int mkem_dec(uint8_t *sk, uint8_t *ct, uint8_t *ss);

#endif
