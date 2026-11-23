#ifndef PQC_AIGIS_OWCPA_H
#define PQC_AIGIS_OWCPA_H

#include <stdint.h>
#include "params.h"

int pqc_aigis_owcpa_keypair(uint8_t *pk, uint8_t *sk, const uint8_t *seeds);
int pqc_aigis_owcpa_enc(uint8_t *ct, const uint8_t *m, const uint8_t *pk,
                       const uint8_t *coins);
void pqc_aigis_owcpa_dec(uint8_t *m, const uint8_t *ct, const uint8_t *sk);

#endif
