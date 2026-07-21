#ifndef KYBER_H
#define KYBER_H

#include <stdint.h>
#include "poly.h"

#define KYBER_N 256
#define KYBER_Q 3329
#define KYBER_QBITS 12
#define KYBER_SYMBYTES 32
#define KYBER_MAX_K 4
#define KYBER_POLYBYTES 384

#define KYBER_MODE_2 2
#define KYBER_MODE_3 3
#define KYBER_MODE_5 5

#define KYBER_MAX_EK_BYTES 1568
#define KYBER_MAX_DK_PKE_BYTES 1536
#define KYBER_MAX_CT_BYTES 1568

typedef struct {
    int mode;
    int k;
    int eta1;
    int eta2;
    int du;
    int dv;
    int ek_bytes;
    int dk_pke_bytes;
    int ct_bytes;
} kyber_param;

void kyber_init_param(kyber_param *param, int mode);
const char *kyber_mode_name(int mode);

int kyber_pke_keygen(uint8_t *ek_pke, uint8_t *dk_pke,
                     const uint8_t d[KYBER_SYMBYTES], kyber_param param);
int kyber_pke_encrypt(uint8_t *ct, const uint8_t *ek_pke,
                      const uint8_t m[KYBER_SYMBYTES],
                      const uint8_t r[KYBER_SYMBYTES], kyber_param param);
int kyber_pke_decrypt(uint8_t m[KYBER_SYMBYTES], const uint8_t *dk_pke,
                      const uint8_t *ct, kyber_param param);

#endif
