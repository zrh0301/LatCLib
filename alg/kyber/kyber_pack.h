#ifndef KYBER_PACK_H
#define KYBER_PACK_H

#include <stdint.h>
#include "kyber.h"

int kyber_param_is_valid(kyber_param param);

lssi_t kyber_reduce(lssi_t a);

lssi_t kyber_compress_one(lssi_t x, int d);
lssi_t kyber_decompress_one(lssi_t y, int d);
void kyber_poly_compress(poly r, const poly a, int d);
void kyber_poly_decompress(poly r, const poly a, int d);

void kyber_byte_encode(uint8_t *out, const poly a, int d);
void kyber_byte_decode(poly a, const uint8_t *in, int d);

void kyber_pack_ek(uint8_t *ek, const poly *t_hat, const uint8_t rho[KYBER_SYMBYTES],
                   kyber_param param);
void kyber_unpack_ek(poly *t_hat, uint8_t rho[KYBER_SYMBYTES], const uint8_t *ek,
                     kyber_param param);

void kyber_pack_dk_pke(uint8_t *dk_pke, const poly *s_hat, kyber_param param);
void kyber_unpack_dk_pke(poly *s_hat, const uint8_t *dk_pke, kyber_param param);

void kyber_pack_ct(uint8_t *ct, const poly *u, const poly v, kyber_param param);
void kyber_unpack_ct(poly *u, poly v, const uint8_t *ct, kyber_param param);

int kyber_check_ek_modulus(const uint8_t *ek, kyber_param param);

#endif
