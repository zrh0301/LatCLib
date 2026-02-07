#ifndef PACK_H
#define PACK_H

#include <stdint.h>
#include "poly.h"

void dili_polyt1_pack(uint8_t *r, const poly p, int n, int d);
void dili_polyt1_unpack(poly p, const uint8_t *a, int n, int d);

void dili_pack_pk(uint8_t *pk_bytes, const uint8_t *rho, const poly *t1, int k);
void dili_unpack_pk(uint8_t *rho, poly *t1, const uint8_t *pk_bytes, int k);

void pack_z(uint8_t* r, const poly z);
void unpack_z(poly z, const uint8_t* r);

void pack_w1(uint8_t* buf, const poly w1);

int pack_h(uint8_t* r, uint8_t h[][256], int k);
void unpack_h(poly* h, const uint8_t* r, int k);

#endif // PACK_H
