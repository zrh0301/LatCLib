#ifndef PQC_GF256_H
#define PQC_GF256_H
#include <stddef.h>
#include <stdint.h>
/* Polynomial basis modulo 0x11b. inv(0)=0. */
uint8_t gf256_mul(uint8_t a,uint8_t b);
uint8_t gf256_inv(uint8_t a);
/* Exact in-place allowed; no partial overlap. n==0 permits NULL. */
void gf256_vector_mul(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t n);
void gf256_vector_madd(uint8_t *acc,const uint8_t *a,uint8_t k,size_t n);
#endif
