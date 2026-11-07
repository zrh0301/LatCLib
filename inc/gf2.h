#ifndef PQC_GF2_H
#define PQC_GF2_H
#include <stddef.h>
#include <stdint.h>
/* n counts 64-bit words. n==0 permits NULL. */
void gf2_vector_xor(uint64_t *out,const uint64_t *a,const uint64_t *b,size_t n);
uint8_t gf2_vector_dot(const uint64_t *a,const uint64_t *b,size_t n);
#endif
