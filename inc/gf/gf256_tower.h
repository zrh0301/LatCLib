#ifndef PQC_GF256_TOWER_H
#define PQC_GF256_TOWER_H
#include <stdint.h>
#include <stddef.h>
/* GF16[y]/(y*y+y+8), byte a0|(a1<<4). Distinct from 0x11b bytes.
 * encode/decode convert external 0x11b bytes to/from tower coordinates. */
uint8_t gf256_tower_encode(uint8_t polynomial_byte);
uint8_t gf256_tower_decode(uint8_t tower_byte);
uint8_t gf256_tower_product(uint8_t a,uint8_t b);
/* Tower coordinates in/out. n counts bytes; exact in-place allowed. */
void gf256_tower_vector_mul(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t n);
#endif
