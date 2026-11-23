#ifndef PQC_GF16_VECTOR_H
#define PQC_GF16_VECTOR_H
#include <stddef.h>
#include <stdint.h>
/* Raw: one GF16 element in the low nibble per byte, modulus 0x13.
 * n counts bytes. Exact in-place output is allowed; partial overlap is not.
 * n==0 permits NULL; otherwise pointers must name n accessible bytes. */
void gf16_vector_mul(uint8_t *out,const uint8_t *a,const uint8_t *b,size_t n);
void gf16_vector_madd(uint8_t *acc,const uint8_t *a,uint8_t k,size_t n);
/* Two independent GF16 elements per byte, low nibble first. n is BYTE count. */
void gf16_packed_mul_const(uint8_t *out,const uint8_t *a,uint8_t k,size_t n);

/* Nibble access is independent of the basis used for field multiplication. */
static inline uint8_t gf16_packed_get(const uint8_t *in, size_t index)
{
    return (uint8_t)((in[index / 2] >> (4 * (index % 2))) & 15u);
}

static inline void gf16_packed_set(uint8_t *out, size_t index, uint8_t value)
{
    const unsigned shift = 4 * (index % 2);
    out[index / 2] = (uint8_t)((out[index / 2] & ~(15u << shift)) |
                               ((value & 15u) << shift));
}
#endif
