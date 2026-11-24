#ifndef PQC_BITS_H
#define PQC_BITS_H

#include <stddef.h>
#include <stdint.h>

/* LSB-first bit packing shared by the multivariate candidates.
 *
 * Bit i of the value is stored in bit (offset + i) of the buffer, so the low
 * bits of the value land in the low bits of the first byte. Both helpers take
 * an explicit bit offset instead of moving a cursor, which lets callers keep
 * their own packing state. Up to 32-bit values are supported.
 *
 * pqc_store_bits_le ORs its bits into the buffer: the destination must already
 * be zeroed (or contain only bits written by earlier calls).
 * pqc_load_bits_le reads only the bytes that the requested bits touch. */

static inline uint32_t pqc_bit_mask(unsigned int nbits)
{
    return nbits >= 32 ? UINT32_MAX : (uint32_t)((UINT32_C(1) << nbits) - 1);
}

static inline void pqc_store_bits_le(uint8_t *dst, size_t bit_offset,
                                     uint32_t value, unsigned int nbits)
{
    const unsigned int shift = (unsigned int)(bit_offset & 7u);
    const size_t index = bit_offset >> 3;
    const uint64_t x = (uint64_t)(value & pqc_bit_mask(nbits)) << shift;

    dst[index] |= (uint8_t)x;
    if (shift + nbits > 8)
        dst[index + 1] |= (uint8_t)(x >> 8);
    if (shift + nbits > 16)
        dst[index + 2] |= (uint8_t)(x >> 16);
    if (shift + nbits > 24)
        dst[index + 3] |= (uint8_t)(x >> 24);
    if (shift + nbits > 32)
        dst[index + 4] |= (uint8_t)(x >> 32);
}

static inline uint32_t pqc_load_bits_le(const uint8_t *src, size_t bit_offset,
                                        unsigned int nbits)
{
    const unsigned int shift = (unsigned int)(bit_offset & 7u);
    const size_t index = bit_offset >> 3;
    uint64_t x = src[index];

    if (shift + nbits > 8)
        x |= (uint64_t)src[index + 1] << 8;
    if (shift + nbits > 16)
        x |= (uint64_t)src[index + 2] << 16;
    if (shift + nbits > 24)
        x |= (uint64_t)src[index + 3] << 24;
    if (shift + nbits > 32)
        x |= (uint64_t)src[index + 4] << 32;
    return (uint32_t)((x >> shift) & pqc_bit_mask(nbits));
}

#endif
