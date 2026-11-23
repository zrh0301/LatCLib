/* Arithmetic for the GF(4)-based tower coordinates used by VDOO.
 * The bit-sliced word kernels follow its reference implementation. */
#include "gf/composite.h"

#include <string.h>

static inline uint8_t gf4_mul_2(uint8_t a) {
    uint8_t r = a << 1;
    r ^= (a >> 1) * 7;
    return r;
}

static inline uint8_t gf4_mul_3(uint8_t a) {
    uint8_t msk = (a - 2) >> 1;
    return (msk & (a * 3)) | ((~msk) & (a - 1));
}

static inline uint8_t gf4_mul(uint8_t a, uint8_t b) {
    uint8_t r = a * (b & 1);
    return r ^ (gf4_mul_2(a) * (b >> 1));
}

static inline uint8_t gf4_squ(uint8_t a) {
    return a ^ (a >> 1);
}

static inline uint32_t gf4v_mul_2_u32(uint32_t a) {
    uint32_t bit0 = a & 0x55555555;
    uint32_t bit1 = a & 0xaaaaaaaa;
    return (bit0 << 1) ^ bit1 ^ (bit1 >> 1);
}

static inline uint32_t gf4v_mul_3_u32(uint32_t a) {
    uint32_t bit0 = a & 0x55555555;
    uint32_t bit1 = a & 0xaaaaaaaa;
    return (bit0 << 1) ^ bit0 ^ (bit1 >> 1);
}

static inline uint32_t gf4v_mul_u32(uint32_t a, uint8_t b) {
    uint32_t bit0_b = ((uint32_t) 0) - ((uint32_t)(b & 1));
    uint32_t bit1_b = ((uint32_t) 0) - ((uint32_t)((b >> 1) & 1));
    return (a & bit0_b) ^ (bit1_b & gf4v_mul_2_u32(a));
}

uint8_t gf16_composite_mul(uint8_t a, uint8_t b) {
    uint8_t a0 = a & 3;
    uint8_t a1 = (a >> 2);
    uint8_t b0 = b & 3;
    uint8_t b1 = (b >> 2);
    uint8_t a0b0 = gf4_mul(a0, b0);
    uint8_t a1b1 = gf4_mul(a1, b1);
    uint8_t a0b1_a1b0 = gf4_mul(a0 ^ a1, b0 ^ b1) ^ a0b0 ^ a1b1;
    uint8_t a1b1_x2 = gf4_mul_2(a1b1);
    return ((a0b1_a1b0 ^ a1b1) << 2) ^ a0b0 ^ a1b1_x2;
}

static inline uint8_t gf16_squ(uint8_t a) {
    uint8_t a0 = a & 3;
    uint8_t a1 = (a >> 2);
    a1 = gf4_squ(a1);
    uint8_t a1squ_x2 = gf4_mul_2(a1);
    return (a1 << 2) ^ a1squ_x2 ^ gf4_squ(a0);
}

uint8_t gf16_composite_inv(uint8_t a) {
    uint8_t a2 = gf16_squ(a);
    uint8_t a4 = gf16_squ(a2);
    uint8_t a8 = gf16_squ(a4);
    uint8_t a6 = gf16_composite_mul(a4, a2);
    return gf16_composite_mul(a8, a6);
}

static inline uint8_t gf16_mul_8(uint8_t a) {
    uint8_t a0 = a & 3;
    uint8_t a1 = a >> 2;
    return (gf4_mul_2(a0 ^ a1) << 2) | gf4_mul_3(a1);
}

static inline uint32_t gf16v_mul_u32(uint32_t a, uint8_t b) {
    uint32_t axb0 = gf4v_mul_u32(a, b);
    uint32_t axb1 = gf4v_mul_u32(a, b >> 2);
    uint32_t a0b1 = (axb1 << 2) & 0xcccccccc;
    uint32_t a1b1 = axb1 & 0xcccccccc;
    uint32_t a1b1_2 = a1b1 >> 2;

    return axb0 ^ a0b1 ^ a1b1 ^ gf4v_mul_2_u32(a1b1_2);
}

static inline uint32_t gf16v_mul_8_u32(uint32_t a) {
    uint32_t a1 = a & 0xcccccccc;
    uint32_t a0 = (a << 2) & 0xcccccccc;
    return gf4v_mul_2_u32(a0 ^ a1) | gf4v_mul_3_u32(a1 >> 2);
}

uint8_t gf256_composite_mul(uint8_t a, uint8_t b) {
    uint8_t a0 = a & 15;
    uint8_t a1 = (a >> 4);
    uint8_t b0 = b & 15;
    uint8_t b1 = (b >> 4);
    uint8_t a0b0 = gf16_composite_mul(a0, b0);
    uint8_t a1b1 = gf16_composite_mul(a1, b1);
    uint8_t a0b1_a1b0 = gf16_composite_mul(a0 ^ a1, b0 ^ b1) ^ a0b0 ^ a1b1;
    uint8_t a1b1_x8 = gf16_mul_8(a1b1);
    return ((a0b1_a1b0 ^ a1b1) << 4) ^ a0b0 ^ a1b1_x8;
}

static inline uint8_t gf256_squ(uint8_t a) {
    uint8_t a0 = a & 15;
    uint8_t a1 = (a >> 4);
    a1 = gf16_squ(a1);
    uint8_t a1squ_x8 = gf16_mul_8(a1);
    return (a1 << 4) ^ a1squ_x8 ^ gf16_squ(a0);
}

uint8_t gf256_composite_inv(uint8_t a) {
    // 128+64+32+16+8+4+2 = 254
    uint8_t a2 = gf256_squ(a);
    uint8_t a4 = gf256_squ(a2);
    uint8_t a8 = gf256_squ(a4);
    uint8_t a4_2 = gf256_composite_mul(a4, a2);
    uint8_t a8_4_2 = gf256_composite_mul(a4_2, a8);
    uint8_t a64_ = gf256_squ(a8_4_2);
    a64_ = gf256_squ(a64_);
    a64_ = gf256_squ(a64_);
    uint8_t a64_2 = gf256_composite_mul(a64_, a8_4_2);
    uint8_t a128_ = gf256_squ(a64_2);
    return gf256_composite_mul(a2, a128_);
}

static inline uint32_t gf256v_mul_u32(uint32_t a, uint8_t b) {
    uint32_t axb0 = gf16v_mul_u32(a, b);
    uint32_t axb1 = gf16v_mul_u32(a, b >> 4);
    uint32_t a0b1 = (axb1 << 4) & 0xf0f0f0f0;
    uint32_t a1b1 = axb1 & 0xf0f0f0f0;
    uint32_t a1b1_4 = a1b1 >> 4;

    return axb0 ^ a0b1 ^ a1b1 ^ gf16v_mul_8_u32(a1b1_4);
}

static void vector_madd(uint8_t *acc, const uint8_t *in, uint8_t scalar,
                        size_t len, int packed)
{
    while (len >= sizeof(uint32_t)) {
        uint32_t x, y;
        memcpy(&x, in, sizeof(x));
        memcpy(&y, acc, sizeof(y));
        y ^= packed ? gf16v_mul_u32(x, scalar) : gf256v_mul_u32(x, scalar);
        memcpy(acc, &y, sizeof(y));
        in += sizeof(x);
        acc += sizeof(y);
        len -= sizeof(x);
    }
    if (len != 0) {
        uint32_t x = 0, y = 0;
        memcpy(&x, in, len);
        memcpy(&y, acc, len);
        y ^= packed ? gf16v_mul_u32(x, scalar) : gf256v_mul_u32(x, scalar);
        memcpy(acc, &y, len);
    }
}

void gf16_composite_packed_madd(uint8_t *acc, const uint8_t *in, uint8_t scalar,
                                size_t len)
{
    vector_madd(acc, in, scalar, len, 1);
}

void gf256_composite_madd(uint8_t *acc, const uint8_t *in, uint8_t scalar,
                          size_t len)
{
    vector_madd(acc, in, scalar, len, 0);
}
