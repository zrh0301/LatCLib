#include "gf16.h"

gf16_t gf16_add(gf16_t a, gf16_t b) {
    return (gf16_t)((a ^ b) & 0x0fu);
}

gf16_t gf16_sub(gf16_t a, gf16_t b) {
    return gf16_add(a, b);
}

gf16_t gf16_mul(gf16_t a, gf16_t b) {
    gf16_t product = 0;
    a &= 0x0fu;
    b &= 0x0fu;

    for (unsigned int i = 0; i < 4; ++i) {
        gf16_t add_mask = (gf16_t)(0u - (b & 1u));
        gf16_t reduce_mask = (gf16_t)(0u - ((a >> 3) & 1u));
        product ^= (gf16_t)(a & add_mask);
        a = (gf16_t)(((a << 1) & 0x0fu) ^ (0x03u & reduce_mask));
        b >>= 1;
    }
    return product;
}

gf16_t gf16_inv(gf16_t a) {
    gf16_t a2 = gf16_mul(a, a);
    gf16_t a4 = gf16_mul(a2, a2);
    gf16_t a8 = gf16_mul(a4, a4);
    gf16_t a12 = gf16_mul(a8, a4);
    return gf16_mul(a12, a2);
}
