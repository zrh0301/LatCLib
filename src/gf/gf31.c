#include "gf/gf31.h"

/* Multiplicative inverses in GF(31); entry 0 is unused. */
const gf31_t pqc_gf31_inv_table[31] = {
    0,  1, 16, 21,  8, 25, 26,  9,
    4,  7, 28, 17, 13, 12, 20, 29,
    2, 11, 19, 18, 14,  3, 24, 27,
    22, 5,  6, 23, 10, 15, 30
};

gf31_t gf31_reduce(int z)
{
    /* 2^5 = 1 (mod 31): fold bits 5 and up onto the low five bits twice,
     * then subtract 31 once if the result reached or passed the modulus. */
    z = (z & 31) + ((z & ~31) >> 5);
    const int carry = (z + 1) & ~31;
    z += carry >> 5;
    z -= carry;
    return (gf31_t)z;
}

gf31_t gf31_add(gf31_t a, gf31_t b)
{
    return gf31_reduce((int)a + (int)b);
}

gf31_t gf31_sub(gf31_t a, gf31_t b)
{
    return gf31_reduce((int)a - (int)b + 31);
}

gf31_t gf31_mul(gf31_t a, gf31_t b)
{
    return gf31_reduce((int)a * (int)b);
}

gf31_t gf31_inv(gf31_t a)
{
    return pqc_gf31_inv_table[a];
}
