#include "poly/ntt16.h"

#include "poly/fp16.h"

#define NTT16_DIM 128

/* Official Kyber twiddles for q = 3329, centered Montgomery form, index 0
 * unused by the forward transform. The inverse table is the one from the
 * reference samplers and already folds in the 1/128 scaling. */
static const int16_t zetas_3329[128] = {
    -1044, -758, -359, -1517, 1493, 1422, 287, 202,
    -171, 622, 1577, 182, 962, -1202, -1474, 1468,
    573, -1325, 264, 383, -829, 1458, -1602, -130,
    -681, 1017, 732, 608, -1542, 411, -205, -1571,
    1223, 652, -552, 1015, -1293, 1491, -282, -1544,
    516, -8, -320, -666, -1618, -1162, 126, 1469,
    -853, -90, -271, 830, 107, -1421, -247, -951,
    -398, 961, -1508, -725, 448, -1065, 677, -1275,
    -1103, 430, 555, 843, -1251, 871, 1550, 105,
    422, 587, 177, -235, -291, -460, 1574, 1653,
    -246, 778, 1159, -147, -777, 1483, -602, 1119,
    -1590, 644, -872, 349, 418, 329, -156, -75,
    817, 1097, 603, 610, 1322, -1285, -1465, 384,
    -1215, -136, 1218, -1335, -874, 220, -1187, -1659,
    -1185, -1530, -1278, 794, -1510, -854, -870, 478,
    -108, -308, 996, 991, 958, -1460, 1522, 1628
};

static const int16_t zetas_inv_3329[128] = {
    -1628, -1522, 1460, -958, -991, -996, 308, 108,
    -478, 870, 854, 1510, -794, 1278, 1530, 1185,
    1659, 1187, -220, 874, 1335, -1218, 136, 1215,
    -384, 1465, 1285, -1322, -610, -603, -1097, -817,
    75, 156, -329, -418, -349, 872, -644, 1590,
    -1119, 602, -1483, 777, 147, -1159, -778, 246,
    -1653, -1574, 460, 291, 235, -177, -587, -422,
    -105, -1550, -871, 1251, -843, -555, -430, 1103,
    1275, -677, 1065, -448, 725, 1508, -961, 398,
    951, 247, 1421, -107, -830, 271, 90, 853,
    -1469, -126, 1162, 1618, 666, 320, 8, -516,
    1544, 282, -1491, 1293, -1015, 552, -652, -1223,
    1571, 205, -411, 1542, -608, -732, -1017, 681,
    130, 1602, -1458, 829, -383, -264, 1325, -573,
    -1468, 1474, 1202, -962, -182, -1577, -622, 171,
    -202, -287, -1422, -1493, 1517, 359, 266
};

static const fp16_ctx ntt16_field = {
    3329, 62209, 9, (uint16_t)(((65536u % 3329) * (65536u % 3329)) % 3329)
};

void pqc_ntt128(int16_t *a)
{
    unsigned k = 1;

    for (int level = 6; level >= 1; --level) {
        const int step = 1 << level;
        for (int start = 0; start < NTT16_DIM; start += step << 1) {
            const int32_t zeta = zetas_3329[k++];
            for (int j = start; j < start + step; ++j) {
                const int16_t t = fp16_montgomery_reduce((int32_t)zeta * a[j + step],
                                                        &ntt16_field);
                a[j + step] = (int16_t)(a[j] - t);
                a[j] = (int16_t)(a[j] + t);
            }
        }
    }
    /* Final layer also reduces both outputs so the result is in (-q, q). */
    for (int start = 0; start < NTT16_DIM; start += 2) {
        const int32_t zeta = zetas_3329[k++];
        const int16_t t = fp16_montgomery_reduce((int32_t)zeta * a[start + 1],
                                                &ntt16_field);
        a[start + 1] = fp16_barrett_reduce((int16_t)(a[start] - t), &ntt16_field);
        a[start] = fp16_barrett_reduce((int16_t)(a[start] + t), &ntt16_field);
    }
}

void pqc_invntt128(int16_t *a)
{
    unsigned k = 0;

    /* Layers 0 and 1 accumulate lazily; layers 1, 3/4 and 5 reduce, which keeps
     * every partial sum inside int16_t for q = 3329. */
    for (int level = 0; level < 2; ++level) {
        const int step = 1 << level;
        for (int start = 0; start < NTT16_DIM; start += step << 1) {
            const int32_t zeta = zetas_inv_3329[k++];
            for (int j = start; j < start + step; ++j) {
                const int16_t t = a[j];
                a[j] = (int16_t)(t + a[j + step]);
                a[j + step] = fp16_montgomery_reduce((int32_t)zeta *
                                                     (int16_t)(t - a[j + step]),
                                                     &ntt16_field);
            }
        }
    }

    for (int start = 0; start < NTT16_DIM; start += 8) {
        const int32_t zeta = zetas_inv_3329[k++];
        for (int j = start; j < start + 4; ++j) {
            const int16_t t = a[j];
            a[j] = fp16_barrett_reduce((int16_t)(t + a[j + 4]), &ntt16_field);
            a[j + 4] = fp16_montgomery_reduce((int32_t)zeta *
                                              (int16_t)(t - a[j + 4]),
                                              &ntt16_field);
        }
    }

    for (int level = 3; level < 5; ++level) {
        const int step = 1 << level;
        for (int start = 0; start < NTT16_DIM; start += step << 1) {
            const int32_t zeta = zetas_inv_3329[k++];
            for (int j = start; j < start + step; ++j) {
                const int16_t t = a[j];
                a[j] = (int16_t)(t + a[j + step]);
                a[j + step] = fp16_montgomery_reduce((int32_t)zeta *
                                                     (int16_t)(t - a[j + step]),
                                                     &ntt16_field);
            }
        }
    }

    for (int start = 0; start < NTT16_DIM; start += 64) {
        const int32_t zeta = zetas_inv_3329[k++];
        for (int j = start; j < start + 32; ++j) {
            const int16_t t = a[j];
            a[j] = fp16_barrett_reduce((int16_t)(t + a[j + 32]), &ntt16_field);
            a[j + 32] = fp16_montgomery_reduce((int32_t)zeta *
                                               (int16_t)(t - a[j + 32]),
                                               &ntt16_field);
        }
    }

    {
        const int32_t zeta = zetas_inv_3329[k++];
        for (int j = 0; j < 64; ++j) {
            const int16_t t = a[j];
            /* 512 = R/128, so this step applies the 1/n scaling. */
            a[j] = fp16_montgomery_reduce((int32_t)512 * (int16_t)(t + a[j + 64]),
                                          &ntt16_field);
            a[j + 64] = fp16_montgomery_reduce((int32_t)zeta *
                                               (int16_t)(t - a[j + 64]),
                                               &ntt16_field);
        }
    }
}

void pqc_ntt128_blocks(int16_t *a, size_t blocks)
{
    for (size_t i = 0; i < blocks; ++i) {
        pqc_ntt128(a + i * NTT16_DIM);
    }
}

void pqc_invntt128_blocks(int16_t *a, size_t blocks)
{
    for (size_t i = 0; i < blocks; ++i) {
        pqc_invntt128(a + i * NTT16_DIM);
    }
}

/* A 128-point transform evaluates y = x^blocks at the 128 roots of
 * y^128 + 1. The final butterfly places opposite roots in adjacent lanes. */
static inline int16_t ntt16_lane_twiddle(size_t lane)
{
    const int16_t zeta = zetas_3329[64 + lane / 2];
    return lane % 2 ? (int16_t)-zeta : zeta;
}

static inline int16_t ntt16_mont(int32_t value)
{
    return fp16_montgomery_reduce(value, &ntt16_field);
}

static inline int16_t ntt16_barrett(int16_t value)
{
    return fp16_barrett_reduce(value, &ntt16_field);
}

/* Karatsuba convolution in each lane, combining chunks in block order. */
static inline void mont_mul2(int32_t out[3], int32_t a0, int32_t a1,
                             int32_t b0, int32_t b1)
{
    out[0] = ntt16_mont(a0 * b0);
    out[2] = ntt16_mont(a1 * b1);
    out[1] = ntt16_mont((a0 + a1) * (b0 + b1)) - out[0] - out[2];
}

static inline void mont_mul4_lane(int32_t out[7], const int16_t *a, const int16_t *b)
{
    int32_t low[3], high[3], cross[3];
    const int32_t a0 = a[0], a1 = a[NTT16_DIM];
    const int32_t a2 = a[2 * NTT16_DIM], a3 = a[3 * NTT16_DIM];
    const int32_t b0 = b[0], b1 = b[NTT16_DIM];
    const int32_t b2 = b[2 * NTT16_DIM], b3 = b[3 * NTT16_DIM];

    mont_mul2(low, a0, a1, b0, b1);
    mont_mul2(high, a2, a3, b2, b3);
    mont_mul2(cross, a0 + a2, a1 + a3, b0 + b2, b1 + b3);

    out[0] = low[0];
    out[1] = low[1];
    out[2] = low[2] + cross[0] - low[0] - high[0];
    out[3] = cross[1] - low[1] - high[1];
    out[4] = high[0] + cross[2] - low[2] - high[2];
    out[5] = high[1];
    out[6] = high[2];
}

static void mont_mul4(int16_t *out, const int16_t *a, const int16_t *b)
{
    for (size_t i = 0; i < NTT16_DIM; ++i) {
        int32_t terms[7];
        mont_mul4_lane(terms, a + i, b + i);
        out[i] = (int16_t)terms[0];
        for (size_t j = 1; j < 6; ++j)
            out[j * NTT16_DIM + i] = ntt16_barrett((int16_t)terms[j]);
        out[6 * NTT16_DIM + i] = (int16_t)terms[6];
    }
}

static inline void mont_mul8(int16_t *r, const int16_t *a, const int16_t *b)
{
    int16_t t1[7 * NTT16_DIM], t2[7 * NTT16_DIM];
    int i;
    for (i = 0; i < 4 * NTT16_DIM; i++)
    {
        t1[i] = a[i] + a[4 * NTT16_DIM + i];
        t2[i] = b[i] + b[4 * NTT16_DIM + i];
    }
    mont_mul4(t1, t1, t2);
    mont_mul4(&r[8 * NTT16_DIM], &a[4 * NTT16_DIM], &b[4 * NTT16_DIM]);
    mont_mul4(r, a, b);

    for (i = 0; i < NTT16_DIM; i++)
        r[7 * NTT16_DIM + i] = ntt16_barrett(t1[3 * NTT16_DIM + i] - r[3 * NTT16_DIM + i] - r[11 * NTT16_DIM + i]);
    for (i = 0; i < 3 * NTT16_DIM; i++)
    {
        t1[i] = ntt16_barrett(r[4 * NTT16_DIM + i] + t1[i] - r[i] - r[8 * NTT16_DIM + i]);
        r[8 * NTT16_DIM + i] = ntt16_barrett(t1[4 * NTT16_DIM + i] - r[4 * NTT16_DIM + i] - r[12 * NTT16_DIM + i] + r[8 * NTT16_DIM + i]);
        r[4 * NTT16_DIM + i] = t1[i];
    }
}

static void mont_product4(int16_t *out, const int16_t *a, const int16_t *b)
{
    for (size_t i = 0; i < NTT16_DIM; ++i) {
        int32_t terms[7];
        const int32_t root = ntt16_lane_twiddle(i);
        mont_mul4_lane(terms, a + i, b + i);
        for (size_t j = 0; j < 3; ++j)
            out[j * NTT16_DIM + i] = (int16_t)(terms[j] + ntt16_mont(terms[j + 4] * root));
        out[3 * NTT16_DIM + i] = (int16_t)terms[3];
    }
}

static void mont_product8(int16_t *r, const int16_t *a, const int16_t *b)
{
    enum { N = 8 * NTT16_DIM };
    int16_t t0[N], t1[N], t2[N];
    int i, j;
    for (i = 0; i < N / 2; i++)
    {
        t0[i] = a[i] + a[N / 2 + i];
        t1[i] = b[i] + b[N / 2 + i];
    }

    mont_mul4(t1, t0, t1);
    mont_mul4(t0, a, b);
    mont_mul4(t2, &a[N / 2], &b[N / 2]);

    for (i = 0; i < 7 * NTT16_DIM; i++)
        t1[i] -= (t0[i] + t2[i]);

    for (i = 0; i < 3 * NTT16_DIM; i++)
    {
        t0[4 * NTT16_DIM + i] = (t0[4 * NTT16_DIM + i] + t1[i]);
        t2[i] = (t1[4 * NTT16_DIM + i] + t2[i]);
    }

    for (i = 0; i < 7; i++)
        for (j = 0; j < NTT16_DIM; j++)
            t2[i * NTT16_DIM + j] = ntt16_mont(t2[i * NTT16_DIM + j] * ntt16_lane_twiddle(j));

    for (i = 0; i < 7 * NTT16_DIM; i++)
        r[i] = t0[i] + t2[i];
    for (i = 0; i < NTT16_DIM; i++)
        r[7 * NTT16_DIM + i] = t1[3 * NTT16_DIM + i];
}

static void mont_product16(int16_t *r, const int16_t *a, const int16_t *b)
{
    enum { N = 16 * NTT16_DIM };
    int16_t t0[N], t1[N], t2[N];
    int i, j;
    for (i = 0; i < N / 2; i++)
    {
        t0[i] = a[i] + a[N / 2 + i];
        t1[i] = b[i] + b[N / 2 + i];
    }

    mont_mul8(t1, t0, t1);
    mont_mul8(t0, a, b);
    mont_mul8(t2, &a[N / 2], &b[N / 2]);

    for (i = 0; i < 15 * NTT16_DIM; i++)
        t1[i] -= (t0[i] + t2[i]);

    for (i = 0; i < 7 * NTT16_DIM; i++)
    {
        t0[8 * NTT16_DIM + i] = (t0[8 * NTT16_DIM + i] + t1[i]);
        t2[i] = (t1[8 * NTT16_DIM + i] + t2[i]);
    }

    for (i = 0; i < 15; i++)
        for (j = 0; j < NTT16_DIM; j++)
            t2[i * NTT16_DIM + j] = ntt16_mont(t2[i * NTT16_DIM + j] * ntt16_lane_twiddle(j));

    for (i = 0; i < 15 * NTT16_DIM; i++)
        r[i] = t0[i] + t2[i];
    for (i = 0; i < NTT16_DIM; i++)
        r[15 * NTT16_DIM + i] = t1[7 * NTT16_DIM + i];
}

int pqc_ntt128_blocks_mont_mul(int16_t *out, const int16_t *a, const int16_t *b,
                                size_t blocks)
{
    switch (blocks) {
    case 4:
        mont_product4(out, a, b);
        return 0;
    case 8:
        mont_product8(out, a, b);
        return 0;
    case 16:
        mont_product16(out, a, b);
        return 0;
    default:
        return -1;
    }
}
