#include "gf/gf31.h"
#include "gf/gf31_matrix.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "gf31 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

/* Independent reference for GF(31) residue arithmetic. */
static unsigned ref_mod(int v)
{
    int r = v % 31;
    return (unsigned)(r < 0 ? r + 31 : r);
}

/* GF(31^2) multiply with x^2 = 3, planar pairs (a0, a1). */
static void ref_t2_mul(const uint8_t *a, const uint8_t *b, uint8_t *out)
{
    unsigned a0 = a[0], a1 = a[1], b0 = b[0], b1 = b[1];
    out[0] = (uint8_t)ref_mod((int)(a0 * b0) + 3 * (int)(a1 * b1));
    out[1] = (uint8_t)ref_mod((int)(a0 * b1) + (int)(a1 * b0));
}

int main(void)
{
    /* Field laws against the reference reduction, including the inverse table. */
    for (int a = 0; a < 31; ++a) {
        CHECK(gf31_reduce(a) == ref_mod(a));
        CHECK(gf31_add((gf31_t)a, 0) == (gf31_t)a);
        CHECK(gf31_mul((gf31_t)a, 1) == (gf31_t)a);
        CHECK(gf31_mul((gf31_t)a, 0) == 0);
        CHECK(gf31_sub((gf31_t)a, (gf31_t)a) == 0);
        for (int b = 0; b < 31; ++b) {
            CHECK(gf31_add((gf31_t)a, (gf31_t)b) == ref_mod(a + b));
            CHECK(gf31_sub((gf31_t)a, (gf31_t)b) == ref_mod(a - b));
            CHECK(gf31_mul((gf31_t)a, (gf31_t)b) == ref_mod(a * b));
        }
        if (a != 0)
            CHECK(gf31_mul((gf31_t)a, gf31_inv((gf31_t)a)) == 1);
    }
    CHECK(gf31_inv(0) == 0);

    /* 2^5 = 1 (mod 31) lets reduce fold the documented input range exactly. */
    CHECK(gf31_reduce(31 * 31 + 30) == 30);
    CHECK(gf31_reduce(1022) == 1022 % 31);
    for (int z = 0; z <= 1022; ++z)
        CHECK(gf31_reduce(z) == ref_mod(z));

    /* Vector operators. */
    gf31_t a[64], b[64], out[64], want[64];
    for (int i = 0; i < 64; ++i) {
        a[i] = (gf31_t)((i * 7 + 3) % 31);
        b[i] = (gf31_t)((i * 11 + 5) % 31);
    }

    memcpy(out, a, sizeof(out));
    gf31_vector_add(out, b, 64);
    for (int i = 0; i < 64; ++i)
        CHECK(out[i] == ref_mod((int)a[i] + (int)b[i]));

    {
        unsigned total = 0;
        for (int i = 0; i < 64; ++i)
            total += (unsigned)a[i] * b[i];
        CHECK(gf31_vector_dot(a, b, 64) == (gf31_t)ref_mod((int)total));
    }

    gf31_vector_mul_scalar(a, 9, 64, out);
    for (int i = 0; i < 64; ++i)
        CHECK(out[i] == ref_mod(9 * (int)a[i]));

    gf31_vector_add_no_mod(a, b, 64, out);
    for (int i = 0; i < 64; ++i)
        CHECK(out[i] == (gf31_t)(a[i] + b[i]));
    gf31_vector_reduce(out, 64, want);
    for (int i = 0; i < 64; ++i)
        CHECK(want[i] == ref_mod((int)a[i] + (int)b[i]));

    gf31_vector_add_mod(a, b, 64, out);
    for (int i = 0; i < 64; ++i)
        CHECK(out[i] == ref_mod((int)a[i] + (int)b[i]));

    gf31_vector_sub_mod(a, b, 64, out);
    for (int i = 0; i < 64; ++i)
        CHECK(out[i] == ref_mod((int)a[i] - (int)b[i]));

    /* Planar coordinate map: (a0, a1) -> (a0 + a1, 8 * a1). */
    {
        uint8_t src[8] = {1, 2, 3, 4, 5, 6, 7, 8}, dst[8] = {0}, exp[8];
        gf31_t2_phi_planar(src, dst, 4);
        for (int i = 0; i < 4; ++i) {
            exp[i] = (uint8_t)ref_mod(src[i] + src[i + 4]);
            exp[i + 4] = (uint8_t)ref_mod(8 * src[i + 4]);
        }
        CHECK(memcmp(dst, exp, sizeof(dst)) == 0);
    }

    /* Tower matrix product, 4x4 with the reference layout: A is row-major with
     * padded rows, B and C are column-major (C without padding). */
    {
        enum { ROWS = 4, COLS = 4, PADDED = 32 };
        uint8_t A[ROWS * 2 * PADDED], B[COLS * 2 * PADDED], C[COLS * 2 * ROWS];
        uint8_t Aref[ROWS][COLS][2], Bref[COLS][COLS][2];

        memset(A, 0, sizeof(A));
        memset(B, 0, sizeof(B));
        for (int i = 0; i < ROWS; ++i)
            for (int j = 0; j < COLS; ++j) {
                Aref[i][j][0] = (uint8_t)((i * 5 + j * 3) % 31);
                Aref[i][j][1] = (uint8_t)((i * 2 + j * 7 + 1) % 31);
                A[i * 2 * PADDED + j] = Aref[i][j][0];
                A[i * 2 * PADDED + PADDED + j] = Aref[i][j][1];
                Bref[i][j][0] = (uint8_t)((i * 11 + j + 2) % 31);
                Bref[i][j][1] = (uint8_t)((i * 3 + j * 9 + 4) % 31);
                B[j * 2 * PADDED + i] = Bref[i][j][0];
                B[j * 2 * PADDED + PADDED + i] = Bref[i][j][1];
            }

        gf31_t2_matrix_mul(A, B, C, ROWS, COLS, COLS, ROWS);

        for (int i = 0; i < ROWS; ++i)
            for (int j = 0; j < COLS; ++j) {
                uint8_t prod[2], sum[2] = {0, 0};
                for (int k = 0; k < COLS; ++k) {
                    ref_t2_mul(Aref[i][k], Bref[k][j], prod);
                    sum[0] = (uint8_t)ref_mod(sum[0] + prod[0]);
                    sum[1] = (uint8_t)ref_mod(sum[1] + prod[1]);
                }
                CHECK(C[j * 2 * ROWS + i] == sum[0]);
                CHECK(C[j * 2 * ROWS + ROWS + i] == sum[1]);
            }
    }

    /* Row vector times matrix, planar vector layout with 32-byte padding. */
    {
        enum { VCOLS = 5, MCOLS = 3, PADDED = 32 };
        uint8_t v[2 * PADDED], M[MCOLS * 2 * PADDED], C[2 * MCOLS], vref[VCOLS][2], Mref[VCOLS][MCOLS][2];

        memset(v, 0, sizeof(v));
        memset(M, 0, sizeof(M));
        for (int k = 0; k < VCOLS; ++k) {
            vref[k][0] = (uint8_t)((k * 4 + 1) % 31);
            vref[k][1] = (uint8_t)((k * 6 + 2) % 31);
            v[k] = vref[k][0];
            v[PADDED + k] = vref[k][1];
        }
        for (int j = 0; j < MCOLS; ++j)
            for (int k = 0; k < VCOLS; ++k) {
                Mref[k][j][0] = (uint8_t)((j * 13 + k * 5 + 3) % 31);
                Mref[k][j][1] = (uint8_t)((j * 17 + k * 7 + 6) % 31);
                M[j * 2 * PADDED + k] = Mref[k][j][0];
                M[j * 2 * PADDED + PADDED + k] = Mref[k][j][1];
            }

        gf31_t2_matrix_vec_mul(v, M, C, VCOLS, MCOLS, MCOLS);

        for (int j = 0; j < MCOLS; ++j) {
            uint8_t prod[2], sum[2] = {0, 0};
            for (int k = 0; k < VCOLS; ++k) {
                ref_t2_mul(vref[k], Mref[k][j], prod);
                sum[0] = (uint8_t)ref_mod(sum[0] + prod[0]);
                sum[1] = (uint8_t)ref_mod(sum[1] + prod[1]);
            }
            CHECK(C[j] == sum[0]);
            CHECK(C[MCOLS + j] == sum[1]);
        }
    }

    puts("gf31: field, tower and vector operators passed");
    return 0;
}
