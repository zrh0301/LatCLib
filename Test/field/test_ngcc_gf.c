#include "gf/byte_matrix.h"
#include "gf/composite_matrix.h"
#include "gf/gf16_vector.h"
#include "gf/gf256_matrix.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

static uint8_t get(const uint8_t *v, size_t i, gf_composite_field field)
{
    return field == GF_COMPOSITE_16 ? gf16_packed_get(v, i) : v[i];
}

static void put(uint8_t *v, size_t i, uint8_t value, gf_composite_field field)
{
    if (field == GF_COMPOSITE_16)
        gf16_packed_set(v, i, value);
    else
        v[i] = value;
}

/* Schoolbook multiplication in the tower, reducing the quadratic term. */
static uint8_t reference_mul(uint8_t a, uint8_t b, unsigned bits)
{
    if (bits == 1)
        return a & b;
    const unsigned half = bits / 2;
    const unsigned mask = (1u << half) - 1u;
    const uint8_t d = bits == 8 ? 8 : bits == 4 ? 2 : 1;
    uint8_t coeff[3] = {0};
    for (unsigned i = 0; i < 2; ++i) {
        for (unsigned j = 0; j < 2; ++j)
            coeff[i + j] ^= reference_mul((a >> (half * i)) & mask,
                                           (b >> (half * j)) & mask, half);
    }
    return (uint8_t)((coeff[0] ^ reference_mul(coeff[2], d, half)) |
                     ((coeff[1] ^ coeff[2]) << half));
}

static int test_composite(gf_composite_field field)
{
    const unsigned bits = field == GF_COMPOSITE_16 ? 4 : 8;
    const unsigned q = 1u << bits;
    uint8_t (*mul)(uint8_t, uint8_t) =
        field == GF_COMPOSITE_16 ? gf16_composite_mul : gf256_composite_mul;
    uint8_t (*inv)(uint8_t) =
        field == GF_COMPOSITE_16 ? gf16_composite_inv : gf256_composite_inv;
    void (*madd)(uint8_t *, const uint8_t *, uint8_t, size_t) =
        field == GF_COMPOSITE_16 ? gf16_composite_packed_madd : gf256_composite_madd;
    for (unsigned a = 0; a < q; ++a) {
        CHECK(mul(a, inv(a)) == (a != 0));
        for (unsigned b = 0; b < q; ++b)
            CHECK(mul(a, b) == reference_mul(a, b, bits));
        const uint8_t in[7] = {0xf0, 0x80, 0xac, 0x01, 0x34, 0x95, 0x77};
        uint8_t out[7];
        memset(out, 0x5a, sizeof(out));
        madd(out, in, a, sizeof(out));
        for (size_t i = 0; i < sizeof(in); ++i) {
            uint8_t expected = reference_mul(in[i] & (q - 1), a, bits);
            if (bits == 4)
                expected |= (uint8_t)(reference_mul(in[i] >> 4, a, bits) << 4);
            CHECK(out[i] == (uint8_t)(expected ^ 0x5a));
        }
    }

    const uint8_t l[] = {2, 3, 1, 4, 5, 6};
    const uint8_t u[] = {1, 2, 3, 4, 5, 7};
    const size_t stride = bits == 4 ? 2 : 3;
    uint8_t il[6], iu[6], m[9], im[9], actual[3];
    gf_composite_triangular_inv(il, l, 3, 0, field);
    gf_composite_triangular_inv(iu, u, 3, 1, field);
    gf_composite_triangular_product(m, l, u, 3, 0, field);
    gf_composite_triangular_product(im, il, iu, 3, 1, field);
    for (size_t col = 0; col < 3; ++col) {
        gf_composite_col_matvec(actual, m, im + stride * col, 3, 3, stride, field);
        for (size_t row = 0; row < 3; ++row)
            CHECK(get(actual, row, field) == (row == col));
        if (bits == 4)
            CHECK((actual[1] & 0xf0) == 0);
    }

    const uint8_t a[] = {0, 1, 2, 1, 2, 3, 0, 0, 1};
    const uint8_t x[] = {2, 3, 4};
    uint8_t b[3] = {0}, got[3], work[12];
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            b[i] ^= reference_mul(a[i * 3 + j], x[j], bits);
    CHECK(gf_composite_matrix_solve(got, a, b, 3, work, field) == 0);
    CHECK(memcmp(got, x, sizeof(got)) == 0);
    const uint8_t singular[9] = {0};
    CHECK(gf_composite_matrix_solve(got, singular, b, 3, work, field) == -1);
    CHECK(memcmp(got, x, sizeof(got)) == 0);

    const uint8_t upper[] = {0x13, 0x25, 0x37, 0x49, 0x5b, 0x6d};
    uint8_t scalars[6] = {0}, out[4] = {1, 2, 3, 4};
    uint8_t trans[4] = {0}, tmp[256], folded[3] = {9, 9, 9}, quad[2];
    const size_t elements = bits == 4 ? 2 : 1;
    for (size_t i = 0; i < 2; ++i)
        put(scalars + i * stride, i, 1, field);
    gf_composite_batch_tri_madd(out, upper, scalars, 2, stride, 2, 1, field);
    CHECK(out[0] == (1 ^ upper[0]));
    CHECK(out[1] == (2 ^ upper[1]));
    CHECK(out[2] == 3);
    CHECK(out[3] == (4 ^ upper[2]));
    gf_composite_batch_transpose_madd(trans, scalars, 2, stride, 2, out, 2, 1, field);
    CHECK(memcmp(trans, out, 4) == 0);
    gf_byte_batch_fold_upper_add(folded, trans, 2, 1);
    CHECK(folded[0] == (9 ^ out[0]));
    CHECK(folded[1] == (9 ^ out[1] ^ out[2]));
    CHECK(folded[2] == (9 ^ out[3]));

    memset(scalars, 0, sizeof(scalars));
    put(scalars, 0, 2, field);
    put(scalars, 1, 3, field);
    gf_composite_batch_quad(quad, upper, scalars, 2, 1, tmp, field);
    for (size_t k = 0; k < elements; ++k) {
        CHECK(get(quad, k, field) ==
              (mul(get(upper, k, field), mul(2, 2)) ^
               mul(get(upper + 1, k, field), mul(2, 3)) ^
               mul(get(upper + 2, k, field), mul(3, 3))));
    }
    return 0;
}

static int test_gf256(void)
{
    const uint8_t matrix[] = {0, 5, 1, 7, 1, 3, 0, 0, 9};
    const uint8_t x[] = {7, 8, 9};
    const uint8_t ident[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    const uint8_t singular[9] = {0};
    uint8_t b[3], inv[9], product[9];

    gf256_matrix_mul_vec(b, matrix, x, 3, 3);
    CHECK(gf256_matrix_mul(product, matrix, x, 3, 3, 1) == 0);
    CHECK(memcmp(product, b, 3) == 0);
    CHECK(gf256_matrix_quadratic(x, matrix, 3) == gf256_vector_dot(x, b, 3));
    CHECK(gf256_matrix_solve(b, matrix, b, 3) == 0);
    CHECK(memcmp(b, x, sizeof(b)) == 0);
    CHECK(gf256_matrix_solve(b, singular, x, 3) == -1);
    CHECK(gf256_matrix_solve(b, matrix, x, SIZE_MAX) == -1);
    CHECK(memcmp(b, x, sizeof(b)) == 0);

    memcpy(inv, matrix, sizeof(inv));
    CHECK(gf256_matrix_inv(inv, inv, 3) == 0);
    CHECK(gf256_matrix_mul(product, matrix, inv, 3, 3, 3) == 0);
    CHECK(memcmp(product, ident, sizeof(ident)) == 0);
    CHECK(gf256_matrix_inv(product, singular, 3) == -1);
    CHECK(gf256_matrix_inv(product, matrix, SIZE_MAX) == -1);
    CHECK(memcmp(product, ident, sizeof(ident)) == 0);

    const uint8_t upper[] = {1, 2, 3, 4, 5, 6};
    const uint8_t square[] = {1, 2, 3, 0, 4, 5, 0, 0, 6};
    const uint8_t transpose[] = {1, 0, 0, 2, 4, 0, 3, 5, 6};
    gf_byte_unpack_upper(inv, upper, 3);
    CHECK(memcmp(inv, square, sizeof(square)) == 0);
    gf_byte_transpose(product, inv, 3, 3);
    CHECK(memcmp(product, transpose, sizeof(transpose)) == 0);
    gf_byte_fold_upper(product, 3);
    CHECK(memcmp(product, square, sizeof(square)) == 0);
    gf_byte_pack_upper(inv, product, 3);
    CHECK(memcmp(inv, upper, sizeof(upper)) == 0);
    return 0;
}

int main(void)
{
    return test_gf256() || test_composite(GF_COMPOSITE_16) ||
           test_composite(GF_COMPOSITE_256);
}
