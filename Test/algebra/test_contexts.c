#include "polyvec.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "algebra check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

static int32_t mod(int64_t a, int32_t q)
{
    int64_t r = a % q;
    return (int32_t)(r < 0 ? r + q : r);
}

static int32_t power(int32_t a, uint32_t exponent, int32_t q)
{
    int64_t result = 1;
    while (exponent != 0) {
        if (exponent & 1) result = result * a % q;
        a = (int32_t)((int64_t)a * a % q);
        exponent >>= 1;
    }
    return (int32_t)result;
}

static unsigned int bitreverse(unsigned int value, unsigned int width)
{
    unsigned int result = 0;
    for (unsigned int i = 0; i < width; ++i) {
        result = 2 * result + (value & 1);
        value >>= 1;
    }
    return result;
}

static uint32_t random_word(void)
{
    static uint32_t state = UINT32_C(0x13375493);
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

/* Returns a value-owned plan whose field and ring locals no longer exist. */
static int make_plan(ntt_plan *out, int32_t q, int layers)
{
    fp32_ctx field;
    ring_ctx ring;
    CHECK(fp32_init(&field, q) == 0);
    CHECK(ring_init_negacyclic(&ring, &field, 256) == 0);
    CHECK(ntt_plan_init_256(out, &ring, layers) == 0);
    return 0;
}

static int check_invalid_contexts(void)
{
    fp32_ctx field;
    ring_ctx ring, invalid_ring;
    unsigned char saved_ring[sizeof(ring)];
    ntt_plan plan;
    unsigned char saved_plan[sizeof(plan)];
    CHECK(fp32_init(&field, 3329) == 0);
    CHECK(ring_init_negacyclic(&ring, &field, 256) == 0);
    memcpy(saved_ring, &ring, sizeof(ring));
    CHECK(ring_init_negacyclic(NULL, &field, 256) == -1);
    CHECK(ring_init_negacyclic(&ring, NULL, 256) == -1);
    CHECK(ring_init_negacyclic(&ring, &field, 0) == -1);
    CHECK(ring_init_negacyclic(&ring, &field, -1) == -1);
    for (int i = 0; i < 4; ++i) {
        fp32_ctx invalid_field = field;
        if (i == 0) invalid_field.q = 9;
        if (i == 1) invalid_field.qinv ^= 1;
        if (i == 2) invalid_field.mont_r ^= 1;
        if (i == 3) invalid_field.mont_r2 ^= 1;
        CHECK(ring_init_negacyclic(&ring, &invalid_field, 256) == -1);
        CHECK(memcmp(&ring, saved_ring, sizeof(ring)) == 0);
    }
    CHECK(ntt_plan_init_256(&plan, &ring, 7) == 0);
    memcpy(saved_plan, &plan, sizeof(plan));
    CHECK(ntt_plan_init_256(NULL, &ring, 7) == -1);
    CHECK(ntt_plan_init_256(&plan, NULL, 7) == -1);
    const int invalid_layers[] = {INT_MIN, -1, 0, 6, 8, 9, INT_MAX};
    for (size_t i = 0; i < sizeof(invalid_layers) / sizeof(invalid_layers[0]); ++i) {
        CHECK(ntt_plan_init_256(&plan, &ring, invalid_layers[i]) == -1);
        CHECK(memcmp(&plan, saved_plan, sizeof(plan)) == 0);
    }
    for (int i = 0; i < 5; ++i) {
        invalid_ring = ring;
        if (i == 0) invalid_ring.n = 128;
        if (i == 1) invalid_ring.modulus = (ring_modulus)0;
        if (i == 2) invalid_ring.field.qinv ^= 1;
        if (i == 3) invalid_ring.field.mont_r ^= 1;
        if (i == 4) invalid_ring.field.mont_r2 ^= 1;
        CHECK(ntt_plan_init_256(&plan, &invalid_ring, 7) == -1);
        CHECK(memcmp(&plan, saved_plan, sizeof(plan)) == 0);
    }
    CHECK(fp32_init(&field, 17) == 0);
    CHECK(ring_init_negacyclic(&invalid_ring, &field, 256) == 0);
    CHECK(ntt_plan_init_256(&plan, &invalid_ring, 7) == -1);
    CHECK(memcmp(&plan, saved_plan, sizeof(plan)) == 0);
    CHECK(ring_init_negacyclic(&ring, &field, 3) == 0);
    CHECK(ring.n == 3); /* The mathematical ring is not restricted to NTT sizes. */
    return 0;
}

static void schoolbook(poly out, const poly a, const poly b, int32_t q)
{
    int64_t sum[256] = {0};
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 256; ++j) {
            int degree = i + j;
            int64_t product = (int64_t)a[i] * b[j];
            if (degree < 256) sum[degree] += product;
            else sum[degree - 256] -= product;
        }
    }
    for (int i = 0; i < 256; ++i) out[i] = mod(sum[i], q);
}

static int equal_mod(const poly a, const poly b, int32_t q)
{
    for (int i = 0; i < 256; ++i) {
        if (mod(a[i], q) != mod(b[i], q)) return 0;
    }
    return 1;
}

static void random_poly(poly a, int32_t q)
{
    for (int i = 0; i < 256; ++i) {
        int32_t value = (int32_t)(random_word() % (uint32_t)q);
        a[i] = (random_word() & 1) ? value : -value;
    }
}

static int check_tables(ntt_plan *plan)
{
    int32_t q = plan->ring.field.q;
    unsigned int width = q == 3329 ? 7 : 8;
    int32_t root = q == 3329 ? 17 : 1753;
    int32_t r = mod(INT64_C(4294967296), q);
    CHECK(plan->block_degree == (256 >> plan->layers));
    CHECK(plan->zetas_count == (size_t)1 << plan->layers);
    CHECK(plan->zetas[0] == 0);
    for (size_t i = 1; i < plan->zetas_count; ++i) {
        int32_t zeta = mod((int64_t)r * power(root, bitreverse((unsigned int)i, width), q), q);
        if (zeta > q / 2) zeta -= q;
        CHECK(plan->zetas[i] == zeta);
    }
    int32_t inv_length = power(1 << plan->layers, (uint32_t)q - 2, q);
    CHECK(plan->inverse_scale == mod((int64_t)mod((int64_t)r * r, q) * inv_length, q));
    return 0;
}

static int check_polynomials(ntt_plan *plan)
{
    int32_t q = plan->ring.field.q;
    poly a, b, ah, bh, actual, expected, alias;
    for (int test = 0; test < 29; ++test) {
        memset(a, 0, sizeof(a));
        memset(b, 0, sizeof(b));
        if (test == 1) a[0] = b[0] = 1;
        if (test == 2) { a[255] = 1; b[1] = 1; }
        if (test == 3 || test == 4) {
            for (int i = 0; i < 256; ++i) {
                a[i] = test == 3 ? q - 1 : (i & 1 ? q / 2 : -q / 2);
                b[i] = test == 3 ? 1 - q : (i & 1 ? 1 - q : q - 1);
            }
        }
        if (test >= 5) { random_poly(a, q); random_poly(b, q); }
        schoolbook(expected, a, b, q);
        memcpy(ah, a, sizeof(a));
        memcpy(bh, b, sizeof(b));
        poly_ntt(ah, plan);
        poly_ntt(bh, plan);
        for (int i = 0; i < 256; ++i) {
            CHECK(ah[i] > -(plan->layers + 1) * q && ah[i] < (plan->layers + 1) * q);
        }
        memcpy(actual, ah, sizeof(actual));
        poly_invntt(actual, plan);
        for (int i = 0; i < 256; ++i) {
            CHECK(mod(actual[i], q) == mod((int64_t)a[i] * plan->ring.field.mont_r, q));
        }
        poly_mul_ntt(actual, ah, bh, plan);
        memcpy(alias, ah, sizeof(alias));
        poly_mul_ntt(alias, alias, bh, plan);
        CHECK(memcmp(alias, actual, sizeof(alias)) == 0);
        memcpy(alias, bh, sizeof(alias));
        poly_mul_ntt(alias, ah, alias, plan);
        CHECK(memcmp(alias, actual, sizeof(alias)) == 0);
        for (int i = 0; i < 256; ++i) CHECK(actual[i] >= 0 && actual[i] < q);
        poly_invntt(actual, plan);
        CHECK(equal_mod(actual, expected, q));
        if (test == 2) CHECK(mod(actual[0], q) == q - 1);
    }
    /* Inverse's normalization accepts all int32_t representatives. */
    for (int i = 0; i < 256; ++i) {
        actual[i] = i & 1 ? INT32_MAX : INT32_MIN;
        expected[i] = mod(actual[i], q);
    }
    poly_invntt(actual, plan);
    poly_invntt(expected, plan);
    CHECK(equal_mod(actual, expected, q));
    return 0;
}

static int check_matrices(ntt_plan *plan)
{
    int32_t q = plan->ring.field.q;
    poly matrix[6], vector[3], left[2];
    poly mh[6], vh[3], lh[2], output[3], expected, term;
    for (int i = 0; i < 6; ++i) random_poly(matrix[i], q);
    for (int i = 0; i < 3; ++i) random_poly(vector[i], q);
    for (int i = 0; i < 2; ++i) random_poly(left[i], q);
    memcpy(mh, matrix, sizeof(mh));
    memcpy(vh, vector, sizeof(vh));
    memcpy(lh, left, sizeof(lh));
    ntt_matrix(mh, 2, 3, plan);
    polyvec_ntt(vh, 3, plan);
    ntt_vector(lh, 2, plan);
    polyvec_matrix_pointwise(output, mh, 2, 3, vh, plan);
    polyvec_invntt(output, 2, plan);
    for (int row = 0; row < 2; ++row) {
        memset(expected, 0, sizeof(expected));
        for (int col = 0; col < 3; ++col) {
            schoolbook(term, matrix[3 * row + col], vector[col], q);
            for (int i = 0; i < 256; ++i) expected[i] = mod((int64_t)expected[i] + term[i], q);
        }
        CHECK(equal_mod(output[row], expected, q));
    }
    polyvec_matrix_transpose_pointwise(output, mh, 2, 3, lh, plan);
    polyvec_invntt(output, 3, plan);
    for (int col = 0; col < 3; ++col) {
        memset(expected, 0, sizeof(expected));
        for (int row = 0; row < 2; ++row) {
            schoolbook(term, matrix[3 * row + col], left[row], q);
            for (int i = 0; i < 256; ++i) expected[i] = mod((int64_t)expected[i] + term[i], q);
        }
        CHECK(equal_mod(output[col], expected, q));
    }
    memset(expected, 0, sizeof(expected));
    polyvec_pointwise_acc(output[0], NULL, NULL, 0, plan);
    CHECK(memcmp(output[0], expected, sizeof(expected)) == 0);
    polyvec_matrix_transpose_pointwise(output, NULL, 0, 3, NULL, plan);
    for (int i = 0; i < 3; ++i) CHECK(memcmp(output[i], expected, sizeof(expected)) == 0);
    return 0;
}

static int check_large_accumulation(ntt_plan *plan)
{
    /* A long dot product must not accumulate lazy residues in int32_t. */
    const int count = 2048;
    int32_t q = plan->ring.field.q;
    poly *u = calloc((size_t)count, sizeof(poly));
    poly *v = calloc((size_t)count, sizeof(poly));
    CHECK(u != NULL && v != NULL);
    poly one, result, transposed[1];
    for (int i = 0; i < 256; ++i) {
        u[0][i] = (int32_t)(random_word() % (uint32_t)q);
        v[0][i] = (int32_t)(random_word() % (uint32_t)q);
    }
    for (int i = 1; i < count; ++i) {
        memcpy(u[i], u[0], sizeof(poly));
        memcpy(v[i], v[0], sizeof(poly));
    }
    poly_mul_ntt(one, u[0], v[0], plan);
    polyvec_pointwise_acc(result, u, v, count, plan);
    polyvec_matrix_transpose_pointwise(transposed, u, count, 1, v, plan);
    for (int i = 0; i < 256; ++i) {
        CHECK(result[i] == mod((int64_t)count * one[i], q));
        CHECK(transposed[0][i] == result[i]);
    }
    free(u);
    free(v);
    return 0;
}

int main(void)
{
    const int32_t primes[] = {3329, 8380417, 8380417};
    const int layers[] = {7, 7, 8};
    CHECK(check_invalid_contexts() == 0);
    for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); ++i) {
        ntt_plan original, copied;
        CHECK(make_plan(&original, primes[i], layers[i]) == 0);
        copied = original;
        memset(&original, 0, sizeof(original));
        CHECK(check_tables(&copied) == 0);
        CHECK(check_polynomials(&copied) == 0);
        CHECK(check_matrices(&copied) == 0);
        CHECK(check_large_accumulation(&copied) == 0);
    }
    puts("algebra: contexts, tables, NTT scaling, aliasing and matrix products passed");
    return 0;
}
