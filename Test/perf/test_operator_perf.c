#include "test.h"
#include "polyvec.h"
#include "keccak.h"
#include "gf16.h"
#include "gf16_matrix.h"

#include <stdio.h>

#define ROUNDS 10000
#define RUN_CASE(id, label, rounds, sink, ...) do { \
    if (strcmp(argv[1], id) == 0) { \
        test_perf_mark start = test_perf_begin(); \
        uint32_t checksum = 0; \
        for (int i = 0; i < (rounds); ++i) { __VA_ARGS__; checksum ^= (uint32_t)(sink); } \
        volatile uint32_t benchmark_output = checksum; \
        (void)benchmark_output; \
        test_perf_print(label, test_perf_end(start), (uint64_t)(rounds)); \
        return 0; \
    } \
} while (0)

static int rng_cb(void *ctx, uint8_t *out, size_t n) {
    (void)ctx;
    memset(out, 0x5a, n);
    return 0;
}

int main(int argc, char **argv) {
    fp32_ctx f;
    ring_ctx r;
    ntt_plan p;
    poly a, b, c, v[4], u[4], m[4 * 4];
    uint8_t bytes[168] = {0}, digest[64];
    KECCAK1600_CTX h;
    gf16_t gf_a[16], gf_b[16], gf_c[16];
    gf16_t gf_a16[256], gf_b16[256], gf_c16[256];

    if (argc != 2) {
        fprintf(stderr, "Pass one operator benchmark name; see BUILDING.txt.\n");
        return 2;
    }
    if (fp32_init(&f, 8380417) != 0 || ring_init_negacyclic(&r, &f, 256) != 0 ||
        ntt_plan_init_256(&p, &r, 8) != 0) {
        fprintf(stderr, "benchmark context initialization failed\n");
        return 1;
    }
    for (int i = 0; i < 256; ++i) {
        a[i] = i;
        b[i] = 255 - i;
        gf_a16[i] = (gf16_t)(i & 0x0f);
        gf_b16[i] = (gf16_t)((i * 7 + 3) & 0x0f);
    }
    for (int i = 0; i < 4; ++i) {
        memcpy(v[i], a, sizeof(a));
        memcpy(u[i], a, sizeof(a));
    }
    for (int i = 0; i < 4 * 4; ++i) memcpy(m[i], b, sizeof(b));
    memcpy(gf_a, gf_a16, sizeof(gf_a));
    memcpy(gf_b, gf_b16, sizeof(gf_b));

    RUN_CASE("gf16_add", "GF16 add", ROUNDS, gf_c[0],
             gf_c[0] = gf16_add(gf_a[i & 15], gf_b[(i + 1) & 15]));
    RUN_CASE("gf16_sub", "GF16 subtract", ROUNDS, gf_c[0],
             gf_c[0] = gf16_sub(gf_a[i & 15], gf_b[(i + 1) & 15]));
    RUN_CASE("gf16_mul", "GF16 multiply", ROUNDS, gf_c[0],
             gf_c[0] = gf16_mul(gf_a[i & 15], gf_b[(i + 1) & 15]));
    RUN_CASE("gf16_inv", "GF16 inverse", ROUNDS, gf_c[0],
             gf_c[0] = gf16_inv(gf_a[(i % 15) + 1]));
    RUN_CASE("gf16_matrix_mul_4", "GF16 matrix multiply (4x4)", 10000, gf_c[0],
             gf16_matrix_mul(gf_c, gf_a, gf_b, 4, 4, 4));
    RUN_CASE("gf16_matrix_mul_16", "GF16 matrix multiply (16x16)", 1000, gf_c16[0],
             gf16_matrix_mul(gf_c16, gf_a16, gf_b16, 16, 16, 16));
    RUN_CASE("gf16_matrix_mul_add_4", "GF16 matrix multiply-add (4x4)", 10000, gf_c[0],
             gf16_matrix_mul_add(gf_c, gf_a, gf_b, 4, 4, 4));
    RUN_CASE("gf16_matrix_det_4", "GF16 determinant (4x4)", 10000, gf_c[0],
             gf16_matrix_determinant(&gf_c[0], gf_a, 4));
    RUN_CASE("gf16_matrix_det_5", "GF16 determinant (5x5)", 1000, gf_c[0],
             gf16_matrix_determinant(&gf_c[0], gf_a16, 5));

    RUN_CASE("fp32_add", "fp32 add", ROUNDS, c[i % 256],
             c[i % 256] = fp32_add(a[i % 256] + i, b[i % 256], &f));
    RUN_CASE("fp32_sub", "fp32 sub", ROUNDS, c[i % 256],
             c[i % 256] = fp32_sub(a[i % 256] + i, b[i % 256], &f));
    RUN_CASE("fp32_mul", "fp32 multiply", ROUNDS, c[i % 256],
             c[i % 256] = fp32_mul(a[i % 256] + i, b[i % 256], &f));
    RUN_CASE("fp32_mont_reduce", "fp32 Montgomery reduce", ROUNDS, c[i % 256],
             c[i % 256] = fp32_mont_reduce((int64_t)(a[i % 256] + i) * b[i % 256], &f));
    RUN_CASE("fp32_inv", "fp32 inverse", ROUNDS, c[i % 256],
             fp32_inv(&c[i % 256], a[(i % 255) + 1] + 1, &f));

    RUN_CASE("poly_add_256", "poly add (256)", ROUNDS, c[0],
             { a[0] = i & 0x7fff; POLY_ADD(c, a, b, 256); });
    RUN_CASE("poly_sub_256", "poly sub (256)", ROUNDS, c[0],
             { a[0] = i & 0x7fff; POLY_SUB(c, a, b, 256); });
    RUN_CASE("poly_reduce_256", "poly reduce (256)", ROUNDS, c[0],
             a[0] = i & 0x7fff, memcpy(c, a, sizeof(a)), poly_reduce_mod(c, 256, f.q));
    RUN_CASE("poly_ntt_256", "poly NTT (256)", 1000, c[0],
             a[0] = i & 0x7fff, memcpy(c, a, sizeof(a)), poly_ntt(c, &p));
    RUN_CASE("poly_invntt_256", "poly inverse NTT (256)", 1000, c[0],
             a[0] = i & 0x7fff, memcpy(c, a, sizeof(a)), poly_ntt(c, &p), poly_invntt(c, &p));
    RUN_CASE("poly_mul_ntt_256", "poly NTT pointwise (256)", 1000, c[0],
             a[0] = i & 0x7fff, poly_mul_ntt(c, a, b, &p));

    RUN_CASE("polyvec_add_4", "polyvec add (4x256)", ROUNDS, v[0][0],
             v[0][0] = i & 0x7fff, polyvec_add(v, (const poly *)u, (const poly *)u, 4));
    RUN_CASE("polyvec_ntt_4", "polyvec NTT (4x256)", 1000, v[0][0],
             v[0][0] = i & 0x7fff, polyvec_ntt(v, 4, &p));
    RUN_CASE("polyvec_inner_4", "polyvec inner product (4)", 1000, c[0],
             v[0][0] = i & 0x7fff, polyvec_pointwise_acc(c, (const poly *)u, (const poly *)u, 4, &p));
    RUN_CASE("polyvec_matrix_4", "polyvec matrix product (4x4)", 1000, v[0][0],
             polyvec_matrix_pointwise(v, (const poly *)m, 4, 4, (const poly *)u, &p));

    RUN_CASE("sample_cbd_256", "CBD sampling (256)", ROUNDS, c[0],
             bytes[i % sizeof(bytes)] = (uint8_t)i, sample_CBD(c, 256, bytes, sizeof(bytes), &p));
    RUN_CASE("randombytes_168", "random bytes (168)", ROUNDS, bytes[i % sizeof(bytes)],
             pqc_randombytes(bytes, sizeof(bytes), rng_cb, NULL));
    RUN_CASE("sha3_256_168", "SHA3-256 (168 bytes)", ROUNDS, digest[0],
             bytes[0] = (uint8_t)i, sha3_init(&h, 0x06, 256), sha3_update(&h, bytes, sizeof(bytes)), sha3_final(digest, &h));

    fprintf(stderr, "unknown operator benchmark: %s\n", argv[1]);
    return 2;
}
