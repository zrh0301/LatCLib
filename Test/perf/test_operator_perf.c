#include "test.h"
#include "poly/fp32.h"
#include "poly/fp16.h"
#include "poly/poly16.h"
#include "poly/ntt16.h"
#include "poly/ring4171.h"
#include "poly/polyvec.h"
#include "sym/keccak.h"
#include "sym/sm3.h"
#include "tree/merkle.h"
#include "util/combinatorics.h"
#include "gf/gf16.h"
#include "gf/gf16_matrix.h"
#include "gf/gf31.h"
#include "gf/gf31_matrix.h"
#include "gf/mod_pow2.h"
#include "gf/byte_matrix.h"
#include "gf/composite_matrix.h"
#include "gf/gf256_matrix.h"
#include "util/bytes.h"
#include "util/u256.h"
#include "util/ct.h"

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

static void benchmark_leaf(uint8_t *out, uint32_t index, void *ctx)
{
    (void)ctx;
    uint8_t encoded[4];
    pqc_store_be(encoded, sizeof(encoded), index);
    sm3_hash(out, encoded, sizeof(encoded));
}

static void benchmark_parent(uint8_t *out, const uint8_t *children,
                              uint32_t height, uint32_t index, void *ctx)
{
    (void)ctx;
    uint8_t input[72];
    pqc_store_be(input, 4, height);
    pqc_store_be(input + 4, 4, index);
    memcpy(input + 8, children, 64);
    sm3_hash(out, input, sizeof(input));
}

int main(int argc, char **argv) {
    fp32_ctx f;
    fp16_ctx f16;
    ring_ctx r;
    ntt_plan p;
    poly a, b, c, v[4], u[4], m[4 * 4];
    uint8_t bytes[168] = {0}, digest[64];
    KECCAK1600_CTX h;
    gf16_t gf_a[16], gf_b[16], gf_c[16];
    gf16_t gf_a16[256], gf_b16[256], gf_c16[256];
    gf31_t q31_a[64], q31_b[64], q31_out[64];
    uint8_t q31_t2_a[4 * 2 * 32] = {0}, q31_t2_b[4 * 2 * 32] = {0};
    uint8_t q31_t2_c[4 * 2 * 4] = {0};
    uint8_t q31_t2_v[2 * 32] = {0}, q31_t2_m[4 * 2 * 32] = {0};
    uint8_t q31_t2_vc[2 * 4] = {0};
    int16_t p16_a[512], p16_b[512], p16_c[512];
    int32_t q41_a[PQC_Q4171_N], q41_b[PQC_Q4171_N], q41_c[PQC_Q4171_N];
    uint8_t cbd_buf[768], packed_poly[768];
    uint8_t gm[16 * 16] = {0}, gi[16 * 16], gs[16], work[16 * 17];
    uint8_t batch[16 * 16 * 8] = {0}, batch_out[16 * 16 * 8] = {0};
    uint8_t quad_work[256 * 8], triangle[16 * 17 / 2];
    pqc_u256 large = {{UINT64_C(0x1234567890abcd), UINT64_MAX, 4, 0}};
    pqc_u256 integer, mod = {{0, 0, 0, 1}};
    uint8_t tree_root[32], tree_auth[12 * 32], proof_leaves[3 * 32];
    pqc_merkle_node tree_positions[12];
    const uint32_t selected[] = {0, 3, 12};
    size_t tree_auth_count;
    pqc_u256 residues[4][256], modulus, rank;
    uint32_t subset[16];

    if (argc != 2) {
        fprintf(stderr, "Pass one operator benchmark name; see BUILDING.txt.\n");
        return 2;
    }
    if (fp32_init(&f, 8380417) != 0 || fp16_init(&f16, 3329) != 0 ||
        ring_init_negacyclic(&r, &f, 256) != 0 ||
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

    for (int i = 0; i < 64; ++i) {
        q31_a[i] = (gf31_t)(i % 31);
        q31_b[i] = (gf31_t)((i * 7 + 3) % 31);
    }
    for (size_t i = 0; i < sizeof(q31_t2_a); ++i) {
        q31_t2_a[i] = (uint8_t)((i * 3 + 1) % 31);
        q31_t2_b[i] = (uint8_t)((i * 5 + 2) % 31);
        q31_t2_m[i] = (uint8_t)((i * 7 + 4) % 31);
    }
    for (int i = 0; i < 4; ++i) {
        q31_t2_v[i] = (uint8_t)(i + 1);
        q31_t2_v[32 + i] = (uint8_t)(i + 2);
    }

    for (int i = 0; i < 512; ++i) {
        p16_a[i] = (int16_t)(i % 3329);
        p16_b[i] = (int16_t)((i * 7 + 3) % 3329);
        q41_a[i] = (int32_t)((i * 97 + 11) % PQC_Q4171);
        q41_b[i] = (int32_t)((i * 53 + 7) % PQC_Q4171);
    }
    for (size_t i = 0; i < sizeof(cbd_buf); ++i) {
        cbd_buf[i] = (uint8_t)(i * 37 + 11);
    }

    for (size_t row = 0; row < 16; ++row) {
        for (size_t col = row; col < 16; ++col) {
            gm[row * 16 + col] = row == col ? 1 : (uint8_t)((row + 3 * col) & 15);
        }
    }
    for (size_t i = 0; i < sizeof(batch); ++i)
        batch[i] = (uint8_t)(i * 37 + 11);
    for (size_t i = 0; i < sizeof(triangle); ++i)
        triangle[i] = (uint8_t)(i % 15 + 1);

    if (pqc_merkle_auth_nodes(tree_positions, &tree_auth_count, 12, selected, 3, 13) != 0 ||
        pqc_merkle_build(tree_root, tree_auth, 32, 13, tree_positions, tree_auth_count,
                          benchmark_leaf, NULL, benchmark_parent, NULL) != 0)
        return 1;
    for (size_t j = 0; j < 3; ++j)
        benchmark_leaf(proof_leaves + j * 32, selected[j], NULL);
    pqc_u256_set(&modulus, 65521);
    uint64_t weight = 1;
    for (size_t position = 4; position-- > 0;) {
        for (unsigned value = 0; value < 256; ++value)
            pqc_u256_set(&residues[position][value], weight * value % 65521);
        weight = weight * 256 % 65521;
    }

    RUN_CASE("merkle_treehash_h4", "Merkle tree and path (16 leaves)", 1000, tree_root[0],
             pqc_merkle_treehash(tree_root, tree_auth, 32, 4, 5, 0,
                                 benchmark_leaf, NULL, benchmark_parent, NULL));
    if (strcmp(argv[1], "merkle_root_path_h4") == 0 &&
        pqc_merkle_treehash(tree_root, tree_auth, 32, 4, 0, 0,
                            benchmark_leaf, NULL, benchmark_parent, NULL) != 0)
        return 1;
    RUN_CASE("merkle_root_path_h4", "Merkle root from path (height 4)", 10000, tree_root[0],
             pqc_merkle_root_from_path(tree_root, proof_leaves, tree_auth, 32, 4, 0, 0,
                                       benchmark_parent, NULL));
    RUN_CASE("merkle_auth_nodes_13", "Merkle multiproof positions (3 of 13)", ROUNDS, tree_auth_count,
             pqc_merkle_auth_nodes(tree_positions, &tree_auth_count, 12, selected, 3, 13));
    RUN_CASE("merkle_build_13", "Merkle shortened tree (13 leaves)", 1000, tree_root[0],
             pqc_merkle_build(tree_root, tree_auth, 32, 13, tree_positions, tree_auth_count,
                               benchmark_leaf, NULL, benchmark_parent, NULL));
    RUN_CASE("merkle_root_proof_13", "Merkle root from multiproof (3 of 13)", 10000, tree_root[0],
             pqc_merkle_root_from_proof(tree_root, selected, proof_leaves, 3, tree_auth,
                                        tree_auth_count, 32, 13, benchmark_parent, NULL));
    RUN_CASE("u256_binomial_4721_27", "256-bit binomial coefficient C(4721,27)", ROUNDS, integer.limb[0],
             pqc_u256_binomial(&integer, 4721, 27));
    RUN_CASE("combination_unrank_128_16", "Combination unranking (16 of 128)", 1000, subset[0],
             pqc_u256_set(&rank, (uint64_t)i),
             pqc_combination_unrank(subset, 128, 16, &rank, NULL));
    RUN_CASE("u256_reduce_bytes_4", "256-bit positional modular reduction (4 bytes)", ROUNDS,
             integer.limb[0], bytes[0] = (uint8_t)i,
             pqc_u256_reduce_bytes_table(&integer, bytes, 4, residues, &modulus));

    RUN_CASE("gf256_dot_256", "GF256 dot product (256)", ROUNDS, gs[0],
             gs[0] = gf256_vector_dot(batch, batch + (i & 63), 256));
    RUN_CASE("gf256_solve_16", "GF256 solve (16)", 2000, gs[0],
             gf256_matrix_solve(gs, gm, batch + (i & 63), 16));
    RUN_CASE("gf256_inv_16", "GF256 matrix inverse (16)", 2000, gi[0],
             gf256_matrix_inv(gi, gm, 16));
    RUN_CASE("gf256_matvec_16", "GF256 matrix-vector product (16)", ROUNDS, gs[0],
             gf256_matrix_mul_vec(gs, gm, batch + (i & 63), 16, 16));
    RUN_CASE("gf256_quadratic_16", "GF256 quadratic form (16)", ROUNDS, gs[0],
             gs[0] = gf256_matrix_quadratic(batch + (i & 63), gm, 16));
    RUN_CASE("gf_byte_upper_16", "GF byte upper-triangle conversion (16)", ROUNDS, gi[0],
             gf_byte_fold_upper(gm, 16), gf_byte_pack_upper(gi, gm, 16),
             gf_byte_unpack_upper(gm, gi, 16), gf_byte_transpose(gi, gm, 16, 16));
    RUN_CASE("gf_byte_batch_upper_16", "GF byte batched upper triangle (16 x 8 bytes)", ROUNDS,
             batch_out[0], gf_byte_batch_fold_upper_add(batch_out, batch, 16, 8));

    RUN_CASE("gf16_composite_mul", "GF16 composite multiply", ROUNDS, gs[0],
             gs[0] = gf16_composite_mul((uint8_t)(i & 15), 7));
    RUN_CASE("gf256_composite_mul", "GF256 composite multiply", ROUNDS, gs[0],
             gs[0] = gf256_composite_mul((uint8_t)i, 123));
    RUN_CASE("gf16_composite_inv", "GF16 composite inverse", ROUNDS, gs[0],
             gs[0] = gf16_composite_inv((uint8_t)(i & 15)));
    RUN_CASE("gf256_composite_inv", "GF256 composite inverse", ROUNDS, gs[0],
             gs[0] = gf256_composite_inv((uint8_t)i));
    RUN_CASE("gf16_composite_madd_256", "GF16 composite packed multiply-add (256 bytes)", ROUNDS,
             batch_out[0], gf16_composite_packed_madd(batch_out, batch, (uint8_t)(i & 15), 256));
    RUN_CASE("gf256_composite_madd_256", "GF256 composite multiply-add (256 bytes)", ROUNDS,
             batch_out[0], gf256_composite_madd(batch_out, batch, (uint8_t)i, 256));
    RUN_CASE("gf16_composite_solve_16", "GF16 composite solve (16)", 2000, gs[0],
             gf_composite_matrix_solve(gs, gm, gf_a16, 16, work, GF_COMPOSITE_16));
    RUN_CASE("gf16_composite_triangular_16", "GF16 composite triangular inverse and product (16)",
             2000, gi[0], gf_composite_triangular_inv(gi, triangle, 16, 1, GF_COMPOSITE_16),
             gf_composite_triangular_product(batch_out, triangle, gi, 16, 0, GF_COMPOSITE_16));
    RUN_CASE("gf16_composite_matvec_16", "GF16 composite packed matrix-vector product (16)",
             ROUNDS, gs[0], gf_composite_col_matvec(gs, gm, batch + (i & 63), 16, 16, 8,
                                                    GF_COMPOSITE_16));
    RUN_CASE("gf16_composite_batch_16", "GF16 composite batched matrix products (16 x 8 bytes)",
             1000, batch_out[0], gf_composite_batch_tri_madd(batch_out, batch, gm, 16, 8, 16,
                                                            8, GF_COMPOSITE_16),
             gf_composite_batch_transpose_madd(batch, gm, 16, 8, 16, batch_out, 16, 8,
                                               GF_COMPOSITE_16));
    RUN_CASE("gf256_composite_quad_16", "GF256 composite batched quadratic form (16 x 8 bytes)",
             ROUNDS, gs[0], gf_composite_batch_quad(gs, batch, batch + (i & 63), 16, 8,
                                                    quad_work, GF_COMPOSITE_256));

    RUN_CASE("u256_mod_add", "256-bit modular addition", ROUNDS, large.limb[0],
             pqc_u256_set(&integer, (uint64_t)i + 1), pqc_u256_add_mod_assign(&large, &integer, &mod));
    RUN_CASE("u256_mul_div_u32", "256-bit multiply and exact divide", ROUNDS, integer.limb[0],
             integer = large, integer.limb[0] ^= (uint64_t)i,
             pqc_u256_mul_u32(&integer, (uint32_t)i + 2),
             pqc_u256_div_u32_exact(&integer, (uint32_t)i + 2));
    RUN_CASE("bytes_be_64", "64-bit big-endian encode and decode", ROUNDS, mod.limb[0],
             large.limb[0] += (uint64_t)i, pqc_store_be(bytes, 8, large.limb[0]),
             mod.limb[0] = pqc_load_be(bytes, 8));

    RUN_CASE("gf16_add", "GF16 add", ROUNDS, gf_c[0],
             gf_c[0] = gf16_add(gf_a[i & 15], gf_b[(i + 1) & 15]));
    RUN_CASE("gf16_sub", "GF16 subtract", ROUNDS, gf_c[0],
             gf_c[0] = gf16_sub(gf_a[i & 15], gf_b[(i + 1) & 15]));
    RUN_CASE("gf16_mul", "GF16 multiply", ROUNDS, gf_c[0],
             gf_c[0] = gf16_mul(gf_a[i & 15], gf_b[(i + 1) & 15]));
    RUN_CASE("gf16_inv", "GF16 inverse", ROUNDS, gf_c[0],
             gf_c[0] = gf16_inv(gf_a[(i % 15) + 1]));
    RUN_CASE("gf31_add", "GF31 add", ROUNDS, q31_out[0],
             q31_out[0] = gf31_add(q31_a[i & 63], q31_b[(i + 1) & 63]));
    RUN_CASE("gf31_mul", "GF31 multiply", ROUNDS, q31_out[0],
             q31_out[0] = gf31_mul(q31_a[i & 63], q31_b[(i + 1) & 63]));
    RUN_CASE("gf31_inv", "GF31 inverse", ROUNDS, q31_out[0],
             q31_out[0] = gf31_inv(q31_a[(i % 30) + 1]));
    RUN_CASE("gf31_vector_add_64", "GF31 vector add (64)", ROUNDS, q31_out[0],
             gf31_vector_add(q31_out, q31_b, 64));
    RUN_CASE("gf31_vector_dot_64", "GF31 vector dot product (64)", ROUNDS, q31_out[0],
             q31_out[0] = gf31_vector_dot(q31_a, q31_b, 64));
    RUN_CASE("gf31_t2_matrix_mul_4", "GF31^2 planar matrix multiply (4x4)",
             10000, q31_t2_c[0],
             gf31_t2_matrix_mul(q31_t2_a, q31_t2_b, q31_t2_c, 4, 4, 4, 4));
    RUN_CASE("gf31_t2_matvec_4", "GF31^2 planar matrix-vector product (4)",
             ROUNDS, q31_t2_vc[0],
             gf31_t2_matrix_vec_mul(q31_t2_v, q31_t2_m, q31_t2_vc, 4, 4, 4));
    RUN_CASE("mod_pow2_1024", "reduce modulo 2^10", ROUNDS, q31_out[0],
             q31_out[0] = (gf31_t)pqc_mod_pow2((uint32_t)i * 1234567u, 10));
    RUN_CASE("mod_pow2_round_1024", "round 2^10 coefficient to 4 bits", ROUNDS, q31_out[0],
             q31_out[0] = (gf31_t)pqc_round_pow2((uint16_t)(i * 37), 10, 4));
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

    RUN_CASE("fp16_mont_reduce", "fp16 Montgomery reduce", ROUNDS, c[i % 256],
             c[i % 256] = fp16_montgomery_reduce((int32_t)(a[i % 256] + i) * b[i % 256], &f16));
    RUN_CASE("fp16_barrett_reduce", "fp16 Barrett reduce", ROUNDS, c[i % 256],
             c[i % 256] = fp16_barrett_reduce((int16_t)((a[i % 256] % 3329) - 3000), &f16));
    RUN_CASE("fp16_caddq2", "fp16 conditional add (twice)", ROUNDS, c[i % 256],
             c[i % 256] = fp16_caddq2((int16_t)((b[i % 256] % 3329) - 6000), &f16));

    RUN_CASE("ct_equal_32", "constant-time equal (32 bytes)", ROUNDS, bytes[0],
             bytes[0] = (uint8_t)i, memcpy(c, a, 32),
             bytes[0] = (uint8_t)pqc_ct_equal((const uint8_t *)c, (const uint8_t *)a, 32));
    RUN_CASE("ct_select_32", "constant-time select (32 bytes)", ROUNDS, bytes[0],
             pqc_ct_select(bytes, (const uint8_t *)a, 32, (uint8_t)(i & 1)));

    RUN_CASE("poly16_add_512", "poly16 add (512)", ROUNDS, p16_c[0],
             p16_a[0] = (int16_t)(i % 3329), fp16_poly_add(p16_c, p16_a, p16_b, 512));
    RUN_CASE("poly16_barrett_512", "poly16 Barrett reduce (512)", ROUNDS, p16_c[0],
             p16_a[0] = (int16_t)(i % 3329), fp16_poly_barrett(p16_a, 512, &f16),
             p16_c[0] = p16_a[0]);
    RUN_CASE("poly16_tomont_512", "poly16 to-Montgomery (512)", ROUNDS, p16_c[0],
             p16_a[0] = (int16_t)(i % 3329), fp16_poly_tomont(p16_a, 512, &f16),
             p16_c[0] = p16_a[0]);
    RUN_CASE("cbd_eta5_512", "CBD eta=5 (512)", ROUNDS, p16_c[0],
             fp16_poly_cbd(p16_c, 512, cbd_buf, 5));
    RUN_CASE("cbd_eta6_512", "CBD eta=6 (512)", ROUNDS, p16_c[0],
             fp16_poly_cbd(p16_c, 512, cbd_buf, 6));

    RUN_CASE("poly16_encode_512", "poly16 encode and decode (512, 12 bits)", ROUNDS, p16_c[0],
             p16_a[0] = (int16_t)(i % 3329), fp16_poly_encode(packed_poly, p16_a, 512, 12),
             fp16_poly_decode(p16_c, packed_poly, 512, 12));
    RUN_CASE("poly16_compress_512", "poly16 compress and decompress (512, 10 bits)", ROUNDS, p16_c[0],
             p16_a[0] = (int16_t)(i % 3329), fp16_poly_compress(packed_poly, p16_a, 512, 10, &f16),
             fp16_poly_decompress(p16_c, packed_poly, 512, 10, &f16));
    RUN_CASE("poly16_repeat_encode_512", "poly16 repeated-bit encode (512)", ROUNDS, p16_c[0],
             cbd_buf[0] = (uint8_t)i, fp16_poly_encode_repeated_bits(p16_c, cbd_buf, 16, 4, 4, &f16));
    RUN_CASE("poly16_repeat_decode_512", "poly16 repeated-bit decode (512)", ROUNDS, packed_poly[0],
             p16_a[0] = (int16_t)(i % 3329),
             fp16_poly_decode_repeated_bits(packed_poly, p16_a, 16, 4, 4, &f16));
    RUN_CASE("poly16_rej_uniform12_512", "poly16 uniform rejection (768 input bytes)", ROUNDS,
             p16_c[0], fp16_poly_rej_uniform12(p16_c, 512, cbd_buf, sizeof(cbd_buf), &f16));

    RUN_CASE("ntt16_128", "128-point NTT", 10000, p16_c[0],
             memcpy(p16_c, p16_a, 128 * sizeof(int16_t)), pqc_ntt128(p16_c));
    RUN_CASE("ntt16_blocks_512", "128-point NTT (4 blocks)", 5000, p16_c[0],
             memcpy(p16_c, p16_a, sizeof(p16_c)), pqc_ntt128_blocks(p16_c, 4));

    RUN_CASE("ntt16_product_512", "128-point NTT split product (4 blocks)", ROUNDS, p16_c[0],
             p16_a[0] = (int16_t)(i % 3329),
             pqc_ntt128_blocks_mont_mul(p16_c, p16_a, p16_b, 4));

    RUN_CASE("q4171_montgomery", "q=4171777 Montgomery reduce", ROUNDS, q41_c[0],
             q41_c[0] = pqc_q4171_montgomery_reduce((int64_t)q41_a[i & 511] * q41_b[i & 511]));
    RUN_CASE("q4171_barrett", "q=4171777 Barrett reduce", ROUNDS, q41_c[0],
             q41_a[0] = (int32_t)i, q41_c[0] = pqc_q4171_barrett_reduce(q41_a[0]));
    RUN_CASE("q4171_ntt_512", "512-point NTT over q=4171777", 5000, q41_c[0],
             memcpy(q41_c, q41_a, sizeof(q41_c)), pqc_q4171_ntt(q41_c));
    RUN_CASE("q4171_invntt_512", "512-point inverse NTT over q=4171777", 5000, q41_c[0],
             memcpy(q41_c, q41_a, sizeof(q41_c)), pqc_q4171_invntt_frominvmont(q41_c));
    RUN_CASE("q4171_pointwise_512", "512-point Montgomery product over q=4171777", ROUNDS, q41_c[0],
             pqc_q4171_pointwise_invmontgomery(q41_c, q41_a, q41_b, PQC_Q4171_N));
    RUN_CASE("q4171_chknorm_512", "512-point norm check over q=4171777", ROUNDS, q41_c[0],
             q41_c[0] = (int32_t)pqc_q4171_chknorm(q41_a, PQC_Q4171_N, 1000));

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
