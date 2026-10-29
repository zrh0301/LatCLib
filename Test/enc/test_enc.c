#include "test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ENC_PERF_ITERATIONS 10000


static void fill_bytes(uint8_t *buf, size_t len, uint8_t seed) {
    for(size_t i = 0; i < len; i++) {
        buf[i] = (uint8_t)(seed + 17u * (uint8_t)i + (uint8_t)(i >> 3));
    }
}

static int check_decrypt(const test_enc_alg *alg, const char *label, const uint8_t *ct, const uint8_t *sk, const uint8_t *msg) {
    uint8_t *out = calloc(1, alg->msg_bytes);
    int ok = 1;
    if(!out) {
        printf("  [FAIL] allocation failed for %s %s decrypt\n", alg->name, label);
        return 0;
    }
    if(alg->decrypt(out, sk, ct, alg->param) != 0) {
        printf("  [FAIL] %s decrypt failed after %s\n", alg->name, label);
        ok = 0;
    }
    if(memcmp(msg, out, alg->msg_bytes) != 0) {
        printf("  [FAIL] %s decrypt mismatch after %s\n", alg->name, label);
        ok = 0;
    }
    free(out);
    return ok;
}

static int test_encrypt_repeatable(const test_enc_alg *alg, const uint8_t *pk, const uint8_t *msg, const uint8_t *coins) {
    uint8_t *ct1 = calloc(1, alg->ct_bytes);
    uint8_t *ct2 = calloc(1, alg->ct_bytes);
    int ok = 1;
    if(!ct1 || !ct2) {
        ok = 0;
        goto cleanup;
    }
    if(alg->encrypt_with_coins(ct1, pk, msg, coins, alg->param) != 0 || alg->encrypt_with_coins(ct2, pk, msg, coins, alg->param) != 0 || memcmp(ct1, ct2, alg->ct_bytes) != 0) {
        printf("  [FAIL] %s encrypt_with_coins is not repeatable\n", alg->name);
        ok = 0;
    }
cleanup:
    free(ct1);
    free(ct2);
    return ok;
}

static int test_one_enc_alg(const test_enc_alg *alg, uint8_t seed_tag) {
    uint8_t *pk = calloc(1, alg->pk_bytes);
    uint8_t *sk = calloc(1, alg->sk_bytes);
    uint8_t *ct = calloc(1, alg->ct_bytes);
    uint8_t *seed = calloc(1, alg->seed_bytes);
    uint8_t *coins = calloc(1, alg->coins_bytes);
    uint8_t *msg = calloc(1, alg->msg_bytes);
    int ok = 1;

    if(!pk || !sk || !ct || !seed || !coins || !msg) {
        printf("  [FAIL] allocation failed for %s\n", alg->name);
        ok = 0;
        goto cleanup;
    }
    if(!alg->keygen_with_seed || !alg->encrypt_with_coins || !alg->decrypt) {
        printf("  [FAIL] %s has incomplete encryption functions\n", alg->name);
        ok = 0;
        goto cleanup;
    }

    fill_bytes(seed, alg->seed_bytes, (uint8_t)(0x20 + seed_tag));
    fill_bytes(coins, alg->coins_bytes, (uint8_t)(0x40 + seed_tag));
    fill_bytes(msg, alg->msg_bytes, (uint8_t)(0x60 + seed_tag));

    if(alg->keygen_with_seed(pk, sk, seed, alg->param) != 0) {
        printf("  [FAIL] %s keygen_with_seed failed\n", alg->name);
        ok = 0;
        goto cleanup;
    }
    if(alg->check_pk && !alg->check_pk(pk, alg->param)) {
        printf("  [FAIL] %s public key check failed\n", alg->name);
        ok = 0;
    }
    if(alg->encrypt_with_coins(ct, pk, msg, coins, alg->param) != 0) {
        printf("  [FAIL] %s encrypt_with_coins failed\n", alg->name);
        ok = 0;
    } else {
        ok &= check_decrypt(alg, "encrypt_with_coins", ct, sk, msg);
        ok &= test_encrypt_repeatable(alg, pk, msg, coins);
    }

cleanup:
    free(pk);
    free(sk);
    free(ct);
    free(seed);
    free(coins);
    free(msg);
    return ok;
}

static int test_enc_correctness(void) {
    int ok = 1;
    printf("\n============================================\n");
    printf("  Test: Encryption Correctness\n");
    printf("============================================\n");
    for(size_t i = 0; i < test_enc_alg_count; i++) {
        int one_ok = test_one_enc_alg(test_enc_algs[i], (uint8_t)i);
        printf("  %-24s %s\n", test_enc_algs[i]->name, one_ok ? "[PASS]" : "[FAIL]");
        ok &= one_ok;
    }
    return ok;
}

static void print_perf(const test_enc_alg *alg, const char *label, test_perf_result result) {
    char module[128];
    snprintf(module, sizeof(module), "%s %s", alg->name, label);
    test_perf_print(module, result, ENC_PERF_ITERATIONS);
}

static void bench_one_enc_alg(const test_enc_alg *alg, uint8_t seed_tag) {
    uint8_t *pk = calloc(1, alg->pk_bytes);
    uint8_t *sk = calloc(1, alg->sk_bytes);
    uint8_t *ct = calloc(1, alg->ct_bytes);
    uint8_t *out = calloc(1, alg->msg_bytes);
    uint8_t *seed = calloc(1, alg->seed_bytes);
    uint8_t *coins = calloc(1, alg->coins_bytes);
    uint8_t *msg = calloc(1, alg->msg_bytes);
    test_perf_mark timer;

    if(!pk || !sk || !ct || !out || !seed || !coins || !msg) {
        printf("  [FAIL] allocation failed for %s perf\n", alg->name);
        goto cleanup;
    }

    fill_bytes(seed, alg->seed_bytes, (uint8_t)(0x10 + seed_tag));
    fill_bytes(coins, alg->coins_bytes, (uint8_t)(0x30 + seed_tag));
    fill_bytes(msg, alg->msg_bytes, (uint8_t)(0x50 + seed_tag));
    alg->keygen_with_seed(pk, sk, seed, alg->param);

    timer = test_perf_begin();
    for(int i = 0; i < ENC_PERF_ITERATIONS; i++) {
        alg->keygen_with_seed(pk, sk, seed, alg->param);
    }
    print_perf(alg, "keygen_with_seed", test_perf_end(timer));

    timer = test_perf_begin();
    for(int i = 0; i < ENC_PERF_ITERATIONS; i++) {
        alg->encrypt_with_coins(ct, pk, msg, coins, alg->param);
    }
    print_perf(alg, "encrypt_with_coins", test_perf_end(timer));

    timer = test_perf_begin();
    for(int i = 0; i < ENC_PERF_ITERATIONS; i++) {
        alg->decrypt(out, sk, ct, alg->param);
    }
    print_perf(alg, "decrypt", test_perf_end(timer));

cleanup:
    free(pk);
    free(sk);
    free(ct);
    free(out);
    free(seed);
    free(coins);
    free(msg);
}

static void test_enc_perf(void) {
    for(size_t i = 0; i < test_enc_alg_count; i++) {
        bench_one_enc_alg(test_enc_algs[i], (uint8_t)i);
    }
}

int main(int argc, char **argv) {
    int ok;
#ifdef PQC_PERF_ONLY
    (void)argc; (void)argv;
    test_enc_perf();
    return 0;
#else
    test_init();
    printf("============================================\n");
    printf("  Test Suite: Encryption Algorithms\n");
    printf("============================================\n");
    ok = test_enc_correctness();
    if (argc == 2 && strcmp(argv[1], "--bench") == 0)
        test_enc_perf();
#endif
    return ok ? 0 : 1;
}
