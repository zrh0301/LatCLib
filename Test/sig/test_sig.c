#include "test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SIG_PERF_ITERATIONS 10000

typedef struct {
    uint8_t seed;
    size_t used;
} test_rng_state;

typedef struct {
    const uint8_t *data;
    size_t size;
    size_t used;
} test_rng_pool;


static size_t sig_len(const test_sig_alg *alg) {
    return alg->sig_bytes;
}

static void fill_bytes(uint8_t *buf, size_t len, uint8_t seed) {
    for(size_t i = 0; i < len; i++) {
        buf[i] = (uint8_t)(seed + 17u * (uint8_t)i + (uint8_t)(i >> 3));
    }
}

static int test_rng(void *rng_ctx, uint8_t *out, size_t out_len) {
    test_rng_state *s = (test_rng_state *)rng_ctx;
    for(size_t i = 0; i < out_len; i++) {
        out[i] = (uint8_t)(s->seed + 31u * (uint8_t)(s->used + i));
    }
    s->used += out_len;
    return 0;
}

static int test_pool_rng(void *rng_ctx, uint8_t *out, size_t out_len) {
    test_rng_pool *pool = (test_rng_pool *)rng_ctx;
    if(pool->used > pool->size || out_len > pool->size - pool->used) {
        return -1;
    }
    memcpy(out, pool->data + pool->used, out_len);
    pool->used += out_len;
    return 0;
}

static int sig_nonzero(const uint8_t *sig, size_t len) {
    for(size_t i = 0; i < len; i++) {
        if(sig[i] != 0) {
            return 1;
        }
    }
    return 0;
}

static int pk_nonzero(const dili_pk *pk) {
    for(int i = 0; i < 32; i++) {
        if(pk->rho[i] != 0) {
            return 1;
        }
    }
    return 0;
}

static int check_signature(const test_sig_alg *alg, const char *label, const uint8_t *sig, size_t slen, const uint8_t *msg, size_t mlen, const dili_pk *pk) {
    uint8_t bad_msg[] = "hello, signature algorithm.";
    int ok = 1;
    if(!sig_nonzero(sig, slen)) {
        printf("  [FAIL] %s %s produced all-zero signature\n", alg->name, label);
        ok = 0;
    }
    if(alg->verify(sig, msg, (int)mlen, pk, alg->param) != 0) {
        printf("  [FAIL] %s verify rejected %s signature\n", alg->name, label);
        ok = 0;
    }
    if(alg->verify(sig, bad_msg, (int)strlen((const char *)bad_msg), pk, alg->param) == 0) {
        printf("  [FAIL] %s verify accepted wrong message for %s\n", alg->name, label);
        ok = 0;
    }
    return ok;
}

static int test_one_sig_alg(const test_sig_alg *alg) {
    uint8_t msg[] = "Hello, signature algorithm.";
    size_t mlen = strlen((const char *)msg);
    size_t slen = sig_len(alg);
    dili_keypair_t *keypair = calloc(1, sizeof(*keypair));
    uint8_t *sig = calloc(1, slen);
    uint8_t *sig2 = calloc(1, slen);
    uint8_t *rnd = alg->rnd_bytes ? calloc(1, alg->rnd_bytes) : NULL;
    const dili_pk *pk;
    const dili_sk *sk;
    test_rng_state rng = {0x43, 0};
    int ok = 1;
    int tested = 0;

    if(!keypair || !sig || !sig2 || (alg->rnd_bytes && !rnd)) {
        printf("  [FAIL] allocation failed for %s\n", alg->name);
        ok = 0;
        goto cleanup;
    }
    alg->keygen(keypair, alg->param);
    pk = &keypair->pk;
    sk = &keypair->sk;
    if(!pk_nonzero(pk)) {
        printf("  [FAIL] %s public key looks empty\n", alg->name);
        ok = 0;
    }

    if(alg->sign) {
        tested++;
        memset(sig, 0, slen);
        if(alg->sign(sig, msg, (int)mlen, sk, alg->param) != 0) {
            printf("  [FAIL] %s sign failed\n", alg->name);
            ok = 0;
        } else {
            ok &= check_signature(alg, "sign", sig, slen, msg, mlen, pk);
        }
    }
    if(alg->sign_with_rnd && alg->rnd_bytes) {
        tested++;
        fill_bytes(rnd, alg->rnd_bytes, 0xA5);
        memset(sig, 0, slen);
        memset(sig2, 0, slen);
        if(alg->sign_with_rnd(sig, msg, (int)mlen, sk, alg->param, rnd) != 0 || alg->sign_with_rnd(sig2, msg, (int)mlen, sk, alg->param, rnd) != 0) {
            printf("  [FAIL] %s sign_with_rnd failed\n", alg->name);
            ok = 0;
        } else {
            ok &= check_signature(alg, "sign_with_rnd", sig, slen, msg, mlen, pk);
            if(memcmp(sig, sig2, slen) != 0) {
                printf("  [FAIL] %s sign_with_rnd is not repeatable\n", alg->name);
                ok = 0;
            }
        }
    }
    if(alg->sign_with_rng) {
        tested++;
        memset(sig, 0, slen);
        rng.used = 0;
        if(alg->sign_with_rng(sig, msg, (int)mlen, sk, alg->param, test_rng, &rng) != 0) {
            printf("  [FAIL] %s sign_with_rng failed\n", alg->name);
            ok = 0;
        } else {
            ok &= check_signature(alg, "sign_with_rng", sig, slen, msg, mlen, pk);
        }
    }
    if(tested == 0) {
        printf("  [FAIL] %s has no sign function\n", alg->name);
        ok = 0;
    }

cleanup:
    free(keypair);
    free(sig);
    free(sig2);
    free(rnd);
    return ok;
}

static int test_sig_correctness(void) {
    int ok = 1;
    printf("\n============================================\n");
    printf("  Test: Signature Correctness\n");
    printf("============================================\n");
    for(size_t i = 0; i < test_sig_alg_count; i++) {
        int one_ok = test_one_sig_alg(test_sig_algs[i]);
        printf("  %-24s %s\n", test_sig_algs[i]->name, one_ok ? "[PASS]" : "[FAIL]");
        ok &= one_ok;
    }
    return ok;
}

static void print_perf(const test_sig_alg *alg, const char *label, test_perf_result result) {
    char module[128];
    snprintf(module, sizeof(module), "%s %s", alg->name, label);
    test_perf_print(module, result, SIG_PERF_ITERATIONS);
}

static void bench_one_sig_alg(const test_sig_alg *alg) {
    uint8_t msg[] = "Hello, signature performance.";
    size_t mlen = strlen((const char *)msg);
    size_t slen = sig_len(alg);
    dili_keypair_t *keypair = calloc(1, sizeof(*keypair));
    uint8_t *sig = calloc(1, slen);
    size_t rnd_pool_size = alg->rnd_bytes * SIG_PERF_ITERATIONS;
    uint8_t *rnd_pool = rnd_pool_size ? calloc(1, rnd_pool_size) : NULL;
    const dili_pk *pk;
    const dili_sk *sk;
    test_rng_pool rng_pool = {rnd_pool, rnd_pool_size, 0};
    test_perf_mark timer;

    if(!keypair || !sig || (rnd_pool_size && !rnd_pool)) {
        printf("  [FAIL] allocation failed for %s perf\n", alg->name);
        goto cleanup;
    }
    if(rnd_pool) {
        test_random_bytes(rnd_pool, rnd_pool_size);
    }
    alg->keygen(keypair, alg->param);
    pk = &keypair->pk;
    sk = &keypair->sk;

    timer = test_perf_begin();
    for(int i = 0; i < SIG_PERF_ITERATIONS; i++) {
        alg->keygen(keypair, alg->param);
    }
    print_perf(alg, "KeyGen", test_perf_end(timer));

    pk = &keypair->pk;
    sk = &keypair->sk;
    if(alg->sign) {
        timer = test_perf_begin();
        for(int i = 0; i < SIG_PERF_ITERATIONS; i++) {
            alg->sign(sig, msg, (int)mlen, sk, alg->param);
        }
        print_perf(alg, "sign", test_perf_end(timer));
    }
    if(alg->sign_with_rnd && rnd_pool) {
        timer = test_perf_begin();
        for(int i = 0; i < SIG_PERF_ITERATIONS; i++) {
            alg->sign_with_rnd(sig, msg, (int)mlen, sk, alg->param, rnd_pool + (size_t)i * alg->rnd_bytes);
        }
        print_perf(alg, "sign_with_rnd", test_perf_end(timer));
    }
    if(alg->sign_with_rng && rnd_pool) {
        timer = test_perf_begin();
        rng_pool.used = 0;
        for(int i = 0; i < SIG_PERF_ITERATIONS; i++) {
            alg->sign_with_rng(sig, msg, (int)mlen, sk, alg->param, test_pool_rng, &rng_pool);
        }
        print_perf(alg, "sign_with_rng", test_perf_end(timer));
    }
    if(alg->sign) {
        alg->sign(sig, msg, (int)mlen, sk, alg->param);
    } else if(alg->sign_with_rnd && rnd_pool) {
        alg->sign_with_rnd(sig, msg, (int)mlen, sk, alg->param, rnd_pool);
    } else if(alg->sign_with_rng && rnd_pool) {
        rng_pool.used = 0;
        alg->sign_with_rng(sig, msg, (int)mlen, sk, alg->param, test_pool_rng, &rng_pool);
    }

    timer = test_perf_begin();
    for(int i = 0; i < SIG_PERF_ITERATIONS; i++) {
        alg->verify(sig, msg, (int)mlen, pk, alg->param);
    }
    print_perf(alg, "Verify", test_perf_end(timer));

cleanup:
    free(keypair);
    free(sig);
    free(rnd_pool);
}

static void test_sig_perf(void) {
    for(size_t i = 0; i < test_sig_alg_count; i++) {
        bench_one_sig_alg(test_sig_algs[i]);
    }
}

int main(int argc, char **argv) {
    int ok;
#ifdef PQC_PERF_ONLY
    (void)argc; (void)argv;
    test_sig_perf();
    return 0;
#else
    test_init();
    printf("============================================\n");
    printf("  Test Suite: Signature Algorithms\n");
    printf("============================================\n");
    ok = test_sig_correctness();
    if (argc == 2 && strcmp(argv[1], "--bench") == 0)
        test_sig_perf();
#endif
    return ok ? 0 : 1;
}
