#include "test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int run_hash_vecs(const test_hash_alg *alg) {
    uint8_t *digest = calloc(1, alg->digest_size);
    void *ctx = calloc(1, alg->ctx_size);
    int ok = 1;

    if(!digest || !ctx) {
        printf("  [FAIL] allocation failed for %s\n", alg->name);
        ok = 0;
        goto cleanup;
    }

    for(size_t i = 0; i < alg->vec_count; i++) {
        const test_hash_vec *vec = &alg->vecs[i];

        alg->digest(digest, vec->msg, vec->msg_len);
        if(memcmp(digest, vec->digest, alg->digest_size) != 0) {
            printf("  [FAIL] %s one-shot vector %zu mismatch\n", alg->name, i);
            ok = 0;
        }

        memset(digest, 0, alg->digest_size);
        memset(ctx, 0, alg->ctx_size);
        alg->starts(ctx);
        for(size_t j = 0; j < vec->msg_len; j++) {
            alg->update(ctx, vec->msg + j, 1);
        }
        alg->finish(ctx, digest);
        if(memcmp(digest, vec->digest, alg->digest_size) != 0) {
            printf("  [FAIL] %s streaming vector %zu mismatch\n", alg->name, i);
            ok = 0;
        }
    }

cleanup:
    free(digest);
    free(ctx);
    return ok;
}

static int run_hash_case(const test_hash_alg *alg) {
    int ok;

    printf("  %-16s", alg->name);
    ok = run_hash_vecs(alg);
    printf("%s\n", ok ? "[PASS]" : "[FAIL]");
    return ok;
}

int main(void) {
    int ok = 1;

    test_init();
    printf("============================================\n");
    printf("  Test Suite: Hash Algorithms\n");
    printf("============================================\n");

    for(size_t i = 0; i < test_hash_alg_count; i++) {
        ok &= run_hash_case(test_hash_algs[i]);
    }

    return ok ? 0 : 1;
}
