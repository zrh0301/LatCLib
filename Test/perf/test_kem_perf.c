/* Unified KEM benchmark: one case per argv[1] so each operation can be
 * registered as its own performance test. */

#include "test.h"
#include "core/pqc_kem.h"

#ifndef PQC_KEM_GETTER
#error "Select the algorithm descriptor"
#endif

int main(int argc, char **argv)
{
    const pqc_kem *a = PQC_KEM_GETTER();
    uint8_t *pk = malloc(a->public_key_bytes);
    uint8_t *sk = malloc(a->secret_key_bytes);
    uint8_t *ct = malloc(a->ciphertext_bytes);
    uint8_t *ss = malloc(a->shared_secret_bytes);
    uint8_t *ss2 = malloc(a->shared_secret_bytes);
    char label[96];
    int rc = 1;

    if (argc != 2) {
        fprintf(stderr, "Pass one of: keypair, encaps, decaps\n");
        return 2;
    }
    if (pk == NULL || sk == NULL || ct == NULL || ss == NULL || ss2 == NULL) {
        fprintf(stderr, "allocation failed\n");
        goto done;
    }
    if (a->keypair(pk, sk) != 0) {
        fprintf(stderr, "keypair failed\n");
        goto done;
    }

    if (strcmp(argv[1], "keypair") == 0) {
        const uint64_t rounds = 10;
        test_perf_mark start = test_perf_begin();
        for (uint64_t i = 0; i < rounds; ++i) {
            if (a->keypair(pk, sk) != 0) {
                fprintf(stderr, "keypair failed\n");
                goto done;
            }
        }
        snprintf(label, sizeof label, "%s keypair", a->name);
        test_perf_print(label, test_perf_end(start), rounds);
    } else if (strcmp(argv[1], "encaps") == 0) {
        const uint64_t rounds = 20;
        test_perf_mark start = test_perf_begin();
        for (uint64_t i = 0; i < rounds; ++i) {
            if (a->encaps(ct, ss, pk) != 0) {
                fprintf(stderr, "encaps failed\n");
                goto done;
            }
        }
        snprintf(label, sizeof label, "%s encaps", a->name);
        test_perf_print(label, test_perf_end(start), rounds);
    } else if (strcmp(argv[1], "decaps") == 0) {
        const uint64_t rounds = 20;
        if (a->encaps(ct, ss, pk) != 0) {
            fprintf(stderr, "encaps failed\n");
            goto done;
        }
        test_perf_mark start = test_perf_begin();
        for (uint64_t i = 0; i < rounds; ++i) {
            if (a->decaps(ss2, ct, sk) != 0) {
                fprintf(stderr, "decaps failed\n");
                goto done;
            }
        }
        snprintf(label, sizeof label, "%s decaps", a->name);
        test_perf_print(label, test_perf_end(start), rounds);
    } else {
        fprintf(stderr, "unknown case: %s\n", argv[1]);
        rc = 2;
        goto done;
    }
    rc = 0;
done:
    if (sk != NULL) {
        memset(sk, 0, a->secret_key_bytes);
    }
    free(pk);
    free(sk);
    free(ct);
    free(ss);
    free(ss2);
    return rc;
}
