/* Unified signature benchmark: one case per argv[1] so each operation can be
 * registered as its own performance test. */

#include "test.h"
#include "core/pqc_signature.h"

#ifndef PQC_SIG_GETTER
#error "Select the algorithm descriptor"
#endif

int main(int argc, char **argv)
{
    const pqc_signature *a = PQC_SIG_GETTER();
    uint8_t msg[32];
    size_t cap = a->signature_bytes + sizeof(msg);
    uint8_t *pk = malloc(a->public_key_bytes);
    uint8_t *sk = malloc(a->secret_key_bytes);
    uint8_t *sm = malloc(cap);
    char label[96];
    uint8_t *out = malloc(cap);
    size_t slen = 0, olen = 0;
    int rc = 1;

    if (argc != 2) {
        fprintf(stderr, "Pass one of: keygen, sign, open\n");
        return 2;
    }
    if (pk == NULL || sk == NULL || sm == NULL || out == NULL) {
        fprintf(stderr, "allocation failed\n");
        goto done;
    }
    for (size_t i = 0; i < sizeof(msg); ++i) {
        msg[i] = (uint8_t)i;
    }
    if (a->keypair(pk, sk) != 0) {
        fprintf(stderr, "keypair failed\n");
        goto done;
    }

    if (strcmp(argv[1], "keygen") == 0) {
        const uint64_t rounds = 5;
        test_perf_mark start = test_perf_begin();
        for (uint64_t i = 0; i < rounds; ++i) {
            if (a->keypair(pk, sk) != 0) {
                fprintf(stderr, "keypair failed\n");
                goto done;
            }
        }
        snprintf(label, sizeof label, "%s keypair", a->name);
        test_perf_print(label, test_perf_end(start), rounds);
    } else if (strcmp(argv[1], "sign") == 0) {
        const uint64_t rounds = 10;
        test_perf_mark start = test_perf_begin();
        for (uint64_t i = 0; i < rounds; ++i) {
            if (a->sign(sm, &slen, msg, sizeof(msg), sk) != 0) {
                fprintf(stderr, "sign failed\n");
                goto done;
            }
        }
        snprintf(label, sizeof label, "%s sign", a->name);
        test_perf_print(label, test_perf_end(start), rounds);
    } else if (strcmp(argv[1], "open") == 0) {
        const uint64_t rounds = 10;
        if (a->sign(sm, &slen, msg, sizeof(msg), sk) != 0) {
            fprintf(stderr, "sign failed\n");
            goto done;
        }
        test_perf_mark start = test_perf_begin();
        for (uint64_t i = 0; i < rounds; ++i) {
            if (a->open(out, &olen, sm, slen, pk) != 0) {
                fprintf(stderr, "open failed\n");
                goto done;
            }
        }
        snprintf(label, sizeof label, "%s open", a->name);
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
    free(sm);
    free(out);
    return rc;
}
