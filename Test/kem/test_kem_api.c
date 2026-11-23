/* Exercises the shared KEM interface against one algorithm descriptor.
 * The caller selects the algorithm with PQC_KEM_GETTER. Besides the round
 * trip, the test requires implicit rejection to actually change the returned
 * secret when the ciphertext is modified (kem-01-1). */

#include "core/pqc_kem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PQC_KEM_GETTER
#error "PQC_KEM_GETTER must name the algorithm accessor"
#endif

static int all_equal(const uint8_t *a, const uint8_t *b, size_t len)
{
    return memcmp(a, b, len) == 0;
}

int main(void)
{
    const pqc_kem *kem = PQC_KEM_GETTER();
    uint8_t *pk, *sk, *ct, *ct_bad, *ss, *ss1, *ss2;
    int failures = 0;

    if (kem == NULL || kem->keypair == NULL || kem->encaps == NULL ||
        kem->decaps == NULL || kem->public_key_bytes == 0 ||
        kem->secret_key_bytes == 0 || kem->ciphertext_bytes == 0 ||
        kem->shared_secret_bytes == 0) {
        fprintf(stderr, "invalid KEM descriptor\n");
        return 1;
    }

    pk = malloc(kem->public_key_bytes);
    sk = malloc(kem->secret_key_bytes);
    ct = malloc(kem->ciphertext_bytes);
    ct_bad = malloc(kem->ciphertext_bytes);
    ss = malloc(kem->shared_secret_bytes);
    ss1 = malloc(kem->shared_secret_bytes);
    ss2 = malloc(kem->shared_secret_bytes);
    if (pk == NULL || sk == NULL || ct == NULL || ct_bad == NULL ||
        ss == NULL || ss1 == NULL || ss2 == NULL) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    if (kem->keypair(pk, sk) != 0) {
        fprintf(stderr, "keypair failed\n");
        return 1;
    }
    if (kem->encaps(ct, ss, pk) != 0) {
        fprintf(stderr, "encaps failed\n");
        failures++;
    }
    if (kem->decaps(ss1, ct, sk) != 0 || !all_equal(ss, ss1, kem->shared_secret_bytes)) {
        fprintf(stderr, "decapsulation does not recover the encapsulated secret\n");
        failures++;
    }

    /* Every ciphertext bit must be bound into the rejection value. */
    memcpy(ct_bad, ct, kem->ciphertext_bytes);
    ct_bad[kem->ciphertext_bytes / 2] ^= 0x01;
    if (kem->decaps(ss1, ct_bad, sk) != 0) {
        fprintf(stderr, "decapsulation of a modified ciphertext failed\n");
        failures++;
    } else if (all_equal(ss, ss1, kem->shared_secret_bytes)) {
        fprintf(stderr, "modified ciphertext returned the real shared secret\n");
        failures++;
    }

    /* Implicit rejection must be deterministic for the same ciphertext. */
    if (kem->decaps(ss2, ct_bad, sk) != 0 ||
        !all_equal(ss1, ss2, kem->shared_secret_bytes)) {
        fprintf(stderr, "implicit rejection is not deterministic\n");
        failures++;
    }

    if (failures == 0)
        printf("%s KEM interface ok\n", kem->name);

    free(pk);
    free(sk);
    free(ct);
    free(ct_bad);
    free(ss);
    free(ss1);
    free(ss2);
    return failures == 0 ? 0 : 1;
}
