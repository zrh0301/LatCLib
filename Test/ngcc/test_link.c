#include "core/pqc_kem.h"
#include "core/pqc_signature.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int signature_roundtrip(const pqc_signature *s)
{
    const uint8_t msg[] = {1, 2, 3, 4};
    uint8_t *pk = malloc(s->public_key_bytes);
    uint8_t *sk = malloc(s->secret_key_bytes);
    uint8_t *sm = malloc(s->signature_bytes + sizeof(msg));
    uint8_t *out = malloc(s->signature_bytes + sizeof(msg));
    size_t slen = 0, mlen = 0;
    int failed = !pk || !sk || !sm || !out ||
                 s->keypair(pk, sk) != 0 ||
                 s->sign(sm, &slen, msg, sizeof(msg), sk) != 0 ||
                 s->open(out, &mlen, sm, slen, pk) != 0 ||
                 mlen != sizeof(msg) || memcmp(out, msg, sizeof(msg)) != 0;
    free(out);
    free(sm);
    free(sk);
    free(pk);
    if (failed)
        fprintf(stderr, "Coexistence: %s signature round trip failed\n", s->name);
    return failed;
}

static int kem_roundtrip(const pqc_kem *k)
{
    uint8_t *pk = malloc(k->public_key_bytes);
    uint8_t *sk = malloc(k->secret_key_bytes);
    uint8_t *ct = malloc(k->ciphertext_bytes);
    uint8_t *ss = malloc(k->shared_secret_bytes);
    uint8_t *out = malloc(k->shared_secret_bytes);
    int failed = !pk || !sk || !ct || !ss || !out ||
                 k->keypair(pk, sk) != 0 ||
                 k->encaps(ct, ss, pk) != 0 ||
                 k->decaps(out, ct, sk) != 0 ||
                 memcmp(out, ss, k->shared_secret_bytes) != 0;
    free(out);
    free(ss);
    free(ct);
    free(sk);
    free(pk);
    if (failed)
        fprintf(stderr, "Coexistence: %s KEM round trip failed\n", k->name);
    return failed;
}

int main(void)
{
    const pqc_signature *algorithms[] = {
        pqc_ngcc_dove_algorithm(),
        pqc_ngcc_vdoo_algorithm(),
        pqc_ngcc_flextree_algorithm(),
        pqc_ngcc_aigis_sig_algorithm(),
        pqc_ngcc_tsuov_algorithm(),
#ifdef PQC_NGCC_LINK_TEST_SNOVA
        pqc_snova_r3_algorithm(),
#endif
    };
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(*algorithms); ++i)
        if (signature_roundtrip(algorithms[i]))
            return 1;
    return kem_roundtrip(pqc_ngcc_aigis_algorithm()) ||
           kem_roundtrip(pqc_ngcc_scloud_algorithm());
}
