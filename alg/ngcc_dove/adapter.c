/* Adapter from the ICCS NGCC signature ABI to the library signature API.
 *
 * Reference sources: ICCS NGCC round-1 submission package, "DOVE" (digital
 * signature track), Reference_Implementation/DOVE_pkc_skc_ref. Only the
 * intra-package include paths, the GF(256) vector layer and one out-of-bounds
 * read were changed; see reference/gf256.h and reference/SIG_AlgorithmInstance.c.
 *
 * ICCS NGCC public review, sign-09: sign-09-3 (out-of-bounds reads) is fixed in
 * the reference sources, and sign-09-1 (unauthenticated salt) can be fixed at
 * build time with PQC_NGCC_DOVE_BIND_SALT. sign-09-2 and sign-09-4 are gaps in
 * the submitted security proof; they have no code-level fix. */

#include "SIG_AlgorithmInstance.h"
#include "core/pqc_signature.h"
#include "ngcc/drng.h"
#include "util/rand.h"
#include <string.h>

/* Referenced by SIG_AlgorithmInstance.c; the KAT driver seeds it directly. */
DRNG_ctx drng_algorithm;

#if DOVE_PARAM_SET == DOVE128
#define PQC_DOVE_NAME "DOVE-128"
#elif DOVE_PARAM_SET == DOVE256
#define PQC_DOVE_NAME "DOVE-256"
#else
#define PQC_DOVE_NAME "DOVE-512"
#endif

static int keypair(uint8_t *pk, uint8_t *sk)
{
    uint8_t seed[32];
    unsigned long long pk_len = sig_get_pk_len_bytes();
    unsigned long long sk_len = sig_get_sk_len_bytes();
    if (pqc_randombytes(seed, sizeof(seed), NULL, NULL) != 0)
        return -1;
    if (init_random_number(&drng_algorithm, seed, sizeof(seed)) != 0)
        return -1;
    return sig_keygen(pk, &pk_len, sk, &sk_len);
}

static int sign_adapter(uint8_t *out, size_t *len, const uint8_t *m, size_t n,
                        const uint8_t *sk)
{
    unsigned long long sn_len = sig_get_sn_len_bytes();
    int r = sig_sign((unsigned char *)sk, sig_get_sk_len_bytes(),
                     (unsigned char *)m, n, out, &sn_len);
    if (r != 0)
        return -1;
    memmove(out + sn_len, m, n);
    *len = (size_t)sn_len + n;
    return 0;
}

static int open_adapter(uint8_t *out, size_t *len, const uint8_t *sm, size_t n,
                        const uint8_t *pk)
{
    const size_t sn_len = (size_t)sig_get_sn_len_bytes();
    if (n < sn_len)
        return -1;
    int r = sig_verify((unsigned char *)pk, sig_get_pk_len_bytes(),
                       (unsigned char *)sm, sn_len,
                       (unsigned char *)(sm + sn_len), n - sn_len);
    if (r != 0)
        return -1;
    *len = n - sn_len;
    memmove(out, sm + sn_len, *len);
    return 0;
}

static const pqc_signature algorithm = {
    PQC_DOVE_NAME,
    DOVE_PK_SEED_BYTES + DOVE_M * UPPER_SIZE(DOVE_O),
    DOVE_SK_SEED_BYTES,
    DOVE_N + DOVE_SALT_BYTES,
    keypair,
    sign_adapter,
    open_adapter,
};

const pqc_signature *pqc_ngcc_dove_algorithm(void)
{
    return &algorithm;
}
