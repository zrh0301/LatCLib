/* Adapter from the ICCS NGCC signature ABI to the library signature API.
 *
 * Reference sources: ICCS NGCC round-1 submission package, "FlexTree" (digital
 * signature track), Reference_Implementation/Flextree-160f. The submission
 * ships no single parameter-independent build: every parameter set carries its
 * own params header and a multi-megabyte precomputed PORS+FP table, so this
 * module exposes the submitted default set, Flextree-160f, only.
 *
 * Local changes to the reference layer:
 *   - hash.c propagates allocation failures instead of signing with an
 *     unwritten randomizer and digest (sign-11-9);
 *   - pors_fp.c rejects non-zero padding in the unused authentication-node
 *     suffix (sign-11-5);
 *   - wots.c checks the specification's digit sum instead of its complement
 *     (sign-11-8, behaviour identical on the shipped set);
 *   - PQC_NGCC_FLEXTREE_STRICT additionally replaces the overlapping
 *     WOTS_ZERO_BITS window with the exact trailing-bit check (sign-11-7).
 * The default build keeps the submitted signer predicate so the published KAT
 * remains reproducible.
 *
 * The remaining public-review findings are protocol-level: sign-11-1 (PORS
 * counter grinding, design), sign-11-3 (the normative verifier omits the OTS
 * encoding checks, while the submitted one performs them) and sign-11-6
 * (specification-only leaf-address reuse). sign-11-2 and sign-11-4 were
 * withdrawn.
 *
 * The reference DRNG is a file-scope context that the KAT driver seeds
 * directly. The adapter seeds it from the system RNG before key generation and
 * before each signature, so library callers never share the KAT stream. */

#include "SIG_AlgorithmInstance.h"
#include "core/pqc_signature.h"
#include "util/rand.h"
#include "ngcc/drng.h"
#include "api.h"
#include <string.h>

/* Referenced by sign.c. The KAT driver
 * seeds this context directly, as the submission's own driver does. */
DRNG_ctx drng_algorithm;

#define PQC_FLEXTREE_NAME "Flextree-160f"

static int seed_drng(void)
{
    uint8_t seed[32];
    if (pqc_randombytes(seed, sizeof(seed), NULL, NULL) != 0)
        return -1;
    if (init_random_number(&drng_algorithm, seed, sizeof(seed)) != 0)
        return -1;
    return 0;
}

static int keypair(uint8_t *pk, uint8_t *sk)
{
    unsigned long long pk_len = CRYPTO_PUBLICKEYBYTES;
    unsigned long long sk_len = CRYPTO_SECRETKEYBYTES;
    if (seed_drng() != 0)
        return -1;
    return sig_keygen(pk, &pk_len, sk, &sk_len);
}

static int sign_adapter(uint8_t *out, size_t *len, const uint8_t *m, size_t n,
                        const uint8_t *sk)
{
    unsigned long long sn_len = CRYPTO_BYTES;
    if (seed_drng() != 0)
        return -1;
    if (sig_sign((unsigned char *)sk, CRYPTO_SECRETKEYBYTES,
                 (unsigned char *)m, n, out, &sn_len) != 0)
        return -1;
    memmove(out + sn_len, m, n);
    *len = (size_t)sn_len + n;
    return 0;
}

static int open_adapter(uint8_t *out, size_t *len, const uint8_t *sm, size_t n,
                        const uint8_t *pk)
{
    if (n < CRYPTO_BYTES)
        return -1;
    if (sig_verify((unsigned char *)pk, CRYPTO_PUBLICKEYBYTES,
                   (unsigned char *)sm, CRYPTO_BYTES,
                   (unsigned char *)(sm + CRYPTO_BYTES),
                   n - CRYPTO_BYTES) != 0)
        return -1;
    *len = n - CRYPTO_BYTES;
    memmove(out, sm + CRYPTO_BYTES, *len);
    return 0;
}

static const pqc_signature algorithm = {
    PQC_FLEXTREE_NAME,
    CRYPTO_PUBLICKEYBYTES,
    CRYPTO_SECRETKEYBYTES,
    CRYPTO_BYTES,
    keypair,
    sign_adapter,
    open_adapter,
};

const pqc_signature *pqc_ngcc_flextree_algorithm(void)
{
    return &algorithm;
}

/* ICCS detached-signature ABI. */
unsigned long long sig_get_pk_len_bytes()
{
    return crypto_sign_publickeybytes();
}

unsigned long long sig_get_sk_len_bytes()
{
    return crypto_sign_secretkeybytes();
}

unsigned long long sig_get_sn_len_bytes()
{
    return crypto_sign_bytes();
}

int sig_keygen(
    unsigned char *pk, unsigned long long *pk_len_bytes,
    unsigned char *sk, unsigned long long *sk_len_bytes)
{
    int ret = crypto_sign_keypair(pk, sk);
    
    if (ret == 0) {
        if (pk_len_bytes) *pk_len_bytes = crypto_sign_publickeybytes();
        if (sk_len_bytes) *sk_len_bytes = crypto_sign_secretkeybytes();
    }
    
    return ret;
}

int sig_sign(
    unsigned char *sk, unsigned long long sk_len_bytes,
    unsigned char *m, unsigned long long m_len_bytes,
    unsigned char *sn, unsigned long long *sn_len_bytes)
{
    if (sk_len_bytes != CRYPTO_SECRETKEYBYTES)
        return -2;

    size_t siglen = 0;
    
    int ret = crypto_sign_signature(sn, &siglen, m, (size_t)m_len_bytes, sk);
    
    if (ret == 0 && sn_len_bytes != NULL) {
        *sn_len_bytes = (unsigned long long)siglen;
    }
    
    return ret;
}

int sig_verify(
    unsigned char *pk, unsigned long long pk_len_bytes,
    unsigned char *sn, unsigned long long sn_len_bytes,
    unsigned char *m, unsigned long long m_len_bytes)
{
    if (pk_len_bytes != CRYPTO_PUBLICKEYBYTES)
        return -2;

    return crypto_sign_verify(sn, (size_t)sn_len_bytes, m, (size_t)m_len_bytes, pk);
}
