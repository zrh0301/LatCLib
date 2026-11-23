/* Adapter from the ICCS NGCC signature ABI to the library signature API.
 *
 * Reference sources: ICCS NGCC round-1 submission package, "VDOO" (digital
 * signature track), Reference_Implementation/vdoo_128. The tower arithmetic
 * and matrix operators are provided by the library GF layer.
 *
 * The submission never seeds its own DRNG: the reference uses a zero-initialised
 * context, which is what makes the published KAT reproducible. The
 * adapter therefore seeds that context from the system RNG before key
 * generation and before each signature, so callers of the library API never
 * run on the fixed all-zero stream. */

#include "SIG_AlgorithmInstance.h"
#include "core/pqc_signature.h"
#include "util/rand.h"
#include "vdoo_config.h"
#include "vdoo_keypair.h"
#include <string.h>

#if defined(VDOO_128)
#define PQC_VDOO_NAME "VDOO-128"
#elif defined(VDOO_256)
#define PQC_VDOO_NAME "VDOO-256"
#else
#define PQC_VDOO_NAME "VDOO-512"
#endif

#define PQC_VDOO_SK_BYTES (sizeof(sk_t))

DRNG_ctx pqc_ngcc_vdoo_drng;

static int seed_drng(void)
{
    uint8_t seed[32];
    if (pqc_randombytes(seed, sizeof(seed), NULL, NULL) != 0)
        return -1;
    return init_random_number(&pqc_ngcc_vdoo_drng, seed, sizeof(seed));
}

static int keypair(uint8_t *pk, uint8_t *sk)
{
    unsigned long long pk_len = sig_get_pk_len_bytes();
    unsigned long long sk_len = sig_get_sk_len_bytes();
    if (seed_drng() != 0)
        return -1;
    return sig_keygen(pk, &pk_len, sk, &sk_len);
}

static int sign_adapter(uint8_t *out, size_t *len, const uint8_t *m, size_t n,
                        const uint8_t *sk)
{
    unsigned long long sn_len = sig_get_sn_len_bytes();
    int r;
    if (seed_drng() != 0)
        return -1;
    r = sig_sign((unsigned char *)sk, sig_get_sk_len_bytes(),
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
    PQC_VDOO_NAME,
    sizeof(pk_t),
    PQC_VDOO_SK_BYTES,
    VDOO_N_BYTE + SALT_BYTES,
    keypair,
    sign_adapter,
    open_adapter,
};

const pqc_signature *pqc_ngcc_vdoo_algorithm(void)
{
    return &algorithm;
}
