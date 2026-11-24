/* Adapter from the ICCS NGCC signature ABI to the library signature API.
 *
 * Reference sources: ICCS NGCC round-1 submission package, "Aigis-Sig+"
 * (signature track), Reference_Implementation/Aigis-Sig+-I. Only the first
 * parameter set is shipped; sets II and III differ in params.h.
 *
 * ICCS NGCC public review, sign-01:
 *   - sign-01-1 (non-canonical hint padding) and sign-01-2 (attacker-controlled
 *     hint counts) are fixed in packing.c, which now bounds-checks the whole
 *     hint against the declared signature length and requires the canonical
 *     re-encoding;
 *   - sign-01-3 (ignored key lengths) is fixed in SIG_AlgorithmInstance.c;
 *   - sign-01-4 (one polynomial written past the mask vector) is fixed in
 *     polyvec.c;
 *   - sign-01-5 (challenge signs collapse to one bit) is fixed only when
 *     PQC_NGCC_AIGIS_SIG_STRICT is enabled, because the fix changes every
 *     signature and therefore the published KAT;
 *   - sign-01-7 (DRBG used as KDF) and sign-01-8 (short CRH for sets II/III)
 *     are accepted: the former is required for interoperability with the
 *     submitted KAT, the latter does not affect set I;
 *   - sign-01-9 (32-bit message-length arithmetic) and sign-01-10 (silently
 *     failed DRBG setup) are fixed by 64-bit lengths and status propagation;
 *   - sign-01-11 (leak and short-signature read) is fixed in sign.c and
 *     packing.c; sign-01-12 affects the AVX2 tree, which is not shipped here.
 * sign-01-6 is a specification-only inconsistency.
 *
 * The reference signature has a variable length. The unified API requires a
 * fixed signature field, so the adapter prefixes the packed signature with a
 * two-byte big-endian length, zero-pads it to SIG_MAX_SIZE_PACKED and appends
 * the message. Verification re-checks the length and every padding byte. */

#include "SIG_AlgorithmInstance.h"
#include "api.h"
#include "core/pqc_signature.h"
#include "ngcc/drng.h"
#include "util/rand.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Referenced by sign.c; the KAT driver seeds it directly. */
DRNG_ctx drng_algorithm;

int randombytes(unsigned char *x, size_t xlen)
{
    return get_random_number(&drng_algorithm, x,
                             (unsigned long long)xlen * 8);
}

#define PQC_AIGIS_SIG_NAME "Aigis-Sig+-I"
#define PQC_AIGIS_SIG_PREFIX 2
#define PQC_AIGIS_SIG_FIELD \
    (PQC_AIGIS_SIG_PREFIX + (size_t)SIG_MAX_SIZE_PACKED)

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
    unsigned long long pk_len = sig_get_pk_len_bytes();
    unsigned long long sk_len = sig_get_sk_len_bytes();
    if (seed_drng() != 0)
        return -1;
    return sig_keygen(pk, &pk_len, sk, &sk_len);
}

static int sign_adapter(uint8_t *out, size_t *len, const uint8_t *m, size_t n,
                        const uint8_t *sk)
{
    unsigned long long sn_len = SIG_MAX_SIZE_PACKED;
    if (n > SIZE_MAX - PQC_AIGIS_SIG_FIELD)
        return -1;
    if (sig_sign((unsigned char *)sk, sig_get_sk_len_bytes(),
                 (unsigned char *)m, n, out + PQC_AIGIS_SIG_PREFIX,
                 &sn_len) != 0)
        return -1;
    if (sn_len < SIG_MIN_SIZE_PACKED || sn_len > SIG_MAX_SIZE_PACKED)
        return -1;
    out[0] = (uint8_t)(sn_len >> 8);
    out[1] = (uint8_t)sn_len;
    memset(out + PQC_AIGIS_SIG_PREFIX + sn_len, 0,
           (size_t)SIG_MAX_SIZE_PACKED - (size_t)sn_len);
    memmove(out + PQC_AIGIS_SIG_FIELD, m, n);
    *len = PQC_AIGIS_SIG_FIELD + n;
    return 0;
}

static int open_adapter(uint8_t *out, size_t *len, const uint8_t *sm, size_t n,
                        const uint8_t *pk)
{
    size_t sn_len, i;
    if (n < PQC_AIGIS_SIG_FIELD)
        return -1;
    sn_len = ((size_t)sm[0] << 8) | sm[1];
    if (sn_len < SIG_MIN_SIZE_PACKED || sn_len > SIG_MAX_SIZE_PACKED)
        return -1;
    for (i = sn_len; i < (size_t)SIG_MAX_SIZE_PACKED; i++)
        if (sm[PQC_AIGIS_SIG_PREFIX + i] != 0)
            return -1;
    if (sig_verify((unsigned char *)pk, sig_get_pk_len_bytes(),
                   (unsigned char *)(sm + PQC_AIGIS_SIG_PREFIX), sn_len,
                   (unsigned char *)(sm + PQC_AIGIS_SIG_FIELD),
                   n - PQC_AIGIS_SIG_FIELD) != 0)
        return -1;
    *len = n - PQC_AIGIS_SIG_FIELD;
    memmove(out, sm + PQC_AIGIS_SIG_FIELD, *len);
    return 0;
}

static const pqc_signature algorithm = {
    PQC_AIGIS_SIG_NAME,
    SIG_PUBLICKEYBYTES,
    SIG_SECRETKEYBYTES,
    PQC_AIGIS_SIG_FIELD,
    keypair,
    sign_adapter,
    open_adapter,
};

const pqc_signature *pqc_ngcc_aigis_sig_algorithm(void)
{
    return &algorithm;
}
