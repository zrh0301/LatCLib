/* Adapter from the ICCS NGCC signature ABI to the library signature API.
 *
 * Reference sources: ICCS NGCC round-1 submission package, "TSUOV"
 * (digital signature track), Reference_Implementation/TSUOV_128. TSUOV is a
 * multivariate construction over GF(31) with L = 2 extension field towers;
 * -DTSUOV_VARIANT selects 128, 256 or 512.
 *
 * ICCS NGCC public review, sign-31:
 *   - sign-31-2 (secret-dependent echelon pivots) is a local timing channel of
 *     the scalar reference elimination. A masked elimination would change the
 *     whole linear-algebra layer and is not attempted here; the reference is
 *     kept and this residual risk is documented.
 *   - sign-31-3 (stale message representative after a failed allocation) and
 *     sign-31-5 (length-addition overflow) are fixed in tsuov_core.c: the
 *     message buffer is sized after an overflow check and a failure is
 *     propagated through TSUOV_Sign/TSUOV_Verify.
 *   - sign-31-4 (ignored declared key/signature lengths) is fixed in
 *     SIG_AlgorithmInstance.c, which now requires exact lengths and therefore
 *     also rejects trailing bytes.
 *   - sign-31-6 (non-canonical 5-bit encodings) is fixed by rejecting the
 *     unused code 31 and non-zero public-key padding bits.
 *   - sign-31-8 (ignored placeholder-XOF failures) is fixed by propagating the
 *     status of pseudoXOF/pseudohash and of every PRG allocation.
 *   - sign-31-1 (512-bit prehash) and sign-31-7 (128-bit salt bound) are
 *     specification/proof issues with no code-level fix.
 *
 * The signature has a fixed length, so the adapter stores the raw signature
 * followed by the message. */

#include "SIG_AlgorithmInstance.h"
#include "core/pqc_signature.h"
#include "ngcc/drng.h"
#include "util/rand.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Referenced by SIG_AlgorithmInstance.c; the KAT driver seeds it directly. */
DRNG_ctx drng_algorithm;

#if TSUOV_VARIANT == 128
#define PQC_TSUOV_NAME "TSUOV-128"
#elif TSUOV_VARIANT == 256
#define PQC_TSUOV_NAME "TSUOV-256"
#elif TSUOV_VARIANT == 512
#define PQC_TSUOV_NAME "TSUOV-512"
#else
#error "unsupported TSUOV_VARIANT"
#endif

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
    unsigned long long sn_len = sig_get_sn_len_bytes();
    if (n > SIZE_MAX - CRYPTO_BYTES)
        return -1;
    if (sig_sign((unsigned char *)sk, sig_get_sk_len_bytes(),
                 (unsigned char *)m, n, out, &sn_len) != 0)
        return -1;
    if (sn_len != CRYPTO_BYTES)
        return -1;
    memmove(out + sn_len, m, n);
    *len = (size_t)sn_len + n;
    return 0;
}

static int open_adapter(uint8_t *out, size_t *len, const uint8_t *sm, size_t n,
                        const uint8_t *pk)
{
    const size_t sn_len = (size_t)CRYPTO_BYTES;
    if (n < sn_len)
        return -1;
    if (sig_verify((unsigned char *)pk, sig_get_pk_len_bytes(),
                   (unsigned char *)sm, sn_len,
                   (unsigned char *)(sm + sn_len), n - sn_len) != 0)
        return -1;
    *len = n - sn_len;
    memmove(out, sm + sn_len, *len);
    return 0;
}

static const pqc_signature algorithm = {
    PQC_TSUOV_NAME,
    CRYPTO_PUBLICKEYBYTES,
    CRYPTO_SECRETKEYBYTES,
    CRYPTO_BYTES, /* fixed signature; the signed message follows it */
    keypair,
    sign_adapter,
    open_adapter,
};

const pqc_signature *pqc_ngcc_tsuov_algorithm(void)
{
    return &algorithm;
}
