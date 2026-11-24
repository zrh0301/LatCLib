/*
The software is provided by the Institute of Commercial Cryptography Standards
(ICCS), and is used for algorithm submissions in the Next-generation Commercial
Cryptographic Algorithms Program (NGCC).

ICCS doesn't represent or warrant that the operation of the software will be
uninterrupted or error-free in all cases. ICCS will take no responsibility for
the use of the software or the results thereof, if the software is used for any
other purposes.
*/

#include "SIG_AlgorithmInstance.h"
#include "tsuov_core.h"
#include "ngcc/drng.h"
#include <string.h>

extern DRNG_ctx drng_algorithm;

static int sig_random_seed(TSUOV_SEED seed)
{
    return get_random_number(&drng_algorithm, seed,
                             (unsigned long long)TSUOV_SEED_LEN * 8);
}

unsigned long long sig_get_pk_len_bytes(void)
{
    return CRYPTO_PUBLICKEYBYTES;
}

unsigned long long sig_get_sk_len_bytes(void)
{
    return CRYPTO_SECRETKEYBYTES;
}

unsigned long long sig_get_sn_len_bytes(void)
{
    return CRYPTO_BYTES;
}

int sig_keygen(unsigned char *pk, unsigned long long *pk_len_bytes,
               unsigned char *sk, unsigned long long *sk_len_bytes)
{
    TSUOV_SEED seed_sk;
    TSUOV_SEED seed_pk;
    static TSUOV_P3 P3;
    int rc;

    if (sig_random_seed(seed_sk) != 0)
        return -1;
    if (sig_random_seed(seed_pk) != 0)
        return -1;

    rc = TSUOV_KeyGen(seed_sk, seed_pk, P3);
    if (rc != 0)
        return -1;

    /* Clear the output buffers so unused padding bits in the encoded public
     * key are zero (sign-31-6). */
    memset(sk, 0, CRYPTO_SECRETKEYBYTES);
    memset(pk, 0, CRYPTO_PUBLICKEYBYTES);

    size_t sk_pool_bits = 0;
    store_TSUOV_SEED(seed_sk, sk, &sk_pool_bits);
    store_TSUOV_SEED(seed_pk, sk, &sk_pool_bits);
    *sk_len_bytes = CRYPTO_SECRETKEYBYTES;

    size_t pk_pool_bits = 0;
    store_TSUOV_SEED(seed_pk, pk, &pk_pool_bits);
    store_TSUOV_P3(P3, pk, &pk_pool_bits);
    *pk_len_bytes = CRYPTO_PUBLICKEYBYTES;

    return 0;
}

int sig_sign(unsigned char *sk, unsigned long long sk_len_bytes,
             unsigned char *m, unsigned long long m_len_bytes,
             unsigned char *sn, unsigned long long *sn_len_bytes)
{
    if (sk_len_bytes != CRYPTO_SECRETKEYBYTES)
        return -1;

    size_t sk_pool_bits = 0;
    TSUOV_SEED seed_sk;
    TSUOV_SEED seed_pk;
    restore_TSUOV_SEED(sk, &sk_pool_bits, seed_sk);
    restore_TSUOV_SEED(sk, &sk_pool_bits, seed_pk);

    TSUOV_SEED seed_v;
    TSUOV_SEED seed_r;
    TSUOV_SEED seed_sol;
    if (sig_random_seed(seed_v) != 0)
        return -1;
    if (sig_random_seed(seed_r) != 0)
        return -1;
    if (sig_random_seed(seed_sol) != 0)
        return -1;

    TSUOV_SIGNATURE sig;
    if (TSUOV_Sign(seed_sk, seed_pk, seed_v, seed_r, seed_sol,
                   m, (size_t)m_len_bytes, sig) != 0)
        return -1;

    size_t sig_pool_bits = 0;
    memset(sn, 0, CRYPTO_BYTES);
    store_TSUOV_SIGNATURE(sig, sn, &sig_pool_bits);
    *sn_len_bytes = CRYPTO_BYTES;

    return 0;
}

int sig_verify(unsigned char *pk, unsigned long long pk_len_bytes,
               unsigned char *sn, unsigned long long sn_len_bytes,
               unsigned char *m, unsigned long long m_len_bytes)
{
    /* sign-31-4: parse the declared lengths instead of ignoring them. */
    if (pk_len_bytes != CRYPTO_PUBLICKEYBYTES ||
        sn_len_bytes != CRYPTO_BYTES)
        return -1;

    size_t pk_pool_bits = 0;
    TSUOV_SEED seed_pk;
    static TSUOV_P3 P3;

    restore_TSUOV_SEED(pk, &pk_pool_bits, seed_pk);
    /* sign-31-6: reject the unused field code 31 and non-zero padding. */
    if (restore_TSUOV_P3(pk, &pk_pool_bits, P3) != 0)
        return -1;

    size_t sig_pool_bits = 0;
    TSUOV_SIGNATURE sig;
    if (restore_TSUOV_SIGNATURE(sn, &sig_pool_bits, sig) != 0)
        return -1;

    int rc = TSUOV_Verify(seed_pk, P3, m, (size_t)m_len_bytes, sig);
    if (rc < 0)
        return -2;
    return rc ? 0 : -1;
}
