/* ICCS NGCC round-1 "VDOO" signature API shim, from the submission's
 * Reference_Implementation/vdoo_128. Local changes: the two shared NGCC
 * headers (drng/auxfunc) are referenced through inc/ngcc, and sig_verify
 * rejects signature buffers whose length differs from CRYPTO_BYTES instead of
 * letting crypto_sign_open reinterpret the tail of the message. */

#include "SIG_AlgorithmInstance.h"
#include "api.h"
#include <stdlib.h>
#include <string.h>

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
    *pk_len_bytes = CRYPTO_PUBLICKEYBYTES;
    *sk_len_bytes = CRYPTO_SECRETKEYBYTES;

    return crypto_sign_keypair(pk, sk);
}

int sig_sign(unsigned char *sk, unsigned long long sk_len_bytes,
             unsigned char *m, unsigned long long m_len_bytes,
             unsigned char *sn, unsigned long long *sn_len_bytes)
{
    unsigned char *sm;
    unsigned long long smlen;
    int ret;

    if (sk_len_bytes != CRYPTO_SECRETKEYBYTES)
        return -2;

    sm = (unsigned char*)malloc(m_len_bytes + CRYPTO_BYTES);
    if (sm == NULL)
        return -1;

    ret = crypto_sign(sm, &smlen, m, m_len_bytes, sk);

    if (ret == 0)
    {
        memcpy(sn, sm + m_len_bytes, CRYPTO_BYTES);
        *sn_len_bytes = CRYPTO_BYTES;
    }

    free(sm);
    return ret;
}

int sig_verify(unsigned char *pk, unsigned long long pk_len_bytes,
               unsigned char *sn, unsigned long long sn_len_bytes,
               unsigned char *m, unsigned long long m_len_bytes)
{
    unsigned char *sm;
    unsigned char *m_out = NULL;
    unsigned long long mlen_out;
    int ret;

    if (pk_len_bytes != CRYPTO_PUBLICKEYBYTES)
        return -2;

    if (sn_len_bytes != CRYPTO_BYTES)
        return -2;

    sm = (unsigned char*)malloc(m_len_bytes + sn_len_bytes);
    if (sm == NULL)
        return -1;

    memcpy(sm, m, m_len_bytes);
    memcpy(sm + m_len_bytes, sn, sn_len_bytes);

    ret = crypto_sign_open(m_out, &mlen_out, sm, m_len_bytes + sn_len_bytes, pk);

    free(sm);
    return ret;
}