#include "ngcc/auxfunc.h"
#include "api.h"

#include "vdoo_keypair.h"
#include "vdoo_sign.h"
#include "vdoo_verif.h"
#include "vdoo_config.h"


int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    unsigned char sk_seed[LEN_SKSEED] = { 0 };
    get_random_number(&pqc_ngcc_vdoo_drng, sk_seed, LEN_SKSEED * 8ULL);

    int r = generate_keypair((pk_t*)pk, (sk_t*)sk, sk_seed);

    for(int i=0;i<LEN_SKSEED;i++)
        sk_seed[i]=0;
    return r;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk) {
    unsigned char digest[HASH_LEN] = {0};
    if (pseudoXOF(HASH_LEN * 8, m, mlen * 8, digest) != 0)
        return -1;

    int r = -1;
    r = vdoo_sign(sm + mlen, (const sk_t *)sk, digest);

    memcpy(sm, m, mlen);
    smlen[0] = mlen + CRYPTO_BYTES;

    return r;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk) {
    if (CRYPTO_BYTES > smlen)
        return -1;

    unsigned long long msg_len = smlen - CRYPTO_BYTES;
    const uint8_t *signature = sm + msg_len;

    unsigned char digest[HASH_LEN] = {0};
    if (pseudoXOF(HASH_LEN * 8, sm, msg_len * 8, digest) != 0)
        return -1;

    int r = vdoo_verify(digest, signature, (const pk_t *)pk);

    if (r == 0 && m != NULL)
    {
        memcpy(m, sm, msg_len);
        *mlen = msg_len;
    }

    return r;
}
