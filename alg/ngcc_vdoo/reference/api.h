#ifndef VDOO_128_API_H
#define VDOO_128_API_H

#include <stddef.h>
#include "vdoo_keypair.h"

#define CRYPTO_SECRETKEYBYTES sizeof(sk_t)
#define CRYPTO_PUBLICKEYBYTES sizeof(pk_t)

#define CRYPTO_ALGNAME "VDOO-128"

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk);

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk);

#endif /* !VDOO_128_API_H */
