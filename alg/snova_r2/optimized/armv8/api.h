#ifndef api_h
#define api_h

#include "snova.h"

#define CRYPTO_SECRETKEYBYTES SEED_LENGTH
#define CRYPTO_PUBLICKEYBYTES BYTES_PK
#define CRYPTO_BYTES BYTES_SIGNATURE

#ifdef AESCTR
#define SNOVA_XOF _AES
#else
#define SNOVA_XOF
#endif

#if SNOVA_l == SNOVA_r
#define SNOVA_JOIN_(a, b, c, d, e) SNOVA_##a##_##b##_##c##_##d##e
#define SNOVA_JOIN(a, b, c, d, e) SNOVA_JOIN_(a, b, c, d, e)
#define CRYPTO_ALGNAME_ SNOVA_JOIN(SNOVA_v, SNOVA_o, SNOVA_q, SNOVA_l, SNOVA_XOF)
#else
#define SNOVA_JOIN_(a, b, c, d, r, e) SNOVA_##a##_##b##_##c##_##d##x##r##e
#define SNOVA_JOIN(a, b, c, d, r, e) SNOVA_JOIN_(a, b, c, d, r, e)
#define CRYPTO_ALGNAME_ SNOVA_JOIN(SNOVA_v, SNOVA_o, SNOVA_q, SNOVA_l, SNOVA_r, SNOVA_XOF)
#endif
#define str(s) #s
#define xstr(s) str(s)
#define CRYPTO_ALGNAME xstr(CRYPTO_ALGNAME_)

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk);

int crypto_sign(unsigned char* sm, unsigned long long* smlen, const unsigned char* m, unsigned long long mlen,
                const unsigned char *sk);

int crypto_sign_open(unsigned char* m, unsigned long long* mlen, const unsigned char* sm, unsigned long long smlen,
                     const unsigned char *pk);

#endif
