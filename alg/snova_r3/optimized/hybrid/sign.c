// SPDX-License-Identifier: MIT

/**
 * Interface to NIST API.
 *
 * This file is the only point where randombytes is used by SNOVA.
 * The snova_* implementations are deterministic.
 *
 * SNOVA Team 2025
 */

#include <string.h>

#if defined(VALGRIND)
#include <valgrind/memcheck.h>
#endif

#include "api.h"
#include "rng.h"
#include "symmetric.h"

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
	uint8_t seed[SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE];

	randombytes(seed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);

#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_UNDEFINED(seed + SEED_LENGTH_PUBLIC, SEED_LENGTH_PRIVATE);
#endif

	int res = SNOVA_NAMESPACE(genkeys)(pk, sk, seed);

#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_DEFINED(pk, BYTES_PK);
#if HASH_PK
	VALGRIND_MAKE_MEM_DEFINED(sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, BYTES_PK_HASH);
#endif
#endif

	SNOVA_SECURE_CLEAR(seed);

	return res;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen, const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk) {
	int res;

#if SNOVA_OPT == 5
	expanded_SK *pskx = (expanded_SK *)sk;
#else
	expanded_SK skx_d;
	expanded_SK *pskx = &skx_d;
	res = SNOVA_NAMESPACE(sk_expand)(pskx, sk);
	if (res) {
		return res;
	}
#endif

	uint8_t salt[BYTES_SALT];
	randombytes(salt, BYTES_SALT);

#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_UNDEFINED(salt, BYTES_SALT);
#endif

	uint8_t digest[BYTES_DIGEST];
	shake256(digest, BYTES_DIGEST, m, mlen);

#if SNOVA_OPT == 5
	memmove(sm + CRYPTO_BYTES, m, mlen);
	*smlen = mlen + CRYPTO_BYTES;
	res = SNOVA_NAMESPACE(sign)(pskx, sm, digest, BYTES_DIGEST, salt);
#else
	uint8_t sig[CRYPTO_BYTES];
	res = SNOVA_NAMESPACE(sign)(pskx, sig, digest, BYTES_DIGEST, salt);

	if (!res) {
		memmove(sm + CRYPTO_BYTES, m, mlen);
		memcpy(sm, sig, CRYPTO_BYTES);
		*smlen = mlen + CRYPTO_BYTES;
	}

	SNOVA_SECURE_CLEARI(skx_d);
#endif

#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_DEFINED(sm, CRYPTO_BYTES);
#endif

	return res;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen, const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk) {
	if (smlen < CRYPTO_BYTES) {
		return -1;
	}

	int res;

#if SNOVA_OPT == 5
	expanded_PK *ppkx = (expanded_PK *)pk;
#else
	expanded_PK pkx;
	expanded_PK *ppkx = &pkx;
	res = SNOVA_NAMESPACE(pk_expand)(ppkx, pk);
	if (res) {
		return -1;
	}
#endif

	uint8_t digest[BYTES_DIGEST];
	shake256(digest, BYTES_DIGEST, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);

	res = SNOVA_NAMESPACE(verify)(ppkx, sm, digest, BYTES_DIGEST);
	if (!res) {
		memmove(m, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
		*mlen = smlen - CRYPTO_BYTES;
	} else {
		return -1;
	}

	return 0;
}
