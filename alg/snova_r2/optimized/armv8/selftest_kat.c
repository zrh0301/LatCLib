#include "api.h"
#include "symmetric.h"

#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define KAT_MAX_MSG 96
#define KAT_DEFAULT_ITERS 32

static uint64_t kat_rng_state;

static void kat_rng_reset(void) {
	kat_rng_state = 0x6a09e667f3bcc909ull;
}

static uint64_t kat_rng_next_u64(void) {
	uint64_t z;
	kat_rng_state += 0x9e3779b97f4a7c15ull;
	z = kat_rng_state;
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
	return z ^ (z >> 31);
}

static void kat_rng_fill(unsigned char *out, size_t n) {
	size_t offset = 0;
	while (offset < n) {
		uint64_t word = kat_rng_next_u64();
		for (size_t i = 0; i < 8 && offset < n; ++i, ++offset) {
			out[offset] = (unsigned char)(word >> (8 * i));
		}
	}
}

void randombytes_init(unsigned char *entropy_input, unsigned char *personalization_string, int security_strength) {
	(void)entropy_input;
	(void)personalization_string;
	(void)security_strength;
	kat_rng_reset();
}

int randombytes(unsigned char *x, unsigned long long xlen) {
	kat_rng_fill(x, (size_t)xlen);
	return 0;
}

static size_t iters_from_env(void) {
	const char *env = getenv("SNOVA_KAT_ITERS");
	if (env == NULL || env[0] == '\0') {
		return KAT_DEFAULT_ITERS;
	}

	char *end = NULL;
	unsigned long parsed = strtoul(env, &end, 10);
	if (end == env || parsed == 0) {
		return KAT_DEFAULT_ITERS;
	}
	return (size_t)parsed;
}

static int exhaustive_from_env(void) {
	const char *env = getenv("SNOVA_KAT_EXHAUSTIVE");
	return env != NULL && env[0] != '\0' && env[0] != '0';
}

int main(void) {
	static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
	static uint8_t sk[CRYPTO_SECRETKEYBYTES];
	static uint8_t sm[CRYPTO_BYTES + KAT_MAX_MSG];
	static uint8_t msg[KAT_MAX_MSG];
	static uint8_t out[CRYPTO_BYTES + KAT_MAX_MSG];

	const size_t iters = iters_from_env();
	const int exhaustive = exhaustive_from_env();
	size_t failures = 0;
	shake_t acc;

	kat_rng_reset();
	shake256_init(&acc);

	for (size_t i = 0; i < iters; ++i) {
		unsigned long long smlen = 0;
		unsigned long long outlen = 0;
		const size_t mlen = 1 + (i * 7) % KAT_MAX_MSG;

		for (size_t j = 0; j < mlen; ++j) {
			msg[j] = (uint8_t)(i * 31u + j * 17u + 5u);
		}

		if (crypto_sign_keypair(pk, sk) != 0) {
			printf("keypair failed at iteration %zu\n", i);
			failures++;
			break;
		}

		if (crypto_sign(sm, &smlen, msg, mlen, sk) != 0) {
			printf("sign failed at iteration %zu\n", i);
			failures++;
			break;
		}

		shake_absorb(&acc, pk, CRYPTO_PUBLICKEYBYTES);
		shake_absorb(&acc, sk, CRYPTO_SECRETKEYBYTES);
		shake_absorb(&acc, sm, (size_t)smlen);

		if (crypto_sign_open(out, &outlen, sm, smlen, pk) != 0) {
			printf("open rejected a valid signature at iteration %zu\n", i);
			failures++;
			continue;
		}
		if (outlen != (unsigned long long)mlen || memcmp(out, msg, mlen) != 0) {
			printf("open returned the wrong message at iteration %zu\n", i);
			failures++;
			continue;
		}

		sm[i % CRYPTO_BYTES] ^= 0x01u;
		if (crypto_sign_open(out, &outlen, sm, smlen, pk) == 0) {
			printf("open accepted a tampered signature at iteration %zu\n", i);
			failures++;
		}
		sm[i % CRYPTO_BYTES] ^= 0x01u;

		sm[CRYPTO_BYTES + (i % mlen)] ^= 0x01u;
		if (crypto_sign_open(out, &outlen, sm, smlen, pk) == 0) {
			printf("open accepted a tampered message at iteration %zu\n", i);
			failures++;
		}
		sm[CRYPTO_BYTES + (i % mlen)] ^= 0x01u;

		pk[i % SEED_LENGTH_PUBLIC] ^= 0x01u;
		if (crypto_sign_open(out, &outlen, sm, smlen, pk) == 0) {
			printf("open accepted a tampered public key at iteration %zu\n", i);
			failures++;
		}
		pk[i % SEED_LENGTH_PUBLIC] ^= 0x01u;

		if (exhaustive && i == 0) {
			for (size_t pos = 0; pos < CRYPTO_BYTES; ++pos) {
				for (unsigned bit = 0; bit < 8; ++bit) {
					sm[pos] ^= (uint8_t)(1u << bit);
					if (crypto_sign_open(out, &outlen, sm, smlen, pk) == 0) {
						printf("open accepted signature mutation pos=%zu bit=%u\n", pos, bit);
						failures++;
					}
					sm[pos] ^= (uint8_t)(1u << bit);
				}
			}
			for (size_t pos = 0; pos < SEED_LENGTH_PUBLIC; ++pos) {
				for (unsigned bit = 0; bit < 8; ++bit) {
					pk[pos] ^= (uint8_t)(1u << bit);
					if (crypto_sign_open(out, &outlen, sm, smlen, pk) == 0) {
						printf("open accepted public-seed mutation pos=%zu bit=%u\n", pos, bit);
						failures++;
					}
					pk[pos] ^= (uint8_t)(1u << bit);
				}
			}
		}

		if (crypto_sign_open(out, &outlen, sm, CRYPTO_BYTES - 1, pk) == 0) {
			printf("open accepted a truncated signed message at iteration %zu\n", i);
			failures++;
		}
	}

	uint8_t digest[32];
	shake_finalize(&acc);
	shake_squeeze(digest, sizeof(digest), &acc);

	printf("alg        : %s\n", CRYPTO_ALGNAME);
	printf("iterations : %zu\n", iters);
	printf("exhaustive : %s\n", exhaustive ? "yes" : "no");
	printf("sizes      : pk=%d sk=%d sig=%d\n", CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, CRYPTO_BYTES);
	printf("digest     : ");
	for (size_t i = 0; i < sizeof(digest); ++i) {
		printf("%02x", digest[i]);
	}
	printf("\n");
	printf("failures   : %zu\n", failures);

	return failures == 0 ? 0 : 1;
}
