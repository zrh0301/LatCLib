/**
 * @file pke.c
 * @brief Shared PKE control flow for Scloud+.
 *
 * Public-key encryption routines shared by KEM key generation, encapsulation,
 * and decapsulation.  Instance parameters come from the checked-in leaf
 * `parameters.h` through the shared parameter glue; matrix generation,
 * sampling, message encoding, and serialization are delegated to the selected
 * shared modules.  This keeps the PKE byte format stable while allowing either
 * CMake or a direct compiler command to select optimized backend kernels.
 */

#include "pke.h"
#include "scloudplus_param_common.h"
#include "encode.h"
#include "matrix.h"
#include "sample.h"
#include "hash.h"
#include "random.h"
#include "pack.h"
#include "util.h"
#include <string.h>

/*
 * Pack10 and the BW/RBW decoder both consume coefficients modulo q by taking
 * the low 10 bits.  `mat_add`/`mat_sub` also reduce their outputs, so the
 * explicit reduce passes used by early per-instance code paths are redundant.
 */
#define SCLOUDPLUS_PKE_KEEP_REDUCE_BEFORE_PACK 0

#if defined(_WIN32) || defined(_WIN64)
#define SCLOUDPLUS_ALIGN_HEADER(N) __declspec(aligned(N))
#define SCLOUDPLUS_ALIGN_FOOTER(N)
#else
#define SCLOUDPLUS_ALIGN_HEADER(N)
#define SCLOUDPLUS_ALIGN_FOOTER(N) __attribute__((aligned(N)))
#endif

#if ((scloudplus_mbar * scloudplus_n) < (scloudplus_m * scloudplus_nbar))
#error "pke_enc reuses the E1/C1 buffer as B; E1 must be at least as large as B"
#endif

/**
 * @brief Generate a PKE public key and short PKE secret key.
 *
 * `alpha` is expanded by F into `SeedA || r1 || r2`. The sampler uses `r1` for
 * S and `r2` for E. Matrix multiplication computes `B = A*S + E`. The public
 * key is `Pack10(B) || SeedA`, and the PKE secret key is `ShortPack(S)`.
 */
int pke_keygen(uint8_t *pk, uint8_t *sk)
{
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t S[(size_t)scloudplus_n * scloudplus_nbar] SCLOUDPLUS_ALIGN_FOOTER(32);
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t E[(size_t)scloudplus_m * scloudplus_nbar] SCLOUDPLUS_ALIGN_FOOTER(32);
	uint8_t alpha[scloudplus_pke_keygen_coins_bytes];
	uint8_t seed[scloudplus_pke_keygen_expand_bytes];
	uint8_t *seedA = seed;
	uint8_t *r1 = seed + scloudplus_seedA_bytes;
	uint8_t *r2 = r1 + scloudplus_pke_keygen_r1_bytes;
	int rc = -1;

	if (randombytes(alpha, scloudplus_pke_keygen_coins_bytes) != 0)
	{
		goto cleanup;
	}
	scloudplus_F(seed, scloudplus_pke_keygen_expand_bytes, alpha,
				 scloudplus_pke_keygen_coins_bytes);
	sample_s(r1, scloudplus_pke_keygen_r1_bytes, S);
	sample_e(r2, scloudplus_pke_keygen_r2_bytes, E);
	mul_as_e(seedA, S, E, E);
	pack_pk(E, pk);
	memcpy(pk + scloudplus_pk - scloudplus_seedA_bytes, seedA,
		   scloudplus_seedA_bytes);
	pack_sk(S, sk);
	rc = 0;

cleanup:
	scloudplus_secure_zeroize(S, sizeof(S));
	scloudplus_secure_zeroize(E, sizeof(E));
	scloudplus_secure_zeroize(alpha, sizeof(alpha));
	scloudplus_secure_zeroize(seed, sizeof(seed));
	return rc;
}

/**
 * @brief Encrypt a message under a serialized public key.
 *
 * The coins supplied by the KEM layer are expanded into seeds for S', E1, and
 * E2. The ciphertext components are serialized as `Pack10(C1) || Pack10(C2)`
 * in row-major order; Pack10 writes the canonical low 10 bits.
 */
void pke_enc(const uint8_t *pk, const uint8_t *m, const uint8_t *coins, uint8_t *ct)
{
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t S1[(size_t)scloudplus_mbar * scloudplus_m] SCLOUDPLUS_ALIGN_FOOTER(32);
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t E1[(size_t)scloudplus_mbar * scloudplus_n] SCLOUDPLUS_ALIGN_FOOTER(32);
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t E2[(size_t)scloudplus_mbar * scloudplus_nbar] SCLOUDPLUS_ALIGN_FOOTER(32);
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t mu0[(size_t)scloudplus_mbar * scloudplus_nbar] SCLOUDPLUS_ALIGN_FOOTER(32);
	uint8_t seed[scloudplus_pke_enc_expand_bytes];
	const uint8_t *seedA = pk + scloudplus_pk - scloudplus_seedA_bytes;
	uint8_t *r1 = seed;
	uint8_t *r2 = seed + scloudplus_pke_enc_r1_bytes;

	scloudplus_F(seed, scloudplus_pke_enc_expand_bytes, coins,
				 scloudplus_pke_enc_coins_bytes);
	sample_sp(r1, scloudplus_pke_enc_r1_bytes, S1);
	sample_e12(r2, scloudplus_pke_enc_r2_bytes, E1, E2);
	msg_encode(m, mu0);
	mul_sa_e(seedA, S1, E1, E1);
#if SCLOUDPLUS_PKE_KEEP_REDUCE_BEFORE_PACK
	reduce_c1(E1, E1);
#endif
	pack_c1(E1, ct);
	unpack_pk(pk, E1);
	mul_sb_e(S1, E1, E2, E2);
	mat_add(E2, mu0, scloudplus_mbar * scloudplus_nbar, E2);
#if SCLOUDPLUS_PKE_KEEP_REDUCE_BEFORE_PACK
	reduce_c2(E2, E2);
#endif
	pack_c2(E2, ct + scloudplus_c1);
	scloudplus_secure_zeroize(S1, sizeof(S1));
	scloudplus_secure_zeroize(E1, sizeof(E1));
	scloudplus_secure_zeroize(E2, sizeof(E2));
	scloudplus_secure_zeroize(mu0, sizeof(mu0));
	scloudplus_secure_zeroize(seed, sizeof(seed));
}

/**
 * @brief Decrypt a PKE ciphertext and decode the embedded message.
 *
 * The secret matrix is recovered from `ShortPack(S)`. Decryption computes
 * `D = C2 - C1*S` modulo q, then passes D to the selected BW/RBW decoder.
 */
void pke_dec(const uint8_t *sk, const uint8_t *ct, uint8_t *m)
{
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t S[(size_t)scloudplus_n * scloudplus_nbar] SCLOUDPLUS_ALIGN_FOOTER(32);
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t C1[(size_t)scloudplus_mbar * scloudplus_n] SCLOUDPLUS_ALIGN_FOOTER(32);
	SCLOUDPLUS_ALIGN_HEADER(32)
	uint16_t D[(size_t)scloudplus_mbar * scloudplus_nbar] SCLOUDPLUS_ALIGN_FOOTER(32);

	unpack_sk(sk, S);
	unpack_c1(ct, C1);
	mul_cs(C1, S, D);
	unpack_c2(ct + scloudplus_c1, C1);
#if SCLOUDPLUS_PKE_KEEP_REDUCE_BEFORE_PACK
	reduce_c2(C1, C1);
#endif
	mat_sub(C1, D, scloudplus_mbar * scloudplus_nbar, D);
#if SCLOUDPLUS_PKE_KEEP_REDUCE_BEFORE_PACK
	reduce_c2(D, D);
#endif
	msg_decode(D, m);
	scloudplus_secure_zeroize(S, sizeof(S));
	scloudplus_secure_zeroize(C1, sizeof(C1));
	scloudplus_secure_zeroize(D, sizeof(D));
}
