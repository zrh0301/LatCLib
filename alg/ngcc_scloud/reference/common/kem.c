/**
 * @file kem.c
 * @brief Fujisaki-Okamoto style KEM transform for Scloud+.
 *
 * This file is shared by all Scloud+ instances.  It wraps the deterministic
 * PKE core with hashing, re-encryption checking, and constant-time fallback
 * key selection.  Instance-specific byte lengths are taken from the leaf
 * parameters.h through scloudplus_param_common.h, so AES/SHAKE instances and
 * SM3 instances can use different hash/seed sizes without changing this
 * control flow.
 */

#include "kem.h"
#include "pke.h"
#include "scloudplus_param_common.h"
#include "hash.h"
#include "util.h"
#include "random.h"
#include <string.h>

/**
 * @brief Generate a KEM key pair.
 *
 * Layout of the secret key:
 *
 *   pke_sk || pk || H(pk) || z
 *
 * `z` is the rejection secret used if decapsulation fails the re-encryption
 * check.  It is sampled independently from the PKE secret so that malformed
 * ciphertexts produce pseudorandom fallback shared secrets.
 */
int scloud_kemkeygen(uint8_t *pk, uint8_t *sk)
{
	uint8_t z[scloudplus_kem_z_bytes];
	int rc = -1;

	if (randombytes(z, scloudplus_kem_z_bytes) != 0 ||
		pke_keygen(pk, sk) != 0)
	{
		goto cleanup;
	}
	memcpy(sk + scloudplus_pke_sk, pk, scloudplus_pk);
	scloudplus_H(sk + scloudplus_pke_sk + scloudplus_pk, pk, scloudplus_pk);
	memcpy(sk + scloudplus_kem_sk - scloudplus_kem_z_bytes, z,
		   scloudplus_kem_z_bytes);
	rc = 0;

cleanup:
	scloudplus_secure_zeroize(z, sizeof(z));
	return rc;
}

/**
 * @brief Encapsulate a shared secret to a public key.
 *
 * The random message m is bound to the public key by hashing H(pk) into the
 * input of G.  G expands m || H(pk) into:
 *
 * - PKE encryption randomness (first hash_bytes);
 * - key-confirmation material used with the ciphertext in K.
 */
int scloud_kemencaps(const uint8_t *pk, uint8_t *ctx, uint8_t *ss)
{
	uint8_t kc[scloudplus_ctx + scloudplus_hash_bytes];
	uint8_t m[scloudplus_ss + scloudplus_hash_bytes];
	uint8_t rk[scloudplus_G_bytes];
	int rc = -1;

	if (randombytes(m, scloudplus_ss) != 0)
	{
		goto cleanup;
	}
	scloudplus_H(m + scloudplus_ss, pk, scloudplus_pk);
	scloudplus_G(rk, m, scloudplus_ss + scloudplus_hash_bytes);
	pke_enc(pk, m, rk, ctx);
	memcpy(kc, rk + scloudplus_hash_bytes, scloudplus_hash_bytes);
	memcpy(kc + scloudplus_hash_bytes, ctx, scloudplus_ctx);
	scloudplus_K(ss, scloudplus_ss, kc, scloudplus_ctx + scloudplus_hash_bytes);
	rc = 0;

cleanup:
	scloudplus_secure_zeroize(kc, sizeof(kc));
	scloudplus_secure_zeroize(m, sizeof(m));
	scloudplus_secure_zeroize(rk, sizeof(rk));
	return rc;
}

/**
 * @brief Decapsulate and select the fallback key in constant time.
 *
 * Decapsulation decrypts m, recomputes the encapsulation ciphertext, and
 * compares it with the received ciphertext.  If comparison fails, the input to
 * K uses z from the secret key instead of the recomputed key material.  The
 * conditional move is branch-free with respect to the failure bit.
 */
int scloud_kemdecaps(const uint8_t *sk, const uint8_t *ctx, uint8_t *ss)
{
	uint8_t m[scloudplus_ss + scloudplus_hash_bytes];
	uint8_t rk[scloudplus_G_bytes];
	uint8_t ctx1[scloudplus_ctx + scloudplus_hash_bytes];
	uint8_t fail_mask;

	pke_dec(sk, ctx, m);
	memcpy(m + scloudplus_ss, sk + scloudplus_pke_sk + scloudplus_pk,
		   scloudplus_hash_bytes);
	scloudplus_G(rk, m, scloudplus_ss + scloudplus_hash_bytes);
	pke_enc(sk + scloudplus_pke_sk, m, rk, ctx1 + scloudplus_hash_bytes);
	fail_mask = scloudplus_verify(ctx, ctx1 + scloudplus_hash_bytes, scloudplus_ctx);
	memcpy(ctx1 + scloudplus_hash_bytes, ctx, scloudplus_ctx);
	scloudplus_cmov(ctx1, rk + scloudplus_hash_bytes,
					sk + scloudplus_kem_sk - scloudplus_kem_z_bytes,
					scloudplus_hash_bytes, fail_mask);
	scloudplus_K(ss, scloudplus_ss, ctx1, scloudplus_ctx + scloudplus_hash_bytes);
	scloudplus_secure_zeroize(m, sizeof(m));
	scloudplus_secure_zeroize(rk, sizeof(rk));
	scloudplus_secure_zeroize(ctx1, sizeof(ctx1));
	return 0;
}
