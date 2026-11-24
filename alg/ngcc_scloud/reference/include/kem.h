/**
 * @file kem.h
 * @brief Internal Scloud+ KEM interface used by the API_PKC wrapper.
 */

#ifndef _SCLOUDPLUS_KEM_H_
#define _SCLOUDPLUS_KEM_H_

#include "scloudplus_param_common.h"
#include <stdint.h>

/**
 * @brief Generate a Scloud+ KEM public key and secret key.
 *
 * The secret-key layout is `ShortPack(S) || pk || H(pk) || z`, with all byte
 * lengths controlled by the checked-in leaf `parameters.h`.
 *
 * @param[out] pk Public key buffer of `scloudplus_pk` bytes.
 * @param[out] sk Secret key buffer of `scloudplus_kem_sk` bytes.
 *
 * @return `0` on success and `-1` if runtime randomness fails.
 */
int scloud_kemkeygen(uint8_t *pk, uint8_t *sk);

/**
 * @brief Encapsulate to a Scloud+ public key.
 *
 * @param[in]  pk  Public key buffer of `scloudplus_pk` bytes.
 * @param[out] ctx Ciphertext buffer of `scloudplus_ctx` bytes.
 * @param[out] ss  Shared-secret buffer of `scloudplus_ss` bytes.
 *
 * @return `0` on success and `-1` if runtime randomness fails.
 */
int scloud_kemencaps(const uint8_t *pk, uint8_t *ctx, uint8_t *ss);

/**
 * @brief Decapsulate a Scloud+ ciphertext.
 *
 * Decapsulation always writes a shared secret. If the re-encryption check
 * fails, the implementation uses the secret fallback value `z` through a
 * constant-time conditional move.
 *
 * @param[in]  sk  Secret key buffer of `scloudplus_kem_sk` bytes.
 * @param[in]  ctx Ciphertext buffer of `scloudplus_ctx` bytes.
 * @param[out] ss  Shared-secret buffer of `scloudplus_ss` bytes.
 *
 * @return Always `0`. Implicit rejection is handled by constant-time selection
 * of the KDF input and is not exposed through a secret-dependent return code.
 */
int scloud_kemdecaps(const uint8_t *sk, const uint8_t *ctx, uint8_t *ss);

#endif
