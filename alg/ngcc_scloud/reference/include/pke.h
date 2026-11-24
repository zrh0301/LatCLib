/**
 * @file pke.h
 * @brief Internal public-key encryption interface used by the KEM transform.
 */

#ifndef SCLOUDPLUS_PKE_H
#define SCLOUDPLUS_PKE_H
#include <stdint.h>

/**
 * @brief Generate one Scloud+ PKE keypair.
 *
 * The public key stores `Pack10(B) || SeedA`; the PKE secret key stores
 * `ShortPack(S)`. The KEM layer appends `pk`, `H(pk)`, and `z`.
 *
 * @return `0` on success and `-1` if runtime randomness fails.
 */
int pke_keygen(uint8_t *pk, uint8_t *sk);

/**
 * @brief Encrypt a message under a Scloud+ PKE public key.
 *
 * The randomness buffer `coins` is expanded into `r1` and `r2` according to
 * the selected family. The ciphertext is `Pack10(C1) || Pack10(C2)`.
 */
void pke_enc(const uint8_t *pk, const uint8_t *m, const uint8_t *coins, uint8_t *ct);

/**
 * @brief Decrypt a Scloud+ PKE ciphertext into the encoded message bytes.
 */
void pke_dec(const uint8_t *sk, const uint8_t *ct, uint8_t *m);
#endif
