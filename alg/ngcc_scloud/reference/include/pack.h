/**
 * @file pack.h
 * @brief Serialization interface for public keys, secret keys, and ciphertexts.
 */

#ifndef SCLOUDPLUS_PACK_H
#define SCLOUDPLUS_PACK_H
#include <stdint.h>

/**
 * @brief Serialize the short secret matrix S with two bits per coefficient.
 *
 * This is the implementation's `ShortPack(S)` format. It is distinct from
 * Pack10 because S only contains coefficients in `{-1,0,1}`.  The two stored
 * bits are the coefficient's low two bits in two's-complement form:
 * `0 -> 00`, `1 -> 01`, and `-1 -> 11`.
 */
void pack_sk(const uint16_t *S, uint8_t *sk);

/**
 * @brief Parse `ShortPack(S)` back into the internal q-ary representation.
 */
void unpack_sk(const uint8_t *sk, uint16_t *S);

/**
 * @brief Serialize the public matrix B in row-major Pack10 order.
 */
void pack_pk(const uint16_t *B, uint8_t *pk);

/**
 * @brief Parse the public matrix B from row-major Pack10 order.
 */
void unpack_pk(const uint8_t *pk, uint16_t *B);

/**
 * @brief Reduce C1 coefficients to canonical representatives modulo q.
 */
void reduce_c1(const uint16_t *C, uint16_t *out);

/**
 * @brief Reduce C2 coefficients to canonical representatives modulo q.
 */
void reduce_c2(const uint16_t *C, uint16_t *out);

/**
 * @brief Serialize ciphertext component C1 in row-major Pack10 order.
 */
void pack_c1(const uint16_t *C, uint8_t *out);

/**
 * @brief Parse ciphertext component C1 from row-major Pack10 order.
 */
void unpack_c1(const uint8_t *in, uint16_t *C);

/**
 * @brief Serialize ciphertext component C2 in row-major Pack10 order.
 */
void pack_c2(const uint16_t *C, uint8_t *out);

/**
 * @brief Parse ciphertext component C2 from row-major Pack10 order.
 */
void unpack_c2(const uint8_t *in, uint16_t *C);
#endif
