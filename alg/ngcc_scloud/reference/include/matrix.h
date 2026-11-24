/**
 * @file matrix.h
 * @brief Matrix arithmetic interface for public matrix, ciphertext, and key operations.
 */

#ifndef SCLOUDPLUS_MATRIX_H
#define SCLOUDPLUS_MATRIX_H
#include <stdint.h>

/**
 * @brief Add two row-major q-ary vectors coefficient-wise modulo `q = 2^10`.
 *
 * @param[in]  lhs Left input vector.
 * @param[in]  rhs Right input vector.
 * @param[in]  len Number of coefficients.
 * @param[out] out Output vector; may alias either input.
 */
void mat_add(uint16_t *lhs, uint16_t *rhs, int len, uint16_t *out);

/**
 * @brief Subtract two row-major q-ary vectors coefficient-wise modulo `q`.
 *
 * @param[in]  lhs Left input vector.
 * @param[in]  rhs Right input vector.
 * @param[in]  len Number of coefficients.
 * @param[out] out Output vector; may alias either input.
 */
void mat_sub(uint16_t *lhs, uint16_t *rhs, int len, uint16_t *out);

/**
 * @brief Compute the public-key matrix `B = A*S + E`.
 *
 * `A` is generated row-by-row from `seedA`. The build-selected backend may
 * use scalar code, AES instructions, AVX2, ARM Crypto, or NEON, but the output
 * is always the same row-major `m x nbar` q-ary matrix. `B` may alias `E` for
 * in-place accumulation by the PKE control flow.
 */
void mul_as_e(const uint8_t *seedA, const uint16_t *S, const uint16_t *E, uint16_t *B);

/**
 * @brief Compute the first ciphertext component `C1 = S'*A + E1`.
 *
 * `C` may alias `E` for in-place accumulation.
 */
void mul_sa_e(const uint8_t *seedA, const uint16_t *S, uint16_t *E, uint16_t *C);

/**
 * @brief Compute the second ciphertext product `S'*B + E2`.
 *
 * `out` may alias `E` for in-place accumulation.
 */
void mul_sb_e(const uint16_t *S, const uint16_t *B, const uint16_t *E, uint16_t *out);

/**
 * @brief Compute `C1*S` during PKE decryption.
 */
void mul_cs(uint16_t *C, uint16_t *S, uint16_t *out);

#endif
