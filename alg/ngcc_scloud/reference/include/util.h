/**
 * @file util.h
 * @brief Utility declarations for byte comparison and constant-time selection.
 */

#ifndef _SCLOUDPLUS_UTIL_H_
#define _SCLOUDPLUS_UTIL_H_
#include <stddef.h>
#include <stdint.h>

/**
 * @brief Constant-time equality test for two byte strings.
 *
 * @return `0x00` if the byte strings are equal, and `0xff` otherwise.
 */
uint8_t scloudplus_verify(const uint8_t *a, const uint8_t *b, size_t len);

/**
 * @brief Constant-time conditional move between two byte strings.
 *
 * If `mask` is `0xff`, bytes from `b` are selected. If `mask` is `0x00`,
 * bytes from `a` are kept. Other mask values are not used by the KEM path.
 */
void scloudplus_cmov(uint8_t *r, const uint8_t *a, const uint8_t *b, size_t len,
					 uint8_t mask);

/**
 * @brief Clear a byte buffer in a way the compiler must keep.
 *
 * This helper is used for stack buffers that contain seeds, secret matrices, or
 * FO transform intermediate material.  GCC/Clang builds use a normal memset
 * followed by a compiler memory barrier, so the wipe can use the platform's
 * optimized clearing sequence without being removed as a dead store.  Other
 * compilers fall back to volatile byte stores.
 */
void scloudplus_secure_zeroize(void *ptr, size_t len);
#endif
