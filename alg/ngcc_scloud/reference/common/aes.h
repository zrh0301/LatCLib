/**
 * @file aes.h
 * @brief AES-128 block and CTR interfaces used for AES-family A generation.
 */

/********************************************************************************************
 * Header defining the APIs for standalone AES
 *********************************************************************************************/

#ifndef __AES_H
#define __AES_H

#include <stdint.h>
#include <stdlib.h>

/**
 * Function to fill a key schedule given an initial key.
 *
 * @param key            Initial Key.
 * @param schedule       Abstract data structure for a key schedule.
 */
void AES128_load_schedule(const uint8_t *key, uint8_t *schedule);

/**
 * Encrypt AES counter blocks whose first 32-bit little-endian word is a
 * consecutive public counter and whose remaining words are zero.
 *
 * This is used by AES-family A generation to avoid materializing the public
 * counter blocks in a temporary buffer before encrypting them.
 *
 * @param counter_start  First public counter value.
 * @param block_count    Number of 16-byte counter blocks to encrypt.
 * @param schedule       Schedule generated with AES128_load_schedule().
 * @param ciphertext     Output buffer of block_count * 16 bytes.
 */
void AES128_CTR_zero_sch(uint32_t counter_start, size_t block_count,
						 const uint8_t *schedule, uint8_t *ciphertext);

/**
 * Function to free a key schedule.
 *
 * @param schedule       Schedule generated with AES128_load_schedule().
 */
void AES128_free_schedule(uint8_t *schedule);

#endif
