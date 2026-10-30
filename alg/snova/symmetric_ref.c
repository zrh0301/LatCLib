// SPDX-License-Identifier: MIT

/**
 * Reference implementation of the symmetric primitives used by SNOVA
 *
 * The SHAKE/SHA3 hash layer is provided by the shared inc/keccak.h of the pqc
 * library.
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <assert.h>

#include "snova.h"
#include "symmetric.h"

/**
 * Reference SNOVA public key expander
 */

#if defined(AESCTR)

void snova_pk_expander_init(snova_pk_expander_t *instance, const uint8_t *seed, size_t input_bytes) {
	(void)input_bytes;
	memcpy(instance, seed, 16);
}

int AES128_CTR(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen);

void snova_pk_expander(uint8_t *data, size_t num_bytes, snova_pk_expander_t *instance) {
	AES128_CTR(data, num_bytes, (uint8_t *)instance, 16);
}

#else

void snova_pk_expander_init(snova_pk_expander_t *ref_instance, const uint8_t *seed, size_t input_bytes) {
	assert(input_bytes <= 152);

	shake_t *instance = (shake_t *)ref_instance;
	shake128_init(instance);
	sha3_update(instance, seed, input_bytes);
}

void snova_pk_expander(uint8_t *data, size_t num_bytes, snova_pk_expander_t *instance) {
	uint64_t index = 0;
	uint64_t block = 0;
	size_t rate = ((shake_t *)instance)->block_size;

	while (index < num_bytes) {
		shake_t block_instance;
		memcpy(&block_instance, instance, sizeof(shake_t));

		// Turn SHAKE128 into SHAKE128 CTR-XOF
		for (int idx = 0; idx < 8; idx++) {
			// Little endian
			uint8_t block_i = (block >> (8 * idx)) & 0xff;
			sha3_update(&block_instance, &block_i, 1);
		}

		shake_xof_final(&block_instance);
		size_t bytes = num_bytes - index;
		if (bytes > rate) {
			bytes = rate;
		}

		shake_squeeze(&block_instance, data, bytes);

		block++;
		data += bytes;
		index += bytes;
	}
}

#endif

void snova_pk_expand(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, in, inlen);
	snova_pk_expander(out, outlen, &instance);
}
