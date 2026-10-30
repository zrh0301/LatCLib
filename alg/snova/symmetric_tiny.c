// SPDX-License-Identifier: MIT

/**
 * Symmetric primitives used by SNOVA
 *
 * The SHAKE/SHA3 hash layer is provided by the shared inc/keccak.h of the pqc
 * library; this file only contains the small-memory SNOVA public key expander.
 * The optimized implementation of snova_pk_expander can be tested against the reference version by generating the KAT files.
 *
 * Copyright (c) 2026 SNOVA TEAM
 */

#include "symmetric_tiny.h"

#include <stdalign.h>
#include <string.h>

#include "snova.h"

/**
 * SNOVA public key expander
 */

#define get_nibble(st, idx) (((idx) & 1) ? st[(idx) / 2] >> 4 : st[(idx) / 2] & 0xf)

#if defined(AESCTR)

void snova_pk_expand(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
	AES128_CTR(out, outlen, in, inlen);
}

#define NUM_BYTES (12 * 16)
typedef struct {
	uint8_t key[16];
	uint8_t states[NUM_BYTES];
	uint32_t block_i;

	uint32_t index;
	uint32_t last_idx;
	int32_t bytes_left;
} snova_aes_impl_t;

_Static_assert(sizeof(snova_aes_impl_t) <= sizeof(snova_pk_expander_t), "snova_aes_impl_t size error");

void snova_pk_expander_init(snova_pk_expander_t *arg, const uint8_t *seed, size_t input_bytes) {
	snova_aes_impl_t *instance = (snova_aes_impl_t *)arg;

	(void)input_bytes;
	memset(instance, 0, sizeof(snova_aes_impl_t));
	memcpy(instance->key, seed, 16);
}

static void snova_aes_expander_block(snova_aes_impl_t *instance) {
	uint8_t in[NUM_BYTES] = {0};

	for (int i = 0; i < NUM_BYTES / 16; i++) {
		for (int j = 0; j < 4; j++) {
			in[i * 16 + 15 - j] = (instance->block_i >> (8 * j)) & 0xff;
		}
		instance->block_i++;
	}

	AES128_ECB(instance->key, in, instance->states, NUM_BYTES);
}

void snova_pk_expander_squeeze(uint16_t *data, size_t num_gf, snova_pk_expander_t *arg) {
	snova_aes_impl_t *instance = (snova_aes_impl_t *)arg;
	uint16_t *data8 = data;
	instance->last_idx += num_gf;

	if (instance->bytes_left > 0) {
		uint8_t *state8 = (uint8_t *)instance->states;
		int state_idx = 2 * NUM_BYTES - instance->bytes_left;

		if (instance->bytes_left >= (int64_t)num_gf) {
			for (size_t i1 = 0; i1 < num_gf; i1++) {
				data8[i1] = get_nibble(state8, state_idx + i1);
			}
			instance->bytes_left -= num_gf;
			instance->index += num_gf;
			return;
		}

		for (int32_t i1 = 0; i1 < instance->bytes_left; i1++) {
			data8[i1] = get_nibble(state8, state_idx + i1);
		}
		instance->index += instance->bytes_left;
		data8 += instance->bytes_left;
	}

	while (instance->index < instance->last_idx) {
		snova_aes_expander_block(instance);
		// Convert to GF16
		uint8_t *state8 = (uint8_t *)instance->states;

		size_t bytes = instance->last_idx - instance->index;
		if (bytes > 2 * NUM_BYTES) {
			bytes = 2 * NUM_BYTES;
		} else {
			instance->bytes_left = 2 * NUM_BYTES - bytes;
		}

		for (size_t i1 = 0; i1 < bytes; i1++) {
			data8[i1] = get_nibble(state8, i1);
		}
		instance->index += bytes;
		data8 += bytes;
	}
}

void snova_pk_expander_goto(snova_pk_expander_t *arg, size_t index) {
	snova_aes_impl_t *instance = (snova_aes_impl_t *)arg;

	instance->block_i = index / 32;
	snova_aes_expander_block(instance);

	instance->last_idx = index;
	instance->bytes_left = 2 * NUM_BYTES - (index % 32);
	instance->index = index;
}

void snova_pk_expander_free(snova_pk_expander_t *arg) {
	(void)arg;
}

#else

#if __AVX512F__
#define PARALLELISM 8
#undef ROL
#include "keccak_avx512.h"
#elif __AVX2__
#define PARALLELISM 4
#include "keccak_avx2.h"
#else
#define PARALLELISM 1
#include "keccak_opt64.h"
#endif

void snova_pk_expand(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
	uint64_t prepared_state[25 * PARALLELISM];
	uint64_t states[50 * PARALLELISM];
	uint64_t rate;
	uint64_t block;
	uint64_t index;
	uint64_t input_bytes;
	uint64_t last_idx;

	uint64_t keccak_instance[25] = {0};
	uint8_t *prepared_state8 = (uint8_t *)prepared_state;

	block = 0;
	index = 0;
	last_idx = 0;
	input_bytes = inlen;
	rate = 168;

	// Align to uint64_t
	memcpy(&keccak_instance[0], in, inlen);

	for (int idx = 0; idx < PARALLELISM; idx++) {
		for (int idx2 = 0; idx2 < 25; idx2++) {
			prepared_state[idx + idx2 * PARALLELISM] = keccak_instance[idx2];
		}
		// SHAKE padding. Use the (uint8_t *)prepared_state8 here
		prepared_state8[idx * 8 + (inlen + 8) * PARALLELISM] ^= 0x1F;
		prepared_state8[idx * 8 + (rate - 8) * PARALLELISM + 7] ^= 0x80;
	}

	uint8_t *data8 = out;
	last_idx += outlen;

	alignas(PARALLELISM * 8) uint64_t buffer[25 * PARALLELISM];
	while (index < last_idx) {
		memcpy(buffer, prepared_state, PARALLELISM * 200);
		for (int idx = 0; idx < PARALLELISM; idx++) {
#if PARALLELISM == 1
			uint8_t *states8 = (uint8_t *)buffer;
			for (int iend = 0; iend < 8; iend++) {
				uint8_t block_i = (block >> (8 * iend)) & 0xff;
				states8[input_bytes * PARALLELISM + idx * 8 + iend] ^= block_i;
			}
#else
			buffer[input_bytes * PARALLELISM / 8 + idx] ^= block;
#endif
			block++;
		}

#if PARALLELISM == 1
		KeccakF1600_StatePermute((void *)buffer);
#elif PARALLELISM == 4
		KeccakP1600times4_PermuteAll_24rounds((void *)buffer);
#elif PARALLELISM == 8
		KeccakP1600times8_PermuteAll_24rounds((void *)buffer);
#endif

		for (size_t idx = 0; idx < PARALLELISM; idx++) {
			for (size_t idx2 = 0; idx2 < (rate / 8); idx2++) {
				states[idx * (rate / 8) + idx2] = buffer[idx + PARALLELISM * idx2];
			}
		}

		size_t bytes = last_idx - index;
		if (bytes > rate * PARALLELISM) {
			bytes = rate * PARALLELISM;
		}

		memcpy(data8, states, bytes);
		index += bytes;
		data8 += bytes;
	}
}

#define STREAM_PAR 1
#define NUM_BYTES (168 * STREAM_PAR)

typedef struct {
	uint64_t states[25 * STREAM_PAR];
	uint8_t seed[SEED_LENGTH_PUBLIC];
	uint32_t block;
	uint32_t index;
	uint32_t last_idx;
	int32_t bytes_left;
} snova_pkx_impl_t;

_Static_assert(sizeof(snova_pkx_impl_t) == sizeof(snova_pk_expander_t), "snova_pk_expander_t size error");

void snova_pk_expander_init(snova_pk_expander_t *arg, const uint8_t *seed, size_t input_bytes) {
	snova_pkx_impl_t *instance = (snova_pkx_impl_t *)arg;

	instance->block = 0;
	instance->index = 0;
	instance->last_idx = 0;
	instance->bytes_left = 0;
	(void)input_bytes;

	memcpy(instance->seed, seed, SEED_LENGTH_PUBLIC);
}

static void snova_pk_expand_gf_block(snova_pkx_impl_t *instance) {
	alignas(STREAM_PAR * 8) uint64_t buffer[25 * STREAM_PAR];

	// Align to uint64_t
	uint8_t *prepared_state8 = (uint8_t *)buffer;
	uint64_t keccak_instance[25] = {0};
	memcpy(&keccak_instance[0], instance->seed, SEED_LENGTH_PUBLIC);

	for (int idx = 0; idx < STREAM_PAR; idx++) {
		for (int idx2 = 0; idx2 < 25; idx2++) {
			buffer[idx + idx2 * STREAM_PAR] = keccak_instance[idx2];
		}
		// SHAKE padding. Use the (uint8_t *)prepared_state8 here
		prepared_state8[idx * 8 + (SEED_LENGTH_PUBLIC + 8) * STREAM_PAR] ^= 0x1F;
		prepared_state8[idx * 8 + (168 - 8) * STREAM_PAR + 7] ^= 0x80;
	}

	for (int idx = 0; idx < STREAM_PAR; idx++) {
#if STREAM_PAR == 1
		uint8_t *states8 = (uint8_t *)buffer;
		for (int iend = 0; iend < 8; iend++) {
			uint8_t block_i = (instance->block >> (8 * iend)) & 0xff;
			states8[SEED_LENGTH_PUBLIC * STREAM_PAR + idx * 8 + iend] ^= block_i;
		}
#else
		buffer[SEED_LENGTH_PUBLIC * STREAM_PAR / 8 + idx] ^= instance->block;
#endif
		instance->block++;
	}

#if STREAM_PAR == 1
	KeccakF1600_StatePermute((void *)buffer);
#elif STREAM_PAR == 4
	KeccakP1600times4_PermuteAll_24rounds((void *)buffer);
#elif STREAM_PAR == 8
	KeccakP1600times8_PermuteAll_24rounds((void *)buffer);
#endif

	for (size_t idx = 0; idx < STREAM_PAR; idx++) {
		for (size_t idx2 = 0; idx2 < (168 / 8); idx2++) {
			instance->states[idx * (168 / 8) + idx2] = buffer[idx + STREAM_PAR * idx2];
		}
	}
}

void snova_pk_expander_squeeze(uint16_t *data, size_t num_gf, snova_pk_expander_t *arg) {
	snova_pkx_impl_t *instance = (snova_pkx_impl_t *)arg;
	uint16_t *data8 = data;
	instance->last_idx += num_gf;

	if (instance->bytes_left > 0) {
		uint8_t *state8 = (uint8_t *)instance->states;
		int state_idx = 2 * NUM_BYTES - instance->bytes_left;

		if (instance->bytes_left >= (int64_t)num_gf) {
			for (size_t i1 = 0; i1 < num_gf; i1++) {
				data8[i1] = get_nibble(state8, state_idx + i1);
			}
			instance->bytes_left -= num_gf;
			instance->index += num_gf;
			return;
		}

		for (int32_t i1 = 0; i1 < instance->bytes_left; i1++) {
			data8[i1] = get_nibble(state8, state_idx + i1);
		}
		instance->index += instance->bytes_left;
		data8 += instance->bytes_left;
	}

	while (instance->index < instance->last_idx) {
		snova_pk_expand_gf_block(instance);
		// Convert to GF16
		uint8_t *state8 = (uint8_t *)instance->states;

		size_t bytes = instance->last_idx - instance->index;
		if (bytes > 2 * NUM_BYTES) {
			bytes = 2 * NUM_BYTES;
		} else {
			instance->bytes_left = 2 * NUM_BYTES - bytes;
		}

		for (size_t i1 = 0; i1 < bytes; i1++) {
			data8[i1] = get_nibble(state8, i1);
		}
		instance->index += bytes;
		data8 += bytes;
	}
}

void snova_pk_expander_goto(snova_pk_expander_t *arg, size_t index) {
	snova_pkx_impl_t *instance = (snova_pkx_impl_t *)arg;

	instance->block = index / (2 * 168);
	snova_pk_expand_gf_block(instance);

	instance->last_idx = index;
	instance->bytes_left = 2 * STREAM_PAR * 168 - (index % (2 * 168));
	instance->index = index;
}

void snova_pk_expander_free(snova_pk_expander_t *instance) {
	(void)instance;
}
#endif
