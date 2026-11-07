#include "symmetric.h"

#include <stdalign.h>
#include <string.h>

#include "snova.h"

#ifdef USE_OPENSSL

#include <openssl/err.h>
#include <openssl/evp.h>

void shake128_init(shake_t* instance) {
	*(EVP_MD_CTX**)instance = EVP_MD_CTX_new();
	EVP_DigestInit_ex(*(EVP_MD_CTX**)instance, EVP_shake128(), NULL);
}

void shake256_init(shake_t* instance) {
	*(EVP_MD_CTX**)instance = EVP_MD_CTX_new();
	EVP_DigestInit_ex(*(EVP_MD_CTX**)instance, EVP_shake256(), NULL);
}

void shake_absorb(shake_t* instance, const uint8_t* in, size_t inlen) {
	EVP_DigestUpdate(*(EVP_MD_CTX**)instance, in, inlen);
}

void shake_finalize(shake_t* instance) {
	(void)instance;
}

void shake_squeeze(uint8_t* out, size_t outlen, shake_t* instance) {
	EVP_DigestFinalXOF(*(EVP_MD_CTX**)instance, out, outlen);
}

void shake_squeeze_keep(uint8_t* out, size_t outlen, shake_t* instance) {
	EVP_DigestSqueeze(*(EVP_MD_CTX**)instance, out, outlen);
}

void shake_release(shake_t* instance) {
	EVP_MD_CTX_free(*(EVP_MD_CTX**)instance);
}

#else

#include "keccak_opt64.h"

void shake128_init(shake_t* instance) {
	memset(instance, 0, sizeof(shake_t));
	instance->state[26] = 168;
}

void shake256_init(shake_t* instance) {
	memset(instance, 0, sizeof(shake_t));
	instance->state[26] = 136;
}

void shake_absorb(shake_t* instance, const uint8_t* in, size_t inlen) {
	keccak_inc_absorb(instance->state, instance->state[26], in, inlen);
}

void shake_finalize(shake_t* instance) {
	keccak_inc_finalize(instance->state, instance->state[26], 0x1F);
}

void shake_squeeze(uint8_t* out, size_t outlen, shake_t* instance) {
	keccak_inc_squeeze(out, outlen, instance->state, instance->state[26]);
}

void shake_squeeze_keep(uint8_t* out, size_t outlen, shake_t* instance) {
	keccak_inc_squeeze(out, outlen, instance->state, instance->state[26]);
}

void shake_release(shake_t* instance) {
	(void)instance;
}
#endif

void shake256(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen) {
	shake_t instance;
	shake256_init(&instance);
	shake_absorb(&instance, in, inlen);
	shake_finalize(&instance);
	shake_squeeze(out, outlen, &instance);
}

typedef struct {
#if __AVX512F__
	uint64_t state[406];
#elif __AVX2__
	uint64_t state[206];
#else
	uint64_t state[56];
#endif
} snova_pk_expander_t;

#if defined(AESCTR)

void snova_pk_expander_init(snova_pk_expander_t* instance, const uint8_t* seed, size_t input_bytes) {
	(void)input_bytes;
	memcpy(instance, seed, 16);
}

#define AES128_CTR SNOVA_NAMESPACE(AES128_CTR)
int AES128_CTR(unsigned char* output, size_t outputByteLen, const unsigned char* input, size_t inputByteLen);

void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* instance) {
	AES128_CTR(data, num_bytes, (uint8_t*)instance, 16);
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

typedef struct {
	uint64_t prepared_state[25 * PARALLELISM];
	uint64_t states[25 * PARALLELISM];
	uint64_t rate;
	uint64_t block;
	uint64_t index;
	uint64_t input_bytes;
	uint64_t last_idx;
	int64_t bytes_left;
} snova_pkx_impl_t;

void snova_pk_expander_init(snova_pk_expander_t* arg, const uint8_t* seed, size_t input_bytes) {
	snova_pkx_impl_t *instance = (snova_pkx_impl_t*)arg;
	uint64_t keccak_instance[25] = {0};
	uint8_t *prepared_state8 = (uint8_t*)instance->prepared_state;

	instance->block = 0;
	instance->index = 0;
	instance->last_idx = 0;
	instance->input_bytes = input_bytes;
	instance->bytes_left = 0;
	instance->rate = 168;

	memcpy(&keccak_instance[0], seed, input_bytes);

	for (int idx = 0; idx < PARALLELISM; idx++) {
		for (int idx2 = 0; idx2 < 25; idx2++) {
			instance->prepared_state[idx + idx2 * PARALLELISM] = keccak_instance[idx2];
		}

		prepared_state8[idx * 8 + (input_bytes + 8) * PARALLELISM] ^= 0x1F;
		prepared_state8[idx * 8 + (instance->rate - 8) * PARALLELISM + 7] ^= 0x80;
	}
}

void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* arg) {
	snova_pkx_impl_t *instance = (snova_pkx_impl_t*)arg;
	uint8_t *data8 = data;
	instance->last_idx += num_bytes;

	if (instance->bytes_left > 0) {
		uint8_t *state8 = (uint8_t*)instance->states;
		state8 += instance->rate * PARALLELISM - instance->bytes_left;

		if (instance->bytes_left > (int64_t)num_bytes) {
			memcpy(data, state8, num_bytes);
			instance->bytes_left -= num_bytes;
			instance->index += num_bytes;
			return;
		}

		memcpy(data8, state8, instance->bytes_left);
		instance->index += instance->bytes_left;
		data8 += instance->bytes_left;
	}

	alignas(PARALLELISM * 8) uint64_t buffer[25 * PARALLELISM];
	while (instance->index < instance->last_idx) {
		memcpy(buffer, instance->prepared_state, PARALLELISM * 200);
		for (int idx = 0; idx < PARALLELISM; idx++) {
#if PARALLELISM == 1
			uint8_t *states8 = (uint8_t*)buffer;
			for (int iend = 0; iend < 8; iend++) {
				uint8_t block_i = (instance->block >> (8 * iend)) & 0xff;
				states8[instance->input_bytes * PARALLELISM + idx * 8 + iend] ^= block_i;
			}
#else
			buffer[instance->input_bytes * PARALLELISM / 8 + idx] ^= instance->block;
#endif
			instance->block++;
		}

#if PARALLELISM == 1
		KeccakF1600_StatePermute((void*)buffer);
#elif PARALLELISM == 4
		KeccakP1600times4_PermuteAll_24rounds((void*)buffer);
#elif PARALLELISM == 8
		KeccakP1600times8_PermuteAll_24rounds((void*)buffer);
#endif

		for (size_t idx = 0; idx < PARALLELISM; idx++) {
			for (size_t idx2 = 0; idx2 < (instance->rate / 8); idx2++) {
				instance->states[idx * (instance->rate / 8) + idx2] = buffer[idx + PARALLELISM * idx2];
			}
		}

		size_t bytes = instance->last_idx - instance->index;
		if (bytes > instance->rate * PARALLELISM) {
			bytes = instance->rate * PARALLELISM;
		} else {
			instance->bytes_left = instance->rate * PARALLELISM - bytes;
		}

		memcpy(data8, instance->states, bytes);
		instance->index += bytes;
		data8 += bytes;
	}
}
#endif

void snova_pk_expand(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen) {
	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, in, inlen);
	snova_pk_expander(out, outlen, &instance);
}
