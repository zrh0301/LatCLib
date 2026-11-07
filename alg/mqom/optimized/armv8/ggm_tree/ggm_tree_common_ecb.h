#ifndef __GGM_TREE_COMMON_ECB_H__
#define __GGM_TREE_COMMON_ECB_H__

#include "common.h"
#include "enc.h"

static inline void SeedDerive_ecb(enc_ctx_ecb *ctx, const uint8_t seed[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed[MQOM2_PARAM_SEED_SIZE]) {
    uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
    LinOrtho(seed, linortho_seed);
    /* Encrypt the LinOrtho seed with the tweaked salt */
    enc_encrypt_ecb(ctx, 1, seed, new_seed);
    /* Xor with LinOrtho seed */
    xor_blocks(new_seed, linortho_seed, new_seed);
	return;
}

static inline void SeedDerive_pub_ecb(enc_ctx_pub_ecb *ctx, const uint8_t seed[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed[MQOM2_PARAM_SEED_SIZE]) {
    uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
    LinOrtho(seed, linortho_seed);
    /* Encrypt the LinOrtho seed with the tweaked salt */
    enc_encrypt_pub_ecb(ctx, 1, seed, new_seed);
    /* Xor with LinOrtho seed */
    xor_blocks(new_seed, linortho_seed, new_seed);
	return;
}

static inline void DeriveSeeds_ecb(enc_ctx_ecb *ctx, const uint8_t in_seeds[][MQOM2_PARAM_SEED_SIZE], uint8_t children[][MQOM2_PARAM_SEED_SIZE], uint32_t nb_in_seeds) {
	uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
	uint32_t i;

	/* Encrypt the LinOrtho seed with the tweaked salt */
	enc_encrypt_ecb(ctx, nb_in_seeds, (const uint8_t*) in_seeds, (uint8_t*) children);
	
	/* Interleave all the derived seeds */
	for (i = nb_in_seeds-1; i > 0; i--) {
		memcpy(children[2*i], children[i], MQOM2_PARAM_SEED_SIZE);
	}

	for (i = 0; i < nb_in_seeds; i++) {
		/* Xor with LinOrtho seed */
		LinOrtho(in_seeds[i], linortho_seed);
		xor_blocks(children[2*i], linortho_seed, children[2*i]);

		/* Derive the sibling node */
		xor_blocks(children[2*i], in_seeds[i], children[2*i+1]);
	}
	return;
}

static inline void DeriveSeeds_pub_ecb(enc_ctx_pub_ecb *ctx, const uint8_t in_seeds[][MQOM2_PARAM_SEED_SIZE], uint8_t children[][MQOM2_PARAM_SEED_SIZE], uint32_t nb_in_seeds) {
	uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
	uint32_t i;

	/* Encrypt the LinOrtho seed with the tweaked salt */
	enc_encrypt_pub_ecb(ctx, nb_in_seeds, (const uint8_t*) in_seeds, (uint8_t*) children);
	
	/* Interleave all the derived seeds */
	for (i = nb_in_seeds-1; i > 0; i--) {
		memcpy(children[2*i], children[i], MQOM2_PARAM_SEED_SIZE);
	}

	for (i = 0; i < nb_in_seeds; i++) {
		/* Xor with LinOrtho seed */
		LinOrtho(in_seeds[i], linortho_seed);
		xor_blocks(children[2*i], linortho_seed, children[2*i]);

		/* Derive the sibling node */
		xor_blocks(children[2*i], in_seeds[i], children[2*i+1]);
	}
	return;
}

#endif /* __GGM_TREE_COMMON_ECB_H__ */
