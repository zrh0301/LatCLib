#ifndef __GGM_TREE_COMMON_H__
#define __GGM_TREE_COMMON_H__

#include "enc.h"

#if !defined(GGMTREE_NB_ENC_CTX_IN_MEMORY)
/* Default to 1 */
#define GGMTREE_NB_ENC_CTX_IN_MEMORY 1
#else
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > (MQOM2_PARAM_NB_EVALS_LOG-2)
#error "GGMTREE_NB_ENC_CTX_IN_MEMORY should be smaller than (or equal to) MQOM2_PARAM_NB_EVALS_LOG-2"
#endif
#endif
#define GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY ((MQOM2_PARAM_NB_EVALS_LOG-1)-GGMTREE_NB_ENC_CTX_IN_MEMORY)


/* SeedDerive variants
 * NOTE: we factorize the key schedule, the tweaked salt is inside the encryption context */
static inline void SeedDerive(enc_ctx *ctx, const uint8_t seed[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed[MQOM2_PARAM_SEED_SIZE]) {
	if (seed) {
		uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed, linortho_seed);
		/* Encrypt the LinOrtho seed with the tweaked salt */
		enc_encrypt(ctx, seed, new_seed);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed, linortho_seed, new_seed);
	}

	return;
}
static inline void SeedDerive_pub(enc_ctx_pub *ctx, const uint8_t seed[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed[MQOM2_PARAM_SEED_SIZE]) {
	if (seed) {
		uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed, linortho_seed);
		/* Encrypt the LinOrtho seed with the tweaked salt */
		enc_encrypt_pub(ctx, seed, new_seed);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed, linortho_seed, new_seed);
	}

	return;
}

static inline void SeedDerive_x2(enc_ctx *ctx1, enc_ctx *ctx2, const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE]) {
	if (seed1 && seed2) {
		uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed1, linortho_seed1);
		LinOrtho(seed2, linortho_seed2);
		/* Encrypt the LinOrtho seed with the tweaked salt */
		enc_encrypt_x2(ctx1, ctx2, seed1, seed2, new_seed1, new_seed2);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed1, linortho_seed1, new_seed1);
		xor_blocks(new_seed2, linortho_seed2, new_seed2);
	} else {
		SeedDerive(ctx1, seed1, new_seed1);
		SeedDerive(ctx2, seed2, new_seed2);
	}

	return;
}
static inline void SeedDerive_x2_pub(enc_ctx_pub *ctx1, enc_ctx_pub *ctx2, const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE]) {
	if (seed1 && seed2) {
		uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed1, linortho_seed1);
		LinOrtho(seed2, linortho_seed2);
		/* Encrypt the LinOrtho seed with the tweaked salt */
		enc_encrypt_x2_pub(ctx1, ctx2, seed1, seed2, new_seed1, new_seed2);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed1, linortho_seed1, new_seed1);
		xor_blocks(new_seed2, linortho_seed2, new_seed2);
	} else {
		SeedDerive_pub(ctx1, seed1, new_seed1);
		SeedDerive_pub(ctx2, seed2, new_seed2);
	}

	return;
}

static inline void SeedDerive_x2_x2(enc_ctx_x2 *ctx, const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE]) {
	uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
	uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
	LinOrtho(seed1, linortho_seed1);
	LinOrtho(seed2, linortho_seed2);
	/* Encrypt the LinOrtho seed with the tweaked salt */
	enc_encrypt_x2_x2(ctx, seed1, seed2, new_seed1, new_seed2);
	/* Xor with LinOrtho seed */
	xor_blocks(new_seed1, linortho_seed1, new_seed1);
	xor_blocks(new_seed2, linortho_seed2, new_seed2);
	return;
}

static inline void SeedDerive_x4(enc_ctx *ctx1, enc_ctx *ctx2, enc_ctx *ctx3, enc_ctx *ctx4,
                                 const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                 const uint8_t seed3[MQOM2_PARAM_SEED_SIZE], const uint8_t seed4[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed3[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed4[MQOM2_PARAM_SEED_SIZE]) {
	if (seed1 && seed2 && seed3 && seed4) {
		uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed3[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed4[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed1, linortho_seed1);
		LinOrtho(seed2, linortho_seed2);
		LinOrtho(seed3, linortho_seed3);
		LinOrtho(seed4, linortho_seed4);
		/* Encrypt the seed with the tweaked salt */
		enc_encrypt_x4(ctx1, ctx2, ctx3, ctx4, seed1, seed2, seed3, seed4, new_seed1, new_seed2, new_seed3, new_seed4);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed1, linortho_seed1, new_seed1);
		xor_blocks(new_seed2, linortho_seed2, new_seed2);
		xor_blocks(new_seed3, linortho_seed3, new_seed3);
		xor_blocks(new_seed4, linortho_seed4, new_seed4);
	} else {
		SeedDerive_x2(ctx1, ctx2, seed1, seed2, new_seed1, new_seed2);
		SeedDerive_x2(ctx3, ctx4, seed3, seed4, new_seed3, new_seed4);
	}
	return;
}
static inline void SeedDerive_x4_pub(enc_ctx_pub *ctx1, enc_ctx_pub *ctx2, enc_ctx_pub *ctx3, enc_ctx_pub *ctx4,
                                     const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                     const uint8_t seed3[MQOM2_PARAM_SEED_SIZE], const uint8_t seed4[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed3[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed4[MQOM2_PARAM_SEED_SIZE]) {
	if (seed1 && seed2 && seed3 && seed4) {
		uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed3[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed4[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed1, linortho_seed1);
		LinOrtho(seed2, linortho_seed2);
		LinOrtho(seed3, linortho_seed3);
		LinOrtho(seed4, linortho_seed4);
		/* Encrypt the seed with the tweaked salt */
		enc_encrypt_x4_pub(ctx1, ctx2, ctx3, ctx4, seed1, seed2, seed3, seed4, new_seed1, new_seed2, new_seed3, new_seed4);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed1, linortho_seed1, new_seed1);
		xor_blocks(new_seed2, linortho_seed2, new_seed2);
		xor_blocks(new_seed3, linortho_seed3, new_seed3);
		xor_blocks(new_seed4, linortho_seed4, new_seed4);
	} else {
		SeedDerive_x2_pub(ctx1, ctx2, seed1, seed2, new_seed1, new_seed2);
		SeedDerive_x2_pub(ctx3, ctx4, seed3, seed4, new_seed3, new_seed4);
	}
	return;
}
static inline void SeedDerive_x4_x4(enc_ctx_x4 *ctx,
                                 const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                 const uint8_t seed3[MQOM2_PARAM_SEED_SIZE], const uint8_t seed4[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed3[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed4[MQOM2_PARAM_SEED_SIZE]) {
	uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
	uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
	uint8_t linortho_seed3[MQOM2_PARAM_SEED_SIZE];
	uint8_t linortho_seed4[MQOM2_PARAM_SEED_SIZE];
	LinOrtho(seed1, linortho_seed1);
	LinOrtho(seed2, linortho_seed2);
	LinOrtho(seed3, linortho_seed3);
	LinOrtho(seed4, linortho_seed4);
	/* Encrypt the seed with the tweaked salt */
	enc_encrypt_x4_x4(ctx, seed1, seed2, seed3, seed4, new_seed1, new_seed2, new_seed3, new_seed4);
	/* Xor with LinOrtho seed */
	xor_blocks(new_seed1, linortho_seed1, new_seed1);
	xor_blocks(new_seed2, linortho_seed2, new_seed2);
	xor_blocks(new_seed3, linortho_seed3, new_seed3);
	xor_blocks(new_seed4, linortho_seed4, new_seed4);
	return;
}

static inline void SeedDerive_x8(enc_ctx *ctx1, enc_ctx *ctx2, enc_ctx *ctx3, enc_ctx *ctx4, enc_ctx *ctx5, enc_ctx *ctx6, enc_ctx *ctx7, enc_ctx *ctx8,
                                 const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                 const uint8_t seed3[MQOM2_PARAM_SEED_SIZE], const uint8_t seed4[MQOM2_PARAM_SEED_SIZE],
                                 const uint8_t seed5[MQOM2_PARAM_SEED_SIZE], const uint8_t seed6[MQOM2_PARAM_SEED_SIZE],
                                 const uint8_t seed7[MQOM2_PARAM_SEED_SIZE], const uint8_t seed8[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed3[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed4[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed5[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed6[MQOM2_PARAM_SEED_SIZE],
                                 uint8_t new_seed7[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed8[MQOM2_PARAM_SEED_SIZE])

{
	if (seed1 && seed2 && seed3 && seed4 && seed5 && seed6 && seed7 && seed8) {
		uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed3[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed4[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed5[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed6[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed7[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed8[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed1, linortho_seed1);
		LinOrtho(seed2, linortho_seed2);
		LinOrtho(seed3, linortho_seed3);
		LinOrtho(seed4, linortho_seed4);
		LinOrtho(seed5, linortho_seed5);
		LinOrtho(seed6, linortho_seed6);
		LinOrtho(seed7, linortho_seed7);
		LinOrtho(seed8, linortho_seed8);
		/* Encrypt the seed with the tweaked salt */
		enc_encrypt_x8(ctx1, ctx2, ctx3, ctx4, ctx5, ctx6, ctx7, ctx8,
		               seed1, seed2, seed3, seed4, seed5, seed6, seed7, seed8,
		               new_seed1, new_seed2, new_seed3, new_seed4, new_seed5, new_seed6, new_seed7, new_seed8);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed1, linortho_seed1, new_seed1);
		xor_blocks(new_seed2, linortho_seed2, new_seed2);
		xor_blocks(new_seed3, linortho_seed3, new_seed3);
		xor_blocks(new_seed4, linortho_seed4, new_seed4);
		xor_blocks(new_seed5, linortho_seed5, new_seed5);
		xor_blocks(new_seed6, linortho_seed6, new_seed6);
		xor_blocks(new_seed7, linortho_seed7, new_seed7);
		xor_blocks(new_seed8, linortho_seed8, new_seed8);
	} else {
		SeedDerive_x4(ctx1, ctx2, ctx3, ctx4, seed1, seed2, seed3, seed4, new_seed1, new_seed2, new_seed3, new_seed4);
		SeedDerive_x4(ctx5, ctx6, ctx7, ctx8, seed5, seed6, seed7, seed8, new_seed5, new_seed6, new_seed7, new_seed8);
	}
	return;
}
static inline void SeedDerive_x8_pub(enc_ctx_pub *ctx1, enc_ctx_pub *ctx2, enc_ctx_pub *ctx3, enc_ctx_pub *ctx4, enc_ctx_pub *ctx5, enc_ctx_pub *ctx6, enc_ctx_pub *ctx7, enc_ctx_pub *ctx8,
                                     const uint8_t seed1[MQOM2_PARAM_SEED_SIZE], const uint8_t seed2[MQOM2_PARAM_SEED_SIZE],
                                     const uint8_t seed3[MQOM2_PARAM_SEED_SIZE], const uint8_t seed4[MQOM2_PARAM_SEED_SIZE],
                                     const uint8_t seed5[MQOM2_PARAM_SEED_SIZE], const uint8_t seed6[MQOM2_PARAM_SEED_SIZE],
                                     const uint8_t seed7[MQOM2_PARAM_SEED_SIZE], const uint8_t seed8[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed1[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed2[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed3[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed4[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed5[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed6[MQOM2_PARAM_SEED_SIZE],
                                     uint8_t new_seed7[MQOM2_PARAM_SEED_SIZE], uint8_t new_seed8[MQOM2_PARAM_SEED_SIZE])

{
	if (seed1 && seed2 && seed3 && seed4 && seed5 && seed6 && seed7 && seed8) {
		uint8_t linortho_seed1[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed2[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed3[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed4[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed5[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed6[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed7[MQOM2_PARAM_SEED_SIZE];
		uint8_t linortho_seed8[MQOM2_PARAM_SEED_SIZE];
		LinOrtho(seed1, linortho_seed1);
		LinOrtho(seed2, linortho_seed2);
		LinOrtho(seed3, linortho_seed3);
		LinOrtho(seed4, linortho_seed4);
		LinOrtho(seed5, linortho_seed5);
		LinOrtho(seed6, linortho_seed6);
		LinOrtho(seed7, linortho_seed7);
		LinOrtho(seed8, linortho_seed8);
		/* Encrypt the seed with the tweaked salt */
		enc_encrypt_x8_pub(ctx1, ctx2, ctx3, ctx4, ctx5, ctx6, ctx7, ctx8,
		                   seed1, seed2, seed3, seed4, seed5, seed6, seed7, seed8,
		                   new_seed1, new_seed2, new_seed3, new_seed4, new_seed5, new_seed6, new_seed7, new_seed8);
		/* Xor with LinOrtho seed */
		xor_blocks(new_seed1, linortho_seed1, new_seed1);
		xor_blocks(new_seed2, linortho_seed2, new_seed2);
		xor_blocks(new_seed3, linortho_seed3, new_seed3);
		xor_blocks(new_seed4, linortho_seed4, new_seed4);
		xor_blocks(new_seed5, linortho_seed5, new_seed5);
		xor_blocks(new_seed6, linortho_seed6, new_seed6);
		xor_blocks(new_seed7, linortho_seed7, new_seed7);
		xor_blocks(new_seed8, linortho_seed8, new_seed8);
	} else {
		SeedDerive_x4_pub(ctx1, ctx2, ctx3, ctx4, seed1, seed2, seed3, seed4, new_seed1, new_seed2, new_seed3, new_seed4);
		SeedDerive_x4_pub(ctx5, ctx6, ctx7, ctx8, seed5, seed6, seed7, seed8, new_seed5, new_seed6, new_seed7, new_seed8);
	}
	return;
}


#endif /* __GGM_TREE_COMMON_H__ */
