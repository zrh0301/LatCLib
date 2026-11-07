#ifndef __BLC_MEMOPT_X1_SEEDCOMMIT_H__
#define __BLC_MEMOPT_X1_SEEDCOMMIT_H__

#include "xof.h"
#include "enc.h"
#include "blc_memopt_x1.h"

typedef struct {
	xof_context xof_ctx;
	uint32_t e;
	const uint8_t* salt;
#ifdef BLC_SEEDCOMMIT_CACHE
	enc_ctx_ecb DECL_VAR(enc_ctx[2]);
#endif
} seedcommit_sign_ctx_t;

typedef struct {
	xof_context xof_ctx;
	uint32_t e;
	const uint8_t* salt;
	uint32_t hidden_index;
	uint8_t out_ls_com[MQOM2_PARAM_DIGEST_SIZE];
#ifdef BLC_SEEDCOMMIT_CACHE
	enc_ctx_pub_ecb DECL_VAR(enc_ctx[2]);
#endif
} seedcommit_verify_ctx_t;

static inline int init_seedcommit_sign(seedcommit_sign_ctx_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], uint32_t e) {
	int ret = -1;

	__BENCHMARK_START__(BS_BLC_XOF);
	ret = xof_init(&ctx->xof_ctx);
	ERR(ret, err);
	ret = xof_update(&ctx->xof_ctx, (const uint8_t*) "\x06", 1);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_XOF);

	ctx->salt = salt;
	ctx->e = e;

#ifdef BLC_SEEDCOMMIT_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	TweakSalt(ctx->salt, tweaked_salt, 0, ctx->e, 0);
	ret = enc_key_sched_ecb(&ctx->enc_ctx[0], tweaked_salt);
	ERR(ret, err);
	tweaked_salt[0] ^= 0x01;
	ret = enc_key_sched_ecb(&ctx->enc_ctx[1], tweaked_salt);
	ERR(ret, err);
#endif

	ret = 0;
err:
	if(ret){
		xof_clean_ctx(&ctx->xof_ctx);
	}
	return ret;
}

static inline int init_seedcommit_verify(seedcommit_verify_ctx_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], uint32_t e, uint32_t hidden_index, const uint8_t out_ls_com[MQOM2_PARAM_DIGEST_SIZE]) {
	int ret = -1;

	__BENCHMARK_START__(BS_BLC_XOF);
	ret = xof_init(&ctx->xof_ctx);
	ERR(ret, err);
	ret = xof_update(&ctx->xof_ctx, (const uint8_t*) "\x06", 1);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_XOF);

	ctx->salt = salt;
	ctx->e = e;

	memcpy(ctx->out_ls_com, out_ls_com, MQOM2_PARAM_DIGEST_SIZE);
	ctx->hidden_index = hidden_index;

#ifdef BLC_SEEDCOMMIT_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	TweakSalt(ctx->salt, tweaked_salt, 0, ctx->e, 0);
	ret = enc_key_sched_pub_ecb(&ctx->enc_ctx[0], tweaked_salt);
	ERR(ret, err);
	tweaked_salt[0] ^= 0x01;
	ret = enc_key_sched_pub_ecb(&ctx->enc_ctx[1], tweaked_salt);
	ERR(ret, err);
#endif

	ret = 0;
err:
	if(ret){
		xof_clean_ctx(&ctx->xof_ctx);
	}
	return ret;
}

static inline void seedcommit_sign_clean_ctx(seedcommit_sign_ctx_t* ctx) {
	xof_clean_ctx(&ctx->xof_ctx);
#ifdef BLC_SEEDCOMMIT_CACHE
	enc_clean_ctx_ecb(&ctx->enc_ctx[0]);
	enc_clean_ctx_ecb(&ctx->enc_ctx[1]);
#endif
}

static inline void seedcommit_verify_clean_ctx(seedcommit_verify_ctx_t* ctx) {
	xof_clean_ctx(&ctx->xof_ctx);
#ifdef BLC_SEEDCOMMIT_CACHE
	enc_clean_ctx_pub_ecb(&ctx->enc_ctx[0]);
	enc_clean_ctx_pub_ecb(&ctx->enc_ctx[1]);
#endif
}

static inline int SeedCommitThenAbsorb_sign(seedcommit_sign_ctx_t* ctx, const uint8_t lseed[][MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	uint32_t i_;

	uint8_t buffer[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	uint8_t ls_com[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_DIGEST_SIZE];
	uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
#ifndef BLC_SEEDCOMMIT_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	enc_ctx_ecb DECL_VAR(enc_ctx_);
	enc_ctx_ecb* enc_ctx_ptr = &enc_ctx_;
#else
	enc_ctx_ecb* enc_ctx_ptr;
#endif

	__BENCHMARK_START__(BS_BLC_SEED_COMMIT);
#ifndef BLC_SEEDCOMMIT_CACHE
	TweakSalt(ctx->salt, tweaked_salt, 0, ctx->e, 0);
	ret = enc_key_sched_ecb(enc_ctx_ptr, tweaked_salt);
	ERR(ret, err);
#else
	enc_ctx_ptr = &ctx->enc_ctx[0];
#endif
	ret = enc_encrypt_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) buffer);
	ERR(ret, err);
	for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		memcpy(&ls_com[i_][0], buffer[i_], MQOM2_PARAM_SEED_SIZE);
	}

#ifndef BLC_SEEDCOMMIT_CACHE
	tweaked_salt[0] ^= 0x01;
	ret = enc_key_sched_ecb(enc_ctx_ptr, tweaked_salt);
	ERR(ret, err);
#else
	enc_ctx_ptr = &ctx->enc_ctx[1];
#endif
	ret = enc_encrypt_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) buffer);
	ERR(ret, err);
	for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		memcpy(&ls_com[i_][MQOM2_PARAM_SEED_SIZE], buffer[i_], MQOM2_PARAM_SEED_SIZE);
	}

	for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		LinOrtho(lseed[i_], linortho_seed);
		xor_blocks(&ls_com[i_][0], linortho_seed, &ls_com[i_][0]);
		xor_blocks(&ls_com[i_][MQOM2_PARAM_SEED_SIZE], linortho_seed, &ls_com[i_][MQOM2_PARAM_SEED_SIZE]);
	}
	__BENCHMARK_STOP__(BS_BLC_SEED_COMMIT);

	__BENCHMARK_START__(BS_BLC_XOF);
	ret = xof_update(&ctx->xof_ctx, (uint8_t*) ls_com, BLC_NB_LEAF_SEEDS_IN_PARALLEL * MQOM2_PARAM_DIGEST_SIZE);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_XOF);

	ret = 0;
err:
#ifndef BLC_SEEDCOMMIT_CACHE
	enc_clean_ctx_ecb(&enc_ctx_);
#endif
	return ret;
}

static inline int SeedCommitThenAbsorb_verify(seedcommit_verify_ctx_t* ctx, uint32_t i, const uint8_t lseed[][MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	uint32_t i_;

	uint8_t buffer[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	uint8_t ls_com[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_DIGEST_SIZE];
	uint8_t linortho_seed[MQOM2_PARAM_SEED_SIZE];
#ifndef BLC_SEEDCOMMIT_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	enc_ctx_pub_ecb DECL_VAR(enc_ctx_);
	enc_ctx_pub_ecb* enc_ctx_ptr = &enc_ctx_;
#else
	enc_ctx_pub_ecb* enc_ctx_ptr;
#endif

	__BENCHMARK_START__(BS_BLC_SEED_COMMIT);
#ifndef BLC_SEEDCOMMIT_CACHE
	TweakSalt(ctx->salt, tweaked_salt, 0, ctx->e, 0);
	ret = enc_key_sched_pub_ecb(enc_ctx_ptr, tweaked_salt);
	ERR(ret, err);
#else
	enc_ctx_ptr = &ctx->enc_ctx[0];
#endif
	ret = enc_encrypt_pub_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) buffer);
	ERR(ret, err);
	for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		memcpy(&ls_com[i_][0], buffer[i_], MQOM2_PARAM_SEED_SIZE);
	}

#ifndef BLC_SEEDCOMMIT_CACHE
	tweaked_salt[0] ^= 0x01;
	ret = enc_key_sched_pub_ecb(enc_ctx_ptr, tweaked_salt);
	ERR(ret, err);
#else
	enc_ctx_ptr = &ctx->enc_ctx[1];
#endif
	ret = enc_encrypt_pub_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) buffer);
	ERR(ret, err);
	for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		memcpy(&ls_com[i_][MQOM2_PARAM_SEED_SIZE], buffer[i_], MQOM2_PARAM_SEED_SIZE);
	}

	for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		LinOrtho(lseed[i_], linortho_seed);
		xor_blocks(&ls_com[i_][0], linortho_seed, &ls_com[i_][0]);
		xor_blocks(&ls_com[i_][MQOM2_PARAM_SEED_SIZE], linortho_seed, &ls_com[i_][MQOM2_PARAM_SEED_SIZE]);
	}
	__BENCHMARK_STOP__(BS_BLC_SEED_COMMIT);

	if(ctx->hidden_index < i + BLC_NB_LEAF_SEEDS_IN_PARALLEL && ctx->hidden_index >= i) {
		memcpy(ls_com[ctx->hidden_index-i], ctx->out_ls_com, MQOM2_PARAM_DIGEST_SIZE);
	}

	__BENCHMARK_START__(BS_BLC_XOF);
	ret = xof_update(&ctx->xof_ctx, (uint8_t*) ls_com, BLC_NB_LEAF_SEEDS_IN_PARALLEL * MQOM2_PARAM_DIGEST_SIZE);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_XOF);

	ret = 0;
err:
#ifndef BLC_SEEDCOMMIT_CACHE
	enc_clean_ctx_pub_ecb(&enc_ctx_);
#endif
	return ret;
}

#endif /* __BLC_MEMOPT_X1_SEEDCOMMIT_H__ */
