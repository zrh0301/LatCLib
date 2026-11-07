#include "ggm_tree_incr_x1_base.h"

int GGMTree_InitIncrementalExpansion_base(ggmtree_ctx_base_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			TweakSalt(salt, ctx->tweaked_salts[j - 1], 2, e, j - 1);
		} else {
			TweakSalt(salt, tweaked_salt, 2, e, j - 1);
			ret = enc_key_sched(&ctx->ctx_enc[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt);
			ERR(ret, err);
		}
	}
	memcpy(ctx->path[0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[1], rseed, MQOM2_PARAM_SEED_SIZE);
	ctx->active = 0;
	ctx->e = e;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeaf_base(ggmtree_ctx_base_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]) {
	uint32_t j;
	int ret = -1;

	if (ctx->active) {
		uint32_t new_num_leaf = ctx->num_leaf + 1;
		uint32_t diff = ctx->num_leaf ^ new_num_leaf;
		ctx->num_leaf = new_num_leaf;
		j = 1;
		while (((diff >> (MQOM2_PARAM_NB_EVALS_LOG - j)) & 0x1) == 0) {
			j++;
		}
		xor_blocks(ctx->path[j - 1], ctx->path[j], ctx->path[j]);
	} else {
		ctx->num_leaf = 0;
		ctx->active = 1;
		j = 1;
	}
	for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			enc_ctx DECL_VAR(ctx_enc);
			ret = enc_key_sched(&ctx_enc, ctx->tweaked_salts[j - 1]);
			ERR(ret, err1);
			SeedDerive(&ctx_enc, ctx->path[j], ctx->path[j + 1]);
err1:
			enc_clean_ctx(&ctx_enc);
			ERR(ret, err);
		} else {
			enc_ctx* ctx_enc_precomputed = &ctx->ctx_enc[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
			SeedDerive(ctx_enc_precomputed, ctx->path[j], ctx->path[j + 1]);
		}
	}
	memcpy(lseed, ctx->path[MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	return ret;
}

int GGMTree_InitIncrementalPartialExpansion_base(ggmtree_ctx_partial_base_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			TweakSalt(salt, ctx->tweaked_salts[j - 1], 2, e, j - 1);
		} else {
			TweakSalt(salt, tweaked_salt, 2, e, j - 1);
			ret = enc_key_sched_pub(&ctx->ctx_enc[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt);
			ERR(ret, err);
		}
	}
	memcpy((uint8_t*) ctx->opening, (uint8_t*) path, sizeof(ctx->opening));
	ctx->i_star = i_star;
	ctx->active = 0;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeafPartial_base(ggmtree_ctx_partial_base_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]) {
	uint32_t j;
	int ret = -1;

	if (ctx->active) {
		uint32_t new_num_leaf = ctx->num_leaf + 1;
		uint32_t diff = ctx->num_leaf ^ new_num_leaf;
		ctx->num_leaf = new_num_leaf;
		j = 1;
		while (((diff >> (MQOM2_PARAM_NB_EVALS_LOG - j)) & 0x1) == 0) {
			j++;
		}
	} else {
		ctx->num_leaf = 0;
		ctx->active = 1;
		j = 1;
	}
	uint32_t diff2 = ctx->num_leaf ^ ctx->i_star;
	if (diff2 == 0) {
		memset(lseed, 0, MQOM2_PARAM_SEED_SIZE);
	} else {
		uint32_t higher = 1;
		while (((diff2 >> (MQOM2_PARAM_NB_EVALS_LOG - higher)) & 0x1) == 0) {
			higher++;
		}
		if (j <= higher) {
			memcpy(ctx->path[higher], ctx->opening[MQOM2_PARAM_NB_EVALS_LOG - higher], MQOM2_PARAM_SEED_SIZE);
			j = higher;
		} else { // j > higher
			xor_blocks(ctx->path[j - 1], ctx->path[j], ctx->path[j]);
		}
		for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
			if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
				enc_ctx_pub DECL_VAR(ctx_enc);
				ret = enc_key_sched_pub(&ctx_enc, ctx->tweaked_salts[j - 1]);
				ERR(ret, err1);
				SeedDerive_pub(&ctx_enc, ctx->path[j], ctx->path[j + 1]);
err1:
				enc_clean_ctx_pub(&ctx_enc);
				ERR(ret, err);
			} else {
				enc_ctx_pub* ctx_enc_precomputed = &ctx->ctx_enc[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
				SeedDerive_pub(ctx_enc_precomputed, ctx->path[j], ctx->path[j + 1]);
			}
		}
		memcpy(lseed, ctx->path[MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
	}

	ret = 0;
err:
	return ret;
}
