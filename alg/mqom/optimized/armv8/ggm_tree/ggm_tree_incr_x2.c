#include "ggm_tree_incr_x2.h"

int GGMTree_InitIncrementalExpansion_x2(ggmtree_ctx_x2_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[2][MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], const uint32_t e[2]) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[2][MQOM2_PARAM_SEED_SIZE];
	for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		TweakSalt(salt, tweaked_salt[0], 2, e[0], j - 1);
		TweakSalt(salt, tweaked_salt[1], 2, e[1], j - 1);
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			memcpy(ctx->tweaked_salts[0][j - 1], tweaked_salt[0], MQOM2_PARAM_SALT_SIZE);
			memcpy(ctx->tweaked_salts[1][j - 1], tweaked_salt[1], MQOM2_PARAM_SALT_SIZE);
		} else {
			ret = enc_key_sched_x2(&ctx->ctx_enc_x2[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt[0], tweaked_salt[1]);
			ERR(ret, err);
		}
	}
	memcpy(ctx->path[0][0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[0][1], rseed[0], MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[1][0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[1][1], rseed[1], MQOM2_PARAM_SEED_SIZE);
	ctx->active = 0;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeaf_x2(ggmtree_ctx_x2_t* ctx, uint8_t lseed[2][MQOM2_PARAM_SEED_SIZE]) {
	enc_ctx_x2 DECL_VAR(ctx_enc_x2);
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
		xor_blocks(ctx->path[0][j - 1], ctx->path[0][j], ctx->path[0][j]);
		xor_blocks(ctx->path[1][j - 1], ctx->path[1][j], ctx->path[1][j]);
	} else {
		ctx->num_leaf = 0;
		ctx->active = 1;
		j = 1;
	}
	enc_ctx_x2* ctx_enc_ptr = &ctx_enc_x2;
	for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			ret = enc_key_sched_x2(ctx_enc_ptr, ctx->tweaked_salts[0][j-1], ctx->tweaked_salts[1][j-1]);
			ERR(ret, err);
		} else {
			ctx_enc_ptr = &ctx->ctx_enc_x2[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
		}
		SeedDerive_x2_x2(ctx_enc_ptr,
		              ctx->path[0][j], ctx->path[1][j],
		              ctx->path[0][j + 1], ctx->path[1][j + 1]);
	}
	memcpy(lseed[0], ctx->path[0][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
	memcpy(lseed[1], ctx->path[1][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx_x2(&ctx_enc_x2);
	return ret;
}

int GGMTree_InitIncrementalPartialExpansion_x2(ggmtree_ctx_partial_x2_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t (*path[2])[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], const uint32_t e[2], const uint32_t i_star[2]) {
	uint32_t i, j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	for (i = 0; i < 2; i++) {
		for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
			if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
				TweakSalt(salt, ctx->tweaked_salts[i][j - 1], 2, e[i], j - 1);
			} else {
				TweakSalt(salt, tweaked_salt, 2, e[i], j - 1);
				ret = enc_key_sched_pub(&ctx->ctx_enc[i][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt);
				ERR(ret, err);
			}
		}
	}
	memcpy((uint8_t*) ctx->opening[0], (uint8_t*) path[0], sizeof(ctx->opening[0]));
	memcpy((uint8_t*) ctx->opening[1], (uint8_t*) path[1], sizeof(ctx->opening[1]));
	ctx->i_star[0] = i_star[0];
	ctx->i_star[1] = i_star[1];
	ctx->active = 0;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeafPartial_x2(ggmtree_ctx_partial_x2_t* ctx, uint8_t lseed[2][MQOM2_PARAM_SEED_SIZE]) {
	enc_ctx_pub DECL_VAR(ctx_enc[2]);
	uint32_t i, j;
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
	uint32_t diffs[2] = {ctx->num_leaf ^ctx->i_star[0], ctx->num_leaf ^ctx->i_star[1]};
	uint8_t node_map[2][MQOM2_PARAM_NB_EVALS_LOG];
	for (i = 0; i < 2; i++) {
		uint32_t higher = 1;
		if (diffs[i]) {
			while (((diffs[i] >> (MQOM2_PARAM_NB_EVALS_LOG - higher)) & 0x1) == 0) {
				higher++;
			}
		} else {
			higher = MQOM2_PARAM_NB_EVALS_LOG;
		}
		for (uint32_t k = j; k < higher; k++) {
			node_map[i][k] = 0;
		}
		for (uint32_t k = higher; k < MQOM2_PARAM_NB_EVALS_LOG; k++) {
			node_map[i][k] = 1;
		}
		if (j > higher) {
			xor_blocks(ctx->path[i][j - 1], ctx->path[i][j], ctx->path[i][j]);
		} else if (j <= higher && diffs[i]) {
			memcpy(ctx->path[i][higher], ctx->opening[i][MQOM2_PARAM_NB_EVALS_LOG - higher], MQOM2_PARAM_SEED_SIZE);
		}
	}
	for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		enc_ctx_pub* ctx_enc0 = NULL;
		enc_ctx_pub* ctx_enc1 = NULL;
		uint8_t *seed0 = node_map[0][j] ? ctx->path[0][j] : NULL;
		uint8_t *seed1 = node_map[1][j] ? ctx->path[1][j] : NULL;
		uint8_t *new_seed0 = ctx->path[0][j + 1];
		uint8_t *new_seed1 = ctx->path[1][j + 1];
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			ctx_enc0 = &ctx_enc[0];
			ctx_enc1 = &ctx_enc[1];
			ret = enc_key_sched_pub(ctx_enc0, ctx->tweaked_salts[0][j - 1]);
			ERR(ret, err);
			ret = enc_key_sched_pub(ctx_enc1, ctx->tweaked_salts[1][j - 1]);
			ERR(ret, err);
		} else {
			ctx_enc0 = &ctx->ctx_enc[0][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
			ctx_enc1 = &ctx->ctx_enc[1][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
		}
		SeedDerive_x2_pub(ctx_enc0, ctx_enc1,
		                  seed0, seed1,
		                  new_seed0, new_seed1);
	}
	for (i = 0; i < 2; i++) {
		if (diffs[i]) {
			memcpy(lseed[i], ctx->path[i][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
		} else {
			memset(lseed[i], 0, MQOM2_PARAM_SEED_SIZE);
		}
	}

	ret = 0;
err:
	for (j = 0; j < 2; j++) {
		enc_clean_ctx_pub(&ctx_enc[j]);
	}

	return ret;
}
