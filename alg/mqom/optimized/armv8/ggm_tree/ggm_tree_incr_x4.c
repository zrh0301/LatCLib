#include "ggm_tree_incr_x4.h"

int GGMTree_InitIncrementalExpansion_x4(ggmtree_ctx_x4_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[4][MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], const uint32_t e[4]) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[4][MQOM2_PARAM_SEED_SIZE];
	for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		TweakSalt(salt, tweaked_salt[0], 2, e[0], j - 1);
		TweakSalt(salt, tweaked_salt[1], 2, e[1], j - 1);
		TweakSalt(salt, tweaked_salt[2], 2, e[2], j - 1);
		TweakSalt(salt, tweaked_salt[3], 2, e[3], j - 1);
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			memcpy(ctx->tweaked_salts[0][j - 1], tweaked_salt[0], MQOM2_PARAM_SALT_SIZE);
			memcpy(ctx->tweaked_salts[1][j - 1], tweaked_salt[1], MQOM2_PARAM_SALT_SIZE);
			memcpy(ctx->tweaked_salts[2][j - 1], tweaked_salt[2], MQOM2_PARAM_SALT_SIZE);
			memcpy(ctx->tweaked_salts[3][j - 1], tweaked_salt[3], MQOM2_PARAM_SALT_SIZE);
		} else {
			ret = enc_key_sched_x4(&ctx->ctx_enc_x4[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt[0], tweaked_salt[1], tweaked_salt[2], tweaked_salt[3]);
			ERR(ret, err);
		}
	}
	memcpy(ctx->path[0][0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[0][1], rseed[0], MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[1][0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[1][1], rseed[1], MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[2][0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[2][1], rseed[2], MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[3][0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[3][1], rseed[3], MQOM2_PARAM_SEED_SIZE);
	ctx->active = 0;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeaf_x4(ggmtree_ctx_x4_t* ctx, uint8_t lseed[4][MQOM2_PARAM_SEED_SIZE]) {
	enc_ctx_x4 DECL_VAR(ctx_enc_x4);
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
		xor_blocks(ctx->path[2][j - 1], ctx->path[2][j], ctx->path[2][j]);
		xor_blocks(ctx->path[3][j - 1], ctx->path[3][j], ctx->path[3][j]);
	} else {
		ctx->num_leaf = 0;
		ctx->active = 1;
		j = 1;
	}
	enc_ctx_x4* ctx_enc_ptr = &ctx_enc_x4;
	for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			ret = enc_key_sched_x4(ctx_enc_ptr, ctx->tweaked_salts[0][j-1], ctx->tweaked_salts[1][j-1], ctx->tweaked_salts[2][j-1], ctx->tweaked_salts[3][j-1]);
			ERR(ret, err);
		} else {
			ctx_enc_ptr = &ctx->ctx_enc_x4[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
		}
		SeedDerive_x4_x4(ctx_enc_ptr,
		              ctx->path[0][j], ctx->path[1][j], ctx->path[2][j], ctx->path[3][j],
		              ctx->path[0][j + 1], ctx->path[1][j + 1], ctx->path[2][j + 1], ctx->path[3][j + 1]);
	}
	memcpy(lseed[0], ctx->path[0][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
	memcpy(lseed[1], ctx->path[1][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
	memcpy(lseed[2], ctx->path[2][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
	memcpy(lseed[3], ctx->path[3][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx_x4(&ctx_enc_x4);
	return ret;
}

int GGMTree_InitIncrementalPartialExpansion_x4(ggmtree_ctx_partial_x4_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t (*path[4])[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], const uint32_t e[4], const uint32_t i_star[4]) {
	uint32_t i, j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	for (i = 0; i < 4; i++) {
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
	memcpy((uint8_t*) ctx->opening[2], (uint8_t*) path[2], sizeof(ctx->opening[2]));
	memcpy((uint8_t*) ctx->opening[3], (uint8_t*) path[3], sizeof(ctx->opening[3]));
	ctx->i_star[0] = i_star[0];
	ctx->i_star[1] = i_star[1];
	ctx->i_star[2] = i_star[2];
	ctx->i_star[3] = i_star[3];
	ctx->active = 0;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeafPartial_x4(ggmtree_ctx_partial_x4_t* ctx, uint8_t lseed[4][MQOM2_PARAM_SEED_SIZE]) {
	enc_ctx_pub DECL_VAR(ctx_enc[4]);
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
	uint32_t diffs[4] = {ctx->num_leaf ^ctx->i_star[0], ctx->num_leaf ^ctx->i_star[1], ctx->num_leaf ^ctx->i_star[2], ctx->num_leaf ^ctx->i_star[3]};
	uint8_t node_map[4][MQOM2_PARAM_NB_EVALS_LOG];
	for (i = 0; i < 4; i++) {
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
		enc_ctx_pub* ctx_enc2 = NULL;
		enc_ctx_pub* ctx_enc3 = NULL;
		uint8_t *seed0 = node_map[0][j] ? ctx->path[0][j] : NULL;
		uint8_t *seed1 = node_map[1][j] ? ctx->path[1][j] : NULL;
		uint8_t *seed2 = node_map[2][j] ? ctx->path[2][j] : NULL;
		uint8_t *seed3 = node_map[3][j] ? ctx->path[3][j] : NULL;
		uint8_t *new_seed0 = ctx->path[0][j + 1];
		uint8_t *new_seed1 = ctx->path[1][j + 1];
		uint8_t *new_seed2 = ctx->path[2][j + 1];
		uint8_t *new_seed3 = ctx->path[3][j + 1];
		if (j < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			ctx_enc0 = &ctx_enc[0];
			ctx_enc1 = &ctx_enc[1];
			ctx_enc2 = &ctx_enc[2];
			ctx_enc3 = &ctx_enc[3];
			ret = enc_key_sched_pub(ctx_enc0, ctx->tweaked_salts[0][j - 1]);
			ERR(ret, err);
			ret = enc_key_sched_pub(ctx_enc1, ctx->tweaked_salts[1][j - 1]);
			ERR(ret, err);
			ret = enc_key_sched_pub(ctx_enc2, ctx->tweaked_salts[2][j - 1]);
			ERR(ret, err);
			ret = enc_key_sched_pub(ctx_enc3, ctx->tweaked_salts[3][j - 1]);
			ERR(ret, err);
		} else {
			ctx_enc0 = &ctx->ctx_enc[0][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
			ctx_enc1 = &ctx->ctx_enc[1][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
			ctx_enc2 = &ctx->ctx_enc[2][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
			ctx_enc3 = &ctx->ctx_enc[3][j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
		}
		SeedDerive_x4_pub(ctx_enc0, ctx_enc1, ctx_enc2, ctx_enc3,
		                  seed0, seed1, seed2, seed3,
		                  new_seed0, new_seed1, new_seed2, new_seed3);
	}
	for (i = 0; i < 4; i++) {
		if (diffs[i]) {
			memcpy(lseed[i], ctx->path[i][MQOM2_PARAM_NB_EVALS_LOG], MQOM2_PARAM_SEED_SIZE);
		} else {
			memset(lseed[i], 0, MQOM2_PARAM_SEED_SIZE);
		}
	}

	ret = 0;
err:
	for (j = 0; j < 4; j++) {
		enc_clean_ctx_pub(&ctx_enc[j]);
	}

	return ret;
}
