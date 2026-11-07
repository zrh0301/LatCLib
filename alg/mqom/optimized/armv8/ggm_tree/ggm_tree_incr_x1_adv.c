#include "ggm_tree_incr_x1_adv.h"
#include "ggm_tree_common_ecb.h"

int GGMTree_InitIncrementalExpansion_adv(ggmtree_ctx_adv_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG > 0
	enc_ctx_ecb DECL_VAR(ctx_enc);
#endif

	for (j = GGMTREE_NB_PARALLEL_DERIVATIONS_LOG+1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j >= GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			TweakSalt(salt, tweaked_salt, 2, e, j - 1);
			ret = enc_key_sched_ecb(&ctx->ctx_enc[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt);
			ERR(ret, err);
		}
	}

	// Derive the first level of the tree with "2*GGMTREE_NB_PARALLEL_DERIVATIONS" nodes.
	memcpy(ctx->path[0][0], rseed, MQOM2_PARAM_SEED_SIZE);
	xor_blocks(rseed, delta, ctx->path[0][1]);

#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG > 0
	uint8_t parent_nodes[GGMTREE_NB_PARALLEL_DERIVATIONS][MQOM2_PARAM_SEED_SIZE];
	for (j = 1; j <= GGMTREE_NB_PARALLEL_DERIVATIONS_LOG; j++) {
		memcpy(parent_nodes, ctx->path[0], (1<<j)*MQOM2_PARAM_SEED_SIZE);
		TweakSalt(salt, tweaked_salt, 2, e, j - 1);
		ret = enc_key_sched_ecb(&ctx_enc, tweaked_salt);
		ERR(ret, err);
		DeriveSeeds_ecb(&ctx_enc, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])parent_nodes, ctx->path[0], (1<<j));
	}
#endif

	ctx->salt = salt;
	ctx->e = e;
	ctx->active = 0;

	ret = 0;
err:
#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG > 0
	enc_clean_ctx_ecb(&ctx_enc);
#endif
	return ret;
}

int GGMTree_GetNextLeaf_adv(ggmtree_ctx_adv_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	enc_ctx_ecb DECL_VAR(ctx_enc);
	enc_ctx_ecb* ctx_enc_ptr = NULL;

	if (ctx->active) {
		uint32_t new_num_leaf = ctx->num_leaf + 1;
		uint32_t diff = (ctx->num_leaf ^ new_num_leaf) >> (GGMTREE_NB_PARALLEL_DERIVATIONS_LOG+1);
		ctx->num_leaf = new_num_leaf;
		j = MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-1;
		while (diff != 0) {
			j--; diff>>=1;
		}
	} else {
		ctx->num_leaf = 0;
		ctx->active = 1;
		j = 0;
	}
	for (; j < MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-1; j++) {
		if (j + GGMTREE_NB_PARALLEL_DERIVATIONS_LOG < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY) {
			ctx_enc_ptr = &ctx_enc;
			TweakSalt(ctx->salt, tweaked_salt, 2, ctx->e, j + GGMTREE_NB_PARALLEL_DERIVATIONS_LOG);
			ret = enc_key_sched_ecb(ctx_enc_ptr, tweaked_salt);
			ERR(ret, err);
		} else {
			ctx_enc_ptr = &ctx->ctx_enc[j + GGMTREE_NB_PARALLEL_DERIVATIONS_LOG - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
		}
		uint8_t is_right = (ctx->num_leaf >> (MQOM2_PARAM_NB_EVALS_LOG - 1 - j)) & 0x1;
		DeriveSeeds_ecb(ctx_enc_ptr, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])&ctx->path[j][GGMTREE_NB_PARALLEL_DERIVATIONS*is_right], ctx->path[j+1], GGMTREE_NB_PARALLEL_DERIVATIONS);
	}
	memcpy(lseed, ctx->path[MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-1][ctx->num_leaf & (2*GGMTREE_NB_PARALLEL_DERIVATIONS-1)], MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx_ecb(&ctx_enc);
	return ret;
}

int GGMTree_InitIncrementalPartialExpansion_adv(ggmtree_ctx_partial_adv_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG > 0
	enc_ctx_pub_ecb DECL_VAR(ctx_enc);
#endif

	memcpy((uint8_t*) ctx->opening, (uint8_t*) path, sizeof(ctx->opening));

	for (j = GGMTREE_NB_PARALLEL_DERIVATIONS_LOG+1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j >= GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY + 1) {
			TweakSalt(salt, tweaked_salt, 2, e, j - 1);
			ret = enc_key_sched_pub_ecb(&ctx->ctx_enc[j - 1 - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY], tweaked_salt);
			ERR(ret, err);
		}
	}

	// Derive the first level of the tree with "2*GGMTREE_NB_PARALLEL_DERIVATIONS" nodes.
	memset(ctx->path[0][0], 0, 2*MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[0][(i_star>>(MQOM2_PARAM_NB_EVALS_LOG-1)) ^ 0x01], ctx->opening[MQOM2_PARAM_NB_EVALS_LOG-1], MQOM2_PARAM_SEED_SIZE);

#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG > 0
	uint8_t parent_nodes[GGMTREE_NB_PARALLEL_DERIVATIONS][MQOM2_PARAM_SEED_SIZE];
	for (j = 1; j <= GGMTREE_NB_PARALLEL_DERIVATIONS_LOG; j++) {
		memcpy(parent_nodes, ctx->path[0], (1<<j)*MQOM2_PARAM_SEED_SIZE);
		TweakSalt(salt, tweaked_salt, 2, e, j - 1);
		ret = enc_key_sched_pub_ecb(&ctx_enc, tweaked_salt);
		ERR(ret, err);
		DeriveSeeds_pub_ecb(&ctx_enc, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])parent_nodes, ctx->path[0], (1<<j));

		// Correct with opening
		uint32_t hidden_node_idx = (i_star>>(MQOM2_PARAM_NB_EVALS_LOG-1-j));
		memcpy(ctx->path[0][hidden_node_idx ^ 0x01], ctx->opening[MQOM2_PARAM_NB_EVALS_LOG-1-j], MQOM2_PARAM_SEED_SIZE);
		memset(ctx->path[0][hidden_node_idx], 0, MQOM2_PARAM_SEED_SIZE);
	}
#endif

	ctx->salt = salt;
	ctx->e = e;
	ctx->active = 0;
	ctx->i_star = i_star;

	ret = 0;
err:
#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG > 0
	enc_clean_ctx_pub_ecb(&ctx_enc);
#endif
	return ret;
}

int GGMTree_GetNextLeafPartial_adv(ggmtree_ctx_partial_adv_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]) {
	uint32_t j;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	enc_ctx_pub_ecb DECL_VAR(ctx_enc);
	enc_ctx_pub_ecb* ctx_enc_ptr = NULL;

	if (ctx->active) {
		uint32_t new_num_leaf = ctx->num_leaf + 1;
		uint32_t diff = (ctx->num_leaf ^ new_num_leaf) >> (GGMTREE_NB_PARALLEL_DERIVATIONS_LOG+1);
		ctx->num_leaf = new_num_leaf;
		j = MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-1;
		while (diff != 0) {
			j--; diff>>=1;
		}
	} else {
		ctx->num_leaf = 0;
		ctx->active = 1;
		j = 0;
	}
	uint32_t diff_with_hidden = ctx->i_star^ctx->num_leaf;
	for (; j < MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-1; j++) {
		if (j + GGMTREE_NB_PARALLEL_DERIVATIONS_LOG < GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY) {
			ctx_enc_ptr = &ctx_enc;
			TweakSalt(ctx->salt, tweaked_salt, 2, ctx->e, j + GGMTREE_NB_PARALLEL_DERIVATIONS_LOG);
			ret = enc_key_sched_pub_ecb(ctx_enc_ptr, tweaked_salt);
			ERR(ret, err);
		} else {
			ctx_enc_ptr = &ctx->ctx_enc[j + GGMTREE_NB_PARALLEL_DERIVATIONS_LOG - GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY];
		}
		uint8_t is_right = (ctx->num_leaf >> (MQOM2_PARAM_NB_EVALS_LOG - 1 - j)) & 0x1;
		DeriveSeeds_pub_ecb(ctx_enc_ptr, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])&ctx->path[j][GGMTREE_NB_PARALLEL_DERIVATIONS*is_right], ctx->path[j+1], GGMTREE_NB_PARALLEL_DERIVATIONS);

		// Correct with opening
		if((diff_with_hidden >> (MQOM2_PARAM_NB_EVALS_LOG-1-j)) == 0) {
			uint32_t hidden_node_idx = ((ctx->i_star>>(MQOM2_PARAM_NB_EVALS_LOG-2-j-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG)) & (2*GGMTREE_NB_PARALLEL_DERIVATIONS-1));
			memcpy(ctx->path[j+1][hidden_node_idx ^ 0x01], ctx->opening[MQOM2_PARAM_NB_EVALS_LOG-2-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-j], MQOM2_PARAM_SEED_SIZE);
			memset(ctx->path[j+1][hidden_node_idx], 0, MQOM2_PARAM_SEED_SIZE);
		}
	}
	memcpy(lseed, ctx->path[MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG-1][ctx->num_leaf & (2*GGMTREE_NB_PARALLEL_DERIVATIONS-1)], MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx_pub_ecb(&ctx_enc);
	return ret;
}

