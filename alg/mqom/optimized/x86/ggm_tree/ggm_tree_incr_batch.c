#include "ggm_tree_incr_batch.h"
#include "ggm_tree_common_ecb.h"

static inline uint32_t get_depth_of_common_ancestor(uint32_t leaf1, uint32_t leaf2) {
	uint32_t diff = leaf1 ^ leaf2;
	uint32_t j = MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG+1;
	while (diff != 0) {
		j--; diff>>=1;
	}
	return j;
}

static inline uint32_t increment_num_leaf(uint32_t* active, uint32_t* num_leaf) {
	uint32_t j;
	if (*active) {
		*num_leaf += 1;
		return get_depth_of_common_ancestor((*num_leaf)-1, *num_leaf);
	} else {
		*num_leaf = 0;
		*active = 1;
		j = 1;
	}
	return j;
}

static inline int precompute_enc_ctx(const uint8_t* salt, uint32_t e, enc_ctx_ecb *ctx_enc_cache) {
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	uint32_t i;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	for (i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++) {
		TweakSalt(salt, tweaked_salt, 2, e, (MQOM2_PARAM_NB_EVALS_LOG-1)-i-1);
		ret = enc_key_sched_ecb(&ctx_enc_cache[i], tweaked_salt);
		ERR(ret, err);
	}
	ret = 0;
err:
	return ret;
#else
	(void) salt;
	(void) e;
	(void) ctx_enc_cache;
	return 0;
#endif
}

static inline int precompute_enc_ctx_pub(const uint8_t* salt, uint32_t e, enc_ctx_pub_ecb *ctx_enc_cache) {
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	uint32_t i;
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];
	for (i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++) {
		TweakSalt(salt, tweaked_salt, 2, e, (MQOM2_PARAM_NB_EVALS_LOG-1)-i-1);
		ret = enc_key_sched_pub_ecb(&ctx_enc_cache[i], tweaked_salt);
		ERR(ret, err);
	}
	ret = 0;
err:
	return ret;
#else
	(void) salt;
	(void) e;
	(void) ctx_enc_cache;
	return 0;
#endif
}

static inline enc_ctx_ecb* recover_enc_ctx(const uint8_t* salt, uint32_t e, uint32_t j, enc_ctx_ecb* ctx_enc_local, enc_ctx_ecb *ctx_enc_cache) {
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];

	// Test if "(MQOM2_PARAM_NB_EVALS_LOG-1) - j < GGMTREE_NB_ENC_CTX_IN_MEMORY"
	if((MQOM2_PARAM_NB_EVALS_LOG-1) < j + GGMTREE_NB_ENC_CTX_IN_MEMORY) {
		return &ctx_enc_cache[(MQOM2_PARAM_NB_EVALS_LOG-1) - j];
	}

	TweakSalt(salt, tweaked_salt, 2, e, j-1);
	ret = enc_key_sched_ecb(ctx_enc_local, tweaked_salt);
	return (ret == 0) ? ctx_enc_local : NULL;
}

static inline enc_ctx_pub_ecb* recover_enc_ctx_pub(const uint8_t* salt, uint32_t e, uint32_t j, enc_ctx_pub_ecb* ctx_enc_local, enc_ctx_pub_ecb *ctx_enc_cache) {
	int ret = -1;
	uint8_t tweaked_salt[MQOM2_PARAM_SEED_SIZE];

	// Test if "(MQOM2_PARAM_NB_EVALS_LOG-1) - j < GGMTREE_NB_ENC_CTX_IN_MEMORY"
	if ((MQOM2_PARAM_NB_EVALS_LOG-1) < j + GGMTREE_NB_ENC_CTX_IN_MEMORY) {
		return &ctx_enc_cache[(MQOM2_PARAM_NB_EVALS_LOG-1) - j];
	}

	TweakSalt(salt, tweaked_salt, 2, e, j-1);
	ret = enc_key_sched_pub_ecb(ctx_enc_local, tweaked_salt);
	return (ret == 0) ? ctx_enc_local : NULL;
}

int GGMTree_InitIncrementalExpansion_batch(ggmtree_ctx_batch_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e) {
	int ret = -1;

	if((ctx == NULL) || (delta == NULL) || (rseed == NULL)){
		goto err;
	}

	ret = precompute_enc_ctx(salt, e, ctx->ctx_enc);
	ERR(ret, err);

	memcpy(ctx->path[0], delta, MQOM2_PARAM_SEED_SIZE);
	memcpy(ctx->path[1], rseed, MQOM2_PARAM_SEED_SIZE);

	ctx->active = 0;
	ctx->num_leaf = 0;
	ctx->e = e;
	ctx->salt = salt;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeafs_batch(ggmtree_ctx_batch_t* ctx, uint8_t lseeds[][MQOM2_PARAM_SEED_SIZE]) {
	enc_ctx_ecb DECL_VAR(ctx_enc);
	enc_ctx_ecb* ctx_enc_ptr = NULL;
	int ret = -1;

	/* Identify the already-computed ancestor */
	uint32_t j;
	j = increment_num_leaf(&ctx->active, &ctx->num_leaf);
	if(ctx->num_leaf > 0) {
		xor_blocks(ctx->path[j - 1], ctx->path[j], ctx->path[j]);
	}

	/* Expand the path to the subtree root */
	for (; j < MQOM2_PARAM_NB_EVALS_LOG - GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG; j++) {
		// Get the right cipher context
		ctx_enc_ptr = recover_enc_ctx(ctx->salt, ctx->e, j, &ctx_enc, ctx->ctx_enc);
		ERR_NULL(ctx_enc_ptr, err);

		// Derive the next node on the root-to-leaf path
		SeedDerive_ecb(ctx_enc_ptr, ctx->path[j], ctx->path[j + 1]);
	}
	memcpy(lseeds[0], ctx->path[MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG], MQOM2_PARAM_SEED_SIZE);

	/* Expand the entire subtree */
	uint8_t parent_nodes[GGMTREE_NB_SIMULTANEOUS_LEAVES/2][MQOM2_PARAM_SEED_SIZE];
	uint32_t nb_nodes = 1;
	for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		// Get the right cipher context
		ctx_enc_ptr = recover_enc_ctx(ctx->salt, ctx->e, j, &ctx_enc, ctx->ctx_enc);
		ERR_NULL(ctx_enc_ptr, err);

		// Derive the next node level
		memcpy(parent_nodes, lseeds, nb_nodes*MQOM2_PARAM_SEED_SIZE);
		DeriveSeeds_ecb(ctx_enc_ptr, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])parent_nodes, lseeds, nb_nodes);
		nb_nodes *= 2;
	}

	ret = 0;
err:
	enc_clean_ctx_ecb(&ctx_enc);
	return ret;
}

int GGMTree_InitIncrementalPartialExpansion_batch(ggmtree_ctx_partial_batch_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star) {
	int ret = -1;

	if((ctx == NULL) || (path == NULL)){
		goto err;
	}

	ret = precompute_enc_ctx_pub(salt, e, ctx->ctx_enc);
	ERR(ret, err);

	memcpy((uint8_t*) ctx->opening, (uint8_t*) path, sizeof(ctx->opening));

	ctx->active = 0;
	ctx->num_leaf = 0;
	ctx->e = e;
	ctx->salt = salt;
	ctx->i_star = i_star;

	ret = 0;
err:
	return ret;
}

int GGMTree_GetNextLeafsPartial_batch(ggmtree_ctx_partial_batch_t* ctx, uint8_t lseeds[][MQOM2_PARAM_SEED_SIZE]) {
	enc_ctx_pub_ecb DECL_VAR(ctx_enc);
	enc_ctx_pub_ecb* ctx_enc_ptr = NULL;
	int ret = -1;

	/* Identify the already-computed ancestor */
	uint32_t j;
	j = increment_num_leaf(&ctx->active, &ctx->num_leaf);

	uint32_t is_hidden_node = (ctx->num_leaf == (ctx->i_star>>GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG));
	if(is_hidden_node) {
		memset(lseeds[0], 0, MQOM2_PARAM_SEED_SIZE);
		j = MQOM2_PARAM_NB_EVALS_LOG - GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG;
	} else {
		uint32_t higher = get_depth_of_common_ancestor(ctx->num_leaf, (ctx->i_star>>GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG));
		if (j <= higher) {
			memcpy(ctx->path[higher], ctx->opening[MQOM2_PARAM_NB_EVALS_LOG - higher], MQOM2_PARAM_SEED_SIZE);
			j = higher;
		} else { // j > higher
			xor_blocks(ctx->path[j - 1], ctx->path[j], ctx->path[j]);
		}
		for (; j < MQOM2_PARAM_NB_EVALS_LOG - GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG; j++) {
			// Get the right cipher context
			ctx_enc_ptr = recover_enc_ctx_pub(ctx->salt, ctx->e, j, &ctx_enc, ctx->ctx_enc);
			ERR_NULL(ctx_enc_ptr, err);

			// Derive the next node on the root-to-leaf path
			SeedDerive_pub_ecb(ctx_enc_ptr, ctx->path[j], ctx->path[j + 1]);
		}
		memcpy(lseeds[0], ctx->path[MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_SIMULTANEOUS_LEAVES_LOG], MQOM2_PARAM_SEED_SIZE);
	}

	/* Expand the entire (partial) subtree */
	uint8_t parent_nodes[GGMTREE_NB_SIMULTANEOUS_LEAVES/2][MQOM2_PARAM_SEED_SIZE];
	uint32_t nb_nodes = 1;
	for (; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		// Get the right cipher context
		ctx_enc_ptr = recover_enc_ctx_pub(ctx->salt, ctx->e, j, &ctx_enc, ctx->ctx_enc);
		ERR_NULL(ctx_enc_ptr, err);

		// Derive the next node level
		memcpy(parent_nodes, lseeds, nb_nodes*MQOM2_PARAM_SEED_SIZE);
		DeriveSeeds_pub_ecb(ctx_enc_ptr, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])parent_nodes, lseeds, nb_nodes);
		nb_nodes *= 2;

		// Correct with opening
		if(is_hidden_node) {
			uint32_t hidden_node_idx = ((ctx->i_star>>((MQOM2_PARAM_NB_EVALS_LOG-1)-j)) & (nb_nodes-1));
			memcpy(lseeds[hidden_node_idx ^ 0x01], ctx->opening[(MQOM2_PARAM_NB_EVALS_LOG-1)-j], MQOM2_PARAM_SEED_SIZE);
			memset(lseeds[hidden_node_idx], 0, MQOM2_PARAM_SEED_SIZE);
		}
	}

	ret = 0;
err:
	enc_clean_ctx_pub_ecb(&ctx_enc);
	return ret;
}
