#include "blc_memopt.h"
#include "benchmark.h"
#include "blc_memopt_common.h"
#include "blc_memopt_x1.h"
#include "blc_memopt_x1_folding.h"
#include "blc_memopt_x1_seedcommit.h"

#ifdef GGM_TREE_NO_BATCHING
#include "ggm_tree_incr.h"
#else
#include "ggm_tree_incr_batch.h"
#if GGMTREE_NB_SIMULTANEOUS_LEAVES % BLC_NB_LEAF_SEEDS_IN_PARALLEL != 0
#error BLC_NB_LEAF_SEEDS_IN_PARALLEL should divide GGMTREE_NB_SIMULTANEOUS_LEAVES.
#endif
#endif

int BLC_Commit_x1_memopt(uint32_t e, const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const field_base_elt x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint8_t com[MQOM2_PARAM_DIGEST_SIZE], uint8_t partial_delta_x[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N)-MQOM2_PARAM_SEED_SIZE], field_ext_elt x0[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], field_ext_elt u0[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], field_ext_elt u1[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)]) {
	int ret = -1;
	uint32_t i, i_;

#ifdef GGM_TREE_NO_BATCHING
	uint8_t lseeds[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	ggmtree_ctx_t DECL_VAR(ggm_tree);
#else
	uint8_t lseeds[GGMTREE_NB_SIMULTANEOUS_LEAVES][MQOM2_PARAM_SEED_SIZE];
	ggmtree_ctx_batch_t DECL_VAR(ggm_tree);
#endif
	folding_sign_t folding;
	seedcommit_sign_ctx_t DECL_VAR(seedcommit_ctx);

	// Initialize the GGM tree
	__BENCHMARK_START__(BS_BLC_EXPAND_TREE);
#ifdef GGM_TREE_NO_BATCHING
	ret = GGMTree_InitIncrementalExpansion(&ggm_tree, salt, rseed, delta, e);
#else
	ret = GGMTree_InitIncrementalExpansion_batch(&ggm_tree, salt, rseed, delta, e);
#endif
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_EXPAND_TREE);

	// Initialize the hash context
	ret = init_seedcommit_sign(&seedcommit_ctx, salt, e);
	ERR(ret, err);

	ret = InitializeFolding_sign(&folding, salt, e);
	ERR(ret, err);
#ifdef GGM_TREE_NO_BATCHING
	for (i = 0; i < MQOM2_PARAM_NB_EVALS; i+= BLC_NB_LEAF_SEEDS_IN_PARALLEL) {

		// Derive the next leaf seeds
		__BENCHMARK_START__(BS_BLC_EXPAND_TREE);
		for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
			ret = GGMTree_GetNextLeaf(&ggm_tree, lseeds[i_]);
			ERR(ret, err);
		}
		__BENCHMARK_STOP__(BS_BLC_EXPAND_TREE);

		// Compute the individual commitments for all the seed leafs,
		// and incrementally hash them.
		ret = SeedCommitThenAbsorb_sign(&seedcommit_ctx, lseeds);
		ERR(ret, err);

		// Expand each seed and accumulate the expanded tapes
		ret = SeedExpandThenAccumulate_sign(&folding, i, lseeds);
		ERR(ret, err);
	}
#else
	for (i = 0; i < MQOM2_PARAM_NB_EVALS; i+= GGMTREE_NB_SIMULTANEOUS_LEAVES) {
		// Derive the next leaf seeds
		__BENCHMARK_START__(BS_BLC_EXPAND_TREE);
		ret = GGMTree_GetNextLeafs_batch(&ggm_tree, lseeds);
		ERR(ret, err);
		__BENCHMARK_STOP__(BS_BLC_EXPAND_TREE);

		for (i_ = 0; i_ < (GGMTREE_NB_SIMULTANEOUS_LEAVES/BLC_NB_LEAF_SEEDS_IN_PARALLEL); i_++) {
			// Compute the individual commitments for all the seed leafs,
			// and incrementally hash them.
			ret = SeedCommitThenAbsorb_sign(&seedcommit_ctx, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])&lseeds[i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL]);
			ERR(ret, err);

			// Expand each seed and accumulate the expanded tapes
			ret = SeedExpandThenAccumulate_sign(&folding, i + i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])&lseeds[i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL]);
			ERR(ret, err);
		}
	}
#endif

	// Finalize the folding to get the committed polynomials
	__BENCHMARK_START__(BS_BLC_ARITH);
	FinalizeFolding_sign(&folding, x, partial_delta_x, x0, u0, u1);
	__BENCHMARK_STOP__(BS_BLC_ARITH);

	// Get the global commitment digest
	__BENCHMARK_START__(BS_BLC_XOF);
	ret = xof_squeeze(&seedcommit_ctx.xof_ctx, com, MQOM2_PARAM_DIGEST_SIZE);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_XOF);

	ret = 0;
err:
	seedcommit_sign_clean_ctx(&seedcommit_ctx);
#ifdef GGM_TREE_NO_BATCHING
	ggmtree_ctx_t_clean(&ggm_tree);
#else
	ggmtree_ctx_batch_t_clean(&ggm_tree);
#endif
	folding_sign_clean_ctx(&folding);
	return ret;
}

int BLC_Eval_x1_memopt(uint32_t e, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_SEED_SIZE*MQOM2_PARAM_NB_EVALS_LOG], const uint8_t out_ls_com[MQOM2_PARAM_DIGEST_SIZE], const uint8_t partial_delta_x[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N)-MQOM2_PARAM_SEED_SIZE], uint16_t i_star, uint8_t com[MQOM2_PARAM_DIGEST_SIZE], field_ext_elt x_eval[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], field_ext_elt u_eval[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)]) {
	int ret = -1;
	uint32_t i, i_;

#ifdef GGM_TREE_NO_BATCHING
	uint8_t lseeds[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	ggmtree_ctx_partial_t DECL_VAR(ggm_tree);
#else
	uint8_t lseeds[GGMTREE_NB_SIMULTANEOUS_LEAVES][MQOM2_PARAM_SEED_SIZE];
	ggmtree_ctx_partial_batch_t DECL_VAR(ggm_tree);
#endif
	folding_verify_t folding;
	seedcommit_verify_ctx_t DECL_VAR(seedcommit_ctx);


	// Initialize the GGM tree
#ifdef GGM_TREE_NO_BATCHING
	ret = GGMTree_InitIncrementalPartialExpansion(&ggm_tree, salt, (const uint8_t(*)[MQOM2_PARAM_SEED_SIZE]) path, e, i_star);
#else
	ret = GGMTree_InitIncrementalPartialExpansion_batch(&ggm_tree, salt, (const uint8_t(*)[MQOM2_PARAM_SEED_SIZE]) path, e, i_star);
#endif
	ERR(ret, err);

	// Initialize the hash context
	ret = init_seedcommit_verify(&seedcommit_ctx, salt, e, i_star, out_ls_com);
	ERR(ret, err);

	ret = InitializeFolding_verify(&folding, salt, e);
	ERR(ret, err);
#ifdef GGM_TREE_NO_BATCHING
	for (i = 0; i < MQOM2_PARAM_NB_EVALS; i+= BLC_NB_LEAF_SEEDS_IN_PARALLEL) {
		for (i_ = 0; i_<BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
			GGMTree_GetNextLeafPartial(&ggm_tree, lseeds[i_]);
		}
		
		ret = SeedCommitThenAbsorb_verify(&seedcommit_ctx, i, lseeds);
		ERR(ret, err);
		ret = SeedExpandThenAccumulate_verify(&folding, i, lseeds, i_star);
		ERR(ret, err);
	}
#else
	for (i = 0; i < MQOM2_PARAM_NB_EVALS; i+= GGMTREE_NB_SIMULTANEOUS_LEAVES) {
		ret = GGMTree_GetNextLeafsPartial_batch(&ggm_tree, lseeds);
		ERR(ret, err);

		for (i_ = 0; i_ < (GGMTREE_NB_SIMULTANEOUS_LEAVES/BLC_NB_LEAF_SEEDS_IN_PARALLEL); i_++) {
			ret = SeedCommitThenAbsorb_verify(&seedcommit_ctx, i + i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])&lseeds[i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL]);
			ERR(ret, err);
			ret = SeedExpandThenAccumulate_verify(&folding, i + i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL, (const uint8_t (*)[MQOM2_PARAM_SEED_SIZE])&lseeds[i_*BLC_NB_LEAF_SEEDS_IN_PARALLEL], i_star);
			ERR(ret, err);
		}
	}
#endif
	FinalizeFolding_verify(&folding, i_star, partial_delta_x, x_eval, u_eval);

	ret = xof_squeeze(&seedcommit_ctx.xof_ctx, com, MQOM2_PARAM_DIGEST_SIZE);
	ERR(ret, err);

	ret = 0;
err:
	seedcommit_verify_clean_ctx(&seedcommit_ctx);
#ifdef GGM_TREE_NO_BATCHING
	ggmtree_ctx_partial_t_clean(&ggm_tree);
#else
	ggmtree_ctx_partial_batch_t_clean(&ggm_tree);
#endif
	folding_verify_clean_ctx(&folding);
	return ret;
}
