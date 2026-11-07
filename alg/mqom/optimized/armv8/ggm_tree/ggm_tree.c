#include "ggm_tree.h"
#include "ggm_tree_common.h"

/* NOTE: in the "node" tree representation, we accept to remove the two first cells to simplify the indices computations.
 * The first cell is used to avoid the "0" index, and the second cell is the root of the correlated tree that is not used.
 * */
int GGMTree_Expand(const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e, uint8_t node[MQOM2_PARAM_FULL_TREE_SIZE + 1][MQOM2_PARAM_SEED_SIZE], uint8_t lseed[MQOM2_PARAM_NB_EVALS][MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	/* j is the level in the tree, k is the index in the nodes array */
	uint32_t j, k;
	enc_ctx DECL_VAR(ctx);
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];

	/* Some sanity check */
	if ((1 << MQOM2_PARAM_NB_EVALS_LOG) != MQOM2_PARAM_NB_EVALS) {
		ret = -1;
		goto err;
	}

	/* The root node is not defined, we begin with the first level */
	memcpy(node[2], rseed, MQOM2_PARAM_SEED_SIZE);
	xor_blocks(node[2], delta, node[3]);

	/* Now deal with the other levels */
	for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		/* Level 1 has 2 derivations, levels > 1 allow for more derivations */
		uint32_t num_derivations;
		if (j == 1) {
			num_derivations = 2;
		}
#ifdef USE_ENC_X8
		else if (j == 2) {
			num_derivations = 4;
		} else {
			num_derivations = 8;
		}
#else
		else {
			num_derivations = 4;
		}
#endif
		/* For a whole level in the tree, we use the same key, which is the tweaked salt, hence
		 * the common key schedule */
		TweakSalt(salt, tweaked_salt, 2, e, j - 1);
		ret = enc_key_sched(&ctx, tweaked_salt);
		ERR(ret, err);
		/* NOTE: hereafter, when using x2 or x4 we are sure that there
		         * are no leftover nodes. The encryptions performed correspong to the
		         * SeedDerive procedure.
		 */
		for (k = ((uint32_t)1 << j); k < ((uint32_t)1 << (j + 1)); k += num_derivations) {
			switch (num_derivations) {
			case 2: {
				uint8_t *pt1 = node[k];
				uint8_t *pt2 = node[k + 1];
				uint8_t *ct1 = node[2 * k];
				uint8_t *ct2 = node[2 * (k + 1)];
				/* We perform an x2 on 2 nodes */
				SeedDerive_x2(&ctx, &ctx,
				              /* Input nodes to encrypt */
				              pt1, pt2,
				              /* Output nodes to encrypt */
				              ct1, ct2);
				/* Compute the corresponding xors */
				uint8_t *rnode_ct1 = node[(2 * k) + 1];
				uint8_t *rnode_ct2 = node[(2 * (k + 1)) + 1];
				xor_blocks(ct1, pt1, rnode_ct1);
				xor_blocks(ct2, pt2, rnode_ct2);
				break;
			}
			case 4: {
				/* We perform an x4 on 4 nodes */
				uint8_t *pt1 = node[k];
				uint8_t *pt2 = node[k + 1];
				uint8_t *pt3 = node[k + 2];
				uint8_t *pt4 = node[k + 3];
				uint8_t *ct1 = node[2 * k];
				uint8_t *ct2 = node[2 * (k + 1)];
				uint8_t *ct3 = node[2 * (k + 2)];
				uint8_t *ct4 = node[2 * (k + 3)];
				/* We perform an x4 on 4 nodes */
				SeedDerive_x4(&ctx, &ctx, &ctx, &ctx,
				              /* Input nodes to encrypt */
				              pt1, pt2, pt3, pt4,
				              /* Output nodes to encrypt */
				              ct1, ct2, ct3, ct4);
				/* Compute the corresponding xors */
				uint8_t *rnode_ct1 = node[(2 * k) + 1];
				uint8_t *rnode_ct2 = node[(2 * (k + 1)) + 1];
				uint8_t *rnode_ct3 = node[(2 * (k + 2)) + 1];
				uint8_t *rnode_ct4 = node[(2 * (k + 3)) + 1];
				xor_blocks(ct1, pt1, rnode_ct1);
				xor_blocks(ct2, pt2, rnode_ct2);
				xor_blocks(ct3, pt3, rnode_ct3);
				xor_blocks(ct4, pt4, rnode_ct4);
				break;
			}
#ifdef USE_ENC_X8
			case 8: {
				/* We perform an x8 on 8 nodes */
				uint8_t *pt1 = node[k];
				uint8_t *pt2 = node[k + 1];
				uint8_t *pt3 = node[k + 2];
				uint8_t *pt4 = node[k + 3];
				uint8_t *pt5 = node[k + 4];
				uint8_t *pt6 = node[k + 5];
				uint8_t *pt7 = node[k + 6];
				uint8_t *pt8 = node[k + 7];
				uint8_t *ct1 = node[2 * k];
				uint8_t *ct2 = node[2 * (k + 1)];
				uint8_t *ct3 = node[2 * (k + 2)];
				uint8_t *ct4 = node[2 * (k + 3)];
				uint8_t *ct5 = node[2 * (k + 4)];
				uint8_t *ct6 = node[2 * (k + 5)];
				uint8_t *ct7 = node[2 * (k + 6)];
				uint8_t *ct8 = node[2 * (k + 7)];
				/* We perform an x8 on 8 nodes */
				SeedDerive_x8(&ctx, &ctx, &ctx, &ctx, &ctx, &ctx, &ctx, &ctx,
				              /* Input nodes to encrypt */
				              pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8,
				              /* Output nodes to encrypt */
				              ct1, ct2, ct3, ct4, ct5, ct6, ct7, ct8);
				/* Compute the corresponding xors */
				uint8_t *rnode_ct1 = node[(2 * k) + 1];
				uint8_t *rnode_ct2 = node[(2 * (k + 1)) + 1];
				uint8_t *rnode_ct3 = node[(2 * (k + 2)) + 1];
				uint8_t *rnode_ct4 = node[(2 * (k + 3)) + 1];
				uint8_t *rnode_ct5 = node[(2 * (k + 4)) + 1];
				uint8_t *rnode_ct6 = node[(2 * (k + 5)) + 1];
				uint8_t *rnode_ct7 = node[(2 * (k + 6)) + 1];
				uint8_t *rnode_ct8 = node[(2 * (k + 7)) + 1];
				xor_blocks(ct1, pt1, rnode_ct1);
				xor_blocks(ct2, pt2, rnode_ct2);
				xor_blocks(ct3, pt3, rnode_ct3);
				xor_blocks(ct4, pt4, rnode_ct4);
				xor_blocks(ct5, pt5, rnode_ct5);
				xor_blocks(ct6, pt6, rnode_ct6);
				xor_blocks(ct7, pt7, rnode_ct7);
				xor_blocks(ct8, pt8, rnode_ct8);
				break;
			}
#endif
			default: {
				ret = -1;
				goto err;
			}
			}
		}
	}

	/* The lseed are the leaves of the tree, copied from it */
	memcpy(lseed, node[MQOM2_PARAM_FULL_TREE_SIZE + 1 - MQOM2_PARAM_NB_EVALS], MQOM2_PARAM_NB_EVALS * MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx(&ctx);
	return ret;
}

int GGMTree_Open(const uint8_t node[MQOM2_PARAM_FULL_TREE_SIZE + 1][MQOM2_PARAM_SEED_SIZE], uint32_t i_star, uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	uint32_t i, j;

	/* Sanity check */
	if (i_star >= MQOM2_PARAM_NB_EVALS) {
		ret = -1;
		goto err;
	}

	i = MQOM2_PARAM_NB_EVALS + i_star;
	for (j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		memcpy(path[j], node[i ^ 1], MQOM2_PARAM_SEED_SIZE);
		i = (i / 2);
	}

	ret = 0;
err:
	return ret;
}

/* XXX: TODO: for now, we allocate the full tree to compute the leaves seeds. It is possible to perform this *in place* with
 * a dedicated index handling, allowing to save memory space as well as a memcpy */
/* XXX: NOTE: we can use public encryption API here as this function is used for verification */
int GGMTree_PartiallyExpand(const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star, uint8_t lseed[MQOM2_PARAM_NB_EVALS][MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	/* j is the level in the tree, k is the index in the nodes array */
	uint32_t i, j, k;
	enc_ctx_pub DECL_VAR(ctx);
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	/* Locally allocate the full tree */
	uint8_t node[MQOM2_PARAM_FULL_TREE_SIZE + 1][MQOM2_PARAM_SEED_SIZE];
	/* We use a shadow map to indicate bot values */
	uint8_t node_map[MQOM2_PARAM_FULL_TREE_SIZE + 1];
	memset(node_map, 0, MQOM2_PARAM_FULL_TREE_SIZE + 1);

	/* Some sanity checks */
	if ((1 << MQOM2_PARAM_NB_EVALS_LOG) != MQOM2_PARAM_NB_EVALS) {
		ret = -1;
		goto err;
	}
	if (i_star >= MQOM2_PARAM_NB_EVALS) {
		ret = -1;
		goto err;
	}

	/* Copy the path */
	i = MQOM2_PARAM_NB_EVALS + i_star;
	for (j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		memcpy(node[i ^ 1], path[j], MQOM2_PARAM_SEED_SIZE);
		node_map[i ^ 1] = 1;
		i = (i / 2);
	}

	/* Compute the other nodes when possible.
	 * */
	for (j = 1; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		/* Level 1 has 2 derivations, levels > 1 allow for more derivations */
		uint32_t num_derivations;
		if (j == 1) {
			num_derivations = 2;
		}
#ifdef USE_ENC_X8
		else if (j == 2) {
			num_derivations = 4;
		} else {
			num_derivations = 8;
		}
#else
		else {
			num_derivations = 4;
		}
#endif
		/* For a whole level in the tree, we use the same key, which is the tweaked salt, hence
		 * the common key schedule */
		TweakSalt(salt, tweaked_salt, 2, e, j - 1);
		ret = enc_key_sched_pub(&ctx, tweaked_salt);
		ERR(ret, err);
		/* NOTE: hereafter, when using x2 or x4 we are sure that there
		         * are no leftover nodes. The encryptions performed correspong to the
		         * SeedDerive procedure.
		 */
		for (k = ((uint32_t)1 << j); k < ((uint32_t)1 << (j + 1)); k += num_derivations) {
			switch (num_derivations) {
			case 2: {
				uint8_t *pt1 = node_map[k]     ? node[k]     : NULL;
				uint8_t *pt2 = node_map[k + 1] ? node[k + 1] : NULL;
				uint8_t *ct1 = node[2 * k];
				uint8_t *ct2 = node[2 * (k + 1)];
				/* We perform an x2 on 2 nodes */
				SeedDerive_x2_pub(&ctx, &ctx,
				                  /* Input nodes to encrypt */
				                  pt1, pt2,
				                  /* Output nodes to encrypt */
				                  ct1, ct2);
				/* Compute the corresponding xors */
				uint8_t *rnode_ct1 = node[(2 * k) + 1];
				uint8_t *rnode_ct2 = node[(2 * (k + 1)) + 1];
				if (pt1) {
					xor_blocks(ct1, pt1, rnode_ct1);
					node_map[(2 * k)] = 1;
					node_map[(2 * k) + 1] = 1;
				}
				if (pt2) {
					xor_blocks(ct2, pt2, rnode_ct2);
					node_map[(2 * (k + 1))] = 1;
					node_map[(2 * (k + 1)) + 1] = 1;
				}
				break;
			}
			case 4: {
				/* We perform an x4 on 4 nodes */
				uint8_t *pt1 = node_map[k]     ? node[k]     : NULL;
				uint8_t *pt2 = node_map[k + 1] ? node[k + 1] : NULL;
				uint8_t *pt3 = node_map[k + 2] ? node[k + 2] : NULL;
				uint8_t *pt4 = node_map[k + 3] ? node[k + 3] : NULL;
				uint8_t *ct1 = node[2 * k];
				uint8_t *ct2 = node[2 * (k + 1)];
				uint8_t *ct3 = node[2 * (k + 2)];
				uint8_t *ct4 = node[2 * (k + 3)];
				/* We perform an x4 on 4 nodes */
				SeedDerive_x4_pub(&ctx, &ctx, &ctx, &ctx,
				                  /* Input nodes to encrypt */
				                  pt1, pt2, pt3, pt4,
				                  /* Output nodes to encrypt */
				                  ct1, ct2, ct3, ct4);
				/* Compute the corresponding xors */
				uint8_t *rnode_ct1 = node[(2 * k) + 1];
				uint8_t *rnode_ct2 = node[(2 * (k + 1)) + 1];
				uint8_t *rnode_ct3 = node[(2 * (k + 2)) + 1];
				uint8_t *rnode_ct4 = node[(2 * (k + 3)) + 1];
				if (pt1) {
					xor_blocks(ct1, pt1, rnode_ct1);
					node_map[(2 * k)] = 1;
					node_map[(2 * k) + 1] = 1;
				}
				if (pt2) {
					xor_blocks(ct2, pt2, rnode_ct2);
					node_map[(2 * (k + 1))] = 1;
					node_map[(2 * (k + 1)) + 1] = 1;
				}
				if (pt3) {
					xor_blocks(ct3, pt3, rnode_ct3);
					node_map[(2 * (k + 2))] = 1;
					node_map[(2 * (k + 2)) + 1] = 1;
				}
				if (pt4) {
					xor_blocks(ct4, pt4, rnode_ct4);
					node_map[(2 * (k + 3))] = 1;
					node_map[(2 * (k + 3)) + 1] = 1;
				}
				break;
			}
#ifdef USE_ENC_X8
			case 8: {
				/* We perform an x8 on 8 nodes */
				uint8_t *pt1 = node_map[k]     ? node[k]     : NULL;
				uint8_t *pt2 = node_map[k + 1] ? node[k + 1] : NULL;
				uint8_t *pt3 = node_map[k + 2] ? node[k + 2] : NULL;
				uint8_t *pt4 = node_map[k + 3] ? node[k + 3] : NULL;
				uint8_t *pt5 = node_map[k + 4] ? node[k + 4] : NULL;
				uint8_t *pt6 = node_map[k + 5] ? node[k + 5] : NULL;
				uint8_t *pt7 = node_map[k + 6] ? node[k + 6] : NULL;
				uint8_t *pt8 = node_map[k + 7] ? node[k + 7] : NULL;
				uint8_t *ct1 = node[2 * k];
				uint8_t *ct2 = node[2 * (k + 1)];
				uint8_t *ct3 = node[2 * (k + 2)];
				uint8_t *ct4 = node[2 * (k + 3)];
				uint8_t *ct5 = node[2 * (k + 4)];
				uint8_t *ct6 = node[2 * (k + 5)];
				uint8_t *ct7 = node[2 * (k + 6)];
				uint8_t *ct8 = node[2 * (k + 7)];
				/* We perform an x4 on 4 nodes */
				SeedDerive_x8_pub(&ctx, &ctx, &ctx, &ctx, &ctx, &ctx, &ctx, &ctx,
				                  /* Input nodes to encrypt */
				                  pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8,
				                  /* Output nodes to encrypt */
				                  ct1, ct2, ct3, ct4, ct5, ct6, ct7, ct8);
				/* Compute the corresponding xors */
				uint8_t *rnode_ct1 = node[(2 * k) + 1];
				uint8_t *rnode_ct2 = node[(2 * (k + 1)) + 1];
				uint8_t *rnode_ct3 = node[(2 * (k + 2)) + 1];
				uint8_t *rnode_ct4 = node[(2 * (k + 3)) + 1];
				uint8_t *rnode_ct5 = node[(2 * (k + 4)) + 1];
				uint8_t *rnode_ct6 = node[(2 * (k + 5)) + 1];
				uint8_t *rnode_ct7 = node[(2 * (k + 6)) + 1];
				uint8_t *rnode_ct8 = node[(2 * (k + 7)) + 1];
				if (pt1) {
					xor_blocks(ct1, pt1, rnode_ct1);
					node_map[(2 * k)] = 1;
					node_map[(2 * k) + 1] = 1;
				}
				if (pt2) {
					xor_blocks(ct2, pt2, rnode_ct2);
					node_map[(2 * (k + 1))] = 1;
					node_map[(2 * (k + 1)) + 1] = 1;
				}
				if (pt3) {
					xor_blocks(ct3, pt3, rnode_ct3);
					node_map[(2 * (k + 2))] = 1;
					node_map[(2 * (k + 2)) + 1] = 1;
				}
				if (pt4) {
					xor_blocks(ct4, pt4, rnode_ct4);
					node_map[(2 * (k + 3))] = 1;
					node_map[(2 * (k + 3)) + 1] = 1;
				}
				if (pt5) {
					xor_blocks(ct5, pt5, rnode_ct5);
					node_map[(2 * (k + 4))] = 1;
					node_map[(2 * (k + 4)) + 1] = 1;
				}
				if (pt6) {
					xor_blocks(ct6, pt6, rnode_ct6);
					node_map[(2 * (k + 5))] = 1;
					node_map[(2 * (k + 5)) + 1] = 1;
				}
				if (pt7) {
					xor_blocks(ct7, pt7, rnode_ct7);
					node_map[(2 * (k + 6))] = 1;
					node_map[(2 * (k + 6)) + 1] = 1;
				}
				if (pt8) {
					xor_blocks(ct8, pt8, rnode_ct8);
					node_map[(2 * (k + 7))] = 1;
					node_map[(2 * (k + 7)) + 1] = 1;
				}
				break;
			}
#endif
			default: {
				ret = -1;
				goto err;
			}
			}
		}

	}

	/* The lseed are the leaves of the tree, copied from it */
	memcpy(lseed, node[MQOM2_PARAM_FULL_TREE_SIZE + 1 - MQOM2_PARAM_NB_EVALS], MQOM2_PARAM_NB_EVALS * MQOM2_PARAM_SEED_SIZE);
	/* Set the hidden leave to zero */
	memset(lseed[i_star], 0, MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx_pub(&ctx);
	return ret;
}

int GGMTree_ExpandPath(const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star, uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint8_t lseed[MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	uint32_t j;
	enc_ctx DECL_VAR(ctx);
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];

	/* Sanity check */
	if (i_star >= MQOM2_PARAM_NB_EVALS) {
		ret = -1;
		goto err;
	}

	uint32_t num_leaf = MQOM2_PARAM_NB_EVALS + i_star;

	uint8_t node[2][MQOM2_PARAM_SEED_SIZE];
	uint8_t parent[MQOM2_PARAM_SEED_SIZE];
	memcpy(parent, delta, MQOM2_PARAM_SEED_SIZE);

	for (j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		if (j == 0) {
			memcpy(node[0], rseed, MQOM2_PARAM_SEED_SIZE);
		} else {
			TweakSalt(salt, tweaked_salt, 2, e, j - 1);
			ret = enc_key_sched(&ctx, tweaked_salt);
			ERR(ret, err);
			SeedDerive(&ctx, parent, node[0]);
		}
		xor_blocks(node[0], parent, node[1]);

		uint32_t bit = (num_leaf >> (MQOM2_PARAM_NB_EVALS_LOG - 1 - j)) & 0x01;
		memcpy(path[MQOM2_PARAM_NB_EVALS_LOG - 1 - j], node[bit ^ 1], MQOM2_PARAM_SEED_SIZE);
		memcpy(parent, node[bit], MQOM2_PARAM_SEED_SIZE);
	}

	memcpy(lseed, parent, MQOM2_PARAM_SEED_SIZE);

	ret = 0;
err:
	enc_clean_ctx(&ctx);
	return ret;
}

