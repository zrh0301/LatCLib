#ifndef __GGM_TREE_INCR_X2_H__
#define __GGM_TREE_INCR_X2_H__

/* MQOM2 parameters */
#include "mqom2_parameters.h"
/* Encryption primitive */
#include "enc.h"
/* Common helpers */
#include "common.h"

#include "ggm_tree_common.h"

/* Deal with namespacing */
#define GGMTree_InitIncrementalExpansion_x2 MQOM_NAMESPACE(GGMTree_InitIncrementalExpansion_x2)
#define GGMTree_GetNextLeaf_x2 MQOM_NAMESPACE(GGMTree_GetNextLeaf_x2)
#define GGMTree_InitIncrementalPartialExpansion_x2 MQOM_NAMESPACE(GGMTree_InitIncrementalPartialExpansion_x2)
#define GGMTree_GetNextLeafPartial_x2 MQOM_NAMESPACE(GGMTree_GetNextLeafPartial_x2)

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint8_t tweaked_salts[2][GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY][MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_x2 *ctx_enc_x2;
#else
	enc_ctx_x2 ctx_enc_x2[GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[2][MQOM2_PARAM_NB_EVALS_LOG + 1][MQOM2_PARAM_SEED_SIZE];
} ggmtree_ctx_x2_t;

/* Cleaning function */
static inline void ggmtree_ctx_x2_t_clean(ggmtree_ctx_x2_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_x2(&ctx->ctx_enc_x2[i]);
	}
#endif
}

int GGMTree_InitIncrementalExpansion_x2(ggmtree_ctx_x2_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[2][MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], const uint32_t e[2]);

int GGMTree_GetNextLeaf_x2(ggmtree_ctx_x2_t* ctx, uint8_t lseed[2][MQOM2_PARAM_SEED_SIZE]);

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint8_t tweaked_salts[2][GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY][MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_pub **ctx_enc;
#else
	enc_ctx_pub ctx_enc[2][GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[2][MQOM2_PARAM_NB_EVALS_LOG + 1][MQOM2_PARAM_SEED_SIZE];
	const uint8_t opening[2][MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE];
	uint32_t i_star[2];
} ggmtree_ctx_partial_x2_t;

/* Cleaning function */
static inline void ggmtree_ctx_partial_x2_t_clean(ggmtree_ctx_partial_x2_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_pub(&ctx->ctx_enc[0][i]);
		enc_clean_ctx_pub(&ctx->ctx_enc[1][i]);
	}
#endif
}

int GGMTree_InitIncrementalPartialExpansion_x2(ggmtree_ctx_partial_x2_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t (*path[2])[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], const uint32_t e[2], const uint32_t i_star[2]);

int GGMTree_GetNextLeafPartial_x2(ggmtree_ctx_partial_x2_t* ctx, uint8_t lseed[2][MQOM2_PARAM_SEED_SIZE]);


#endif /* __GGM_TREE_INCR_X2_H__ */
