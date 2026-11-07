#ifndef __GGM_TREE_INCR_X4_H__
#define __GGM_TREE_INCR_X4_H__

/* MQOM2 parameters */
#include "mqom2_parameters.h"
/* Encryption primitive */
#include "enc.h"
/* Common helpers */
#include "common.h"

#include "ggm_tree_common.h"

/* Deal with namespacing */
#define GGMTree_InitIncrementalExpansion_x4 MQOM_NAMESPACE(GGMTree_InitIncrementalExpansion_x4)
#define GGMTree_GetNextLeaf_x4 MQOM_NAMESPACE(GGMTree_GetNextLeaf_x4)
#define GGMTree_InitIncrementalPartialExpansion_x4 MQOM_NAMESPACE(GGMTree_InitIncrementalPartialExpansion_x4)
#define GGMTree_GetNextLeafPartial_x4 MQOM_NAMESPACE(GGMTree_GetNextLeafPartial_x4)

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint8_t tweaked_salts[4][GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY][MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_x4 *ctx_enc_x4;
#else
	enc_ctx_x4 ctx_enc_x4[GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[4][MQOM2_PARAM_NB_EVALS_LOG + 1][MQOM2_PARAM_SEED_SIZE];
} ggmtree_ctx_x4_t;

/* Cleaning function */
static inline void ggmtree_ctx_x4_t_clean(ggmtree_ctx_x4_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_x4(&ctx->ctx_enc_x4[i]);
	}
#endif
}

int GGMTree_InitIncrementalExpansion_x4(ggmtree_ctx_x4_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[4][MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], const uint32_t e[4]);

int GGMTree_GetNextLeaf_x4(ggmtree_ctx_x4_t* ctx, uint8_t lseed[4][MQOM2_PARAM_SEED_SIZE]);

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint8_t tweaked_salts[4][GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY][MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_pub **ctx_enc;
#else
	enc_ctx_pub ctx_enc[4][GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[4][MQOM2_PARAM_NB_EVALS_LOG + 1][MQOM2_PARAM_SEED_SIZE];
	const uint8_t opening[4][MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE];
	uint32_t i_star[4];
} ggmtree_ctx_partial_x4_t;

/* Cleaning function */
static inline void ggmtree_ctx_partial_x4_t_clean(ggmtree_ctx_partial_x4_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_pub(&ctx->ctx_enc[0][i]);
		enc_clean_ctx_pub(&ctx->ctx_enc[1][i]);
		enc_clean_ctx_pub(&ctx->ctx_enc[2][i]);
		enc_clean_ctx_pub(&ctx->ctx_enc[3][i]);
	}
#endif
}

int GGMTree_InitIncrementalPartialExpansion_x4(ggmtree_ctx_partial_x4_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t (*path[4])[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], const uint32_t e[4], const uint32_t i_star[4]);

int GGMTree_GetNextLeafPartial_x4(ggmtree_ctx_partial_x4_t* ctx, uint8_t lseed[4][MQOM2_PARAM_SEED_SIZE]);

#endif /* __GGM_TREE_INCR_X4_H__ */
