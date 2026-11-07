#ifndef __GGM_TREE_INCR_X1_BASE_H__
#define __GGM_TREE_INCR_X1_BASE_H__

/* MQOM2 parameters */
#include "mqom2_parameters.h"
/* Encryption primitive */
#include "enc.h"
/* Common helpers */
#include "common.h"

#include "ggm_tree_common.h"

/* Deal with namespacing */
#define GGMTree_InitIncrementalExpansion_base MQOM_NAMESPACE(GGMTree_InitIncrementalExpansion_base)
#define GGMTree_GetNextLeaf_base MQOM_NAMESPACE(GGMTree_GetNextLeaf_base)
#define GGMTree_InitIncrementalPartialExpansion_base MQOM_NAMESPACE(GGMTree_InitIncrementalPartialExpansion_base)
#define GGMTree_GetNextLeafPartial_base MQOM_NAMESPACE(GGMTree_GetNextLeafPartial_base)

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint32_t e;
	uint8_t tweaked_salts[GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY][MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx *ctx_enc;
#else
	enc_ctx ctx_enc[GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[MQOM2_PARAM_NB_EVALS_LOG + 1][MQOM2_PARAM_SEED_SIZE];
} ggmtree_ctx_base_t;

/* Cleaning function */
static inline void ggmtree_ctx_base_t_clean(ggmtree_ctx_base_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx(&ctx->ctx_enc[i]);
	}
#endif
}
int GGMTree_InitIncrementalExpansion_base(ggmtree_ctx_base_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e);

int GGMTree_GetNextLeaf_base(ggmtree_ctx_base_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]);

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint8_t tweaked_salts[GGMTREE_NB_TWEAKED_SALTS_IN_MEMORY][MQOM2_PARAM_SEED_SIZE];
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_pub *ctx_enc;
#else
	enc_ctx_pub ctx_enc[GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[MQOM2_PARAM_NB_EVALS_LOG + 1][MQOM2_PARAM_SEED_SIZE];
	uint8_t opening[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE];
	uint32_t i_star;
} ggmtree_ctx_partial_base_t;

/* Cleaning function */
static inline void ggmtree_ctx_partial_base_t_clean(ggmtree_ctx_partial_base_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_pub(&ctx->ctx_enc[i]);
	}
#endif
}

int GGMTree_InitIncrementalPartialExpansion_base(ggmtree_ctx_partial_base_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star);

int GGMTree_GetNextLeafPartial_base(ggmtree_ctx_partial_base_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]);


#endif /* __GGM_TREE_INCR_X1_BASE_H__ */
