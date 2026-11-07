#ifndef __GGM_TREE_INCR_X1_ADV_H__
#define __GGM_TREE_INCR_X1_ADV_H__

/* MQOM2 parameters */
#include "mqom2_parameters.h"
/* Encryption primitive */
#include "enc.h"
/* Common helpers */
#include "common.h"

#include "ggm_tree_common.h"

/* Deal with namespacing */
#define GGMTree_InitIncrementalExpansion_adv MQOM_NAMESPACE(GGMTree_InitIncrementalExpansion_adv)
#define GGMTree_GetNextLeaf_adv MQOM_NAMESPACE(GGMTree_GetNextLeaf_adv)
#define GGMTree_InitIncrementalPartialExpansion_adv MQOM_NAMESPACE(GGMTree_InitIncrementalPartialExpansion_adv)
#define GGMTree_GetNextLeafPartial_adv MQOM_NAMESPACE(GGMTree_GetNextLeafPartial_adv)

#ifndef GGMTREE_NB_PARALLEL_DERIVATIONS_LOG
#define GGMTREE_NB_PARALLEL_DERIVATIONS_LOG 0
#endif

#if GGMTREE_NB_PARALLEL_DERIVATIONS_LOG >= MQOM2_PARAM_NB_EVALS_LOG
#error GGMTREE_NB_PARALLEL_DERIVATIONS_LOG should be smaller than MQOM2_PARAM_NB_EVALS_LOG
#endif

#define GGMTREE_NB_PARALLEL_DERIVATIONS (1<<GGMTREE_NB_PARALLEL_DERIVATIONS_LOG)

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint32_t e;
	const uint8_t* salt;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_ecb *ctx_enc;
#else
	enc_ctx_ecb ctx_enc[GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG][2*GGMTREE_NB_PARALLEL_DERIVATIONS][MQOM2_PARAM_SEED_SIZE];
} ggmtree_ctx_adv_t;

/* Cleaning function */
static inline void ggmtree_ctx_adv_t_clean(ggmtree_ctx_adv_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_ecb(&ctx->ctx_enc[i]);
	}
#endif
}

int GGMTree_InitIncrementalExpansion_adv(ggmtree_ctx_adv_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e);

int GGMTree_GetNextLeaf_adv(ggmtree_ctx_adv_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]);

typedef struct {
	uint32_t active;
	uint32_t num_leaf;
	uint32_t e;
	const uint8_t* salt;
	uint32_t i_star;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY == 0
	/* Dummy value, not used */
	enc_ctx_pub_ecb *ctx_enc;
#else
	enc_ctx_pub_ecb ctx_enc[GGMTREE_NB_ENC_CTX_IN_MEMORY];
#endif
	uint8_t path[MQOM2_PARAM_NB_EVALS_LOG-GGMTREE_NB_PARALLEL_DERIVATIONS_LOG][2*GGMTREE_NB_PARALLEL_DERIVATIONS][MQOM2_PARAM_SEED_SIZE];
	uint8_t opening[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE];
} ggmtree_ctx_partial_adv_t;

/* Cleaning function */
static inline void ggmtree_ctx_partial_adv_t_clean(ggmtree_ctx_partial_adv_t *ctx) {
	(void)ctx;
#if GGMTREE_NB_ENC_CTX_IN_MEMORY > 0
	unsigned int i;
	for(i = 0; i < GGMTREE_NB_ENC_CTX_IN_MEMORY; i++){
		enc_clean_ctx_pub_ecb(&ctx->ctx_enc[i]);
	}
#endif
}

int GGMTree_InitIncrementalPartialExpansion_adv(ggmtree_ctx_partial_adv_t* ctx, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star);

int GGMTree_GetNextLeafPartial_adv(ggmtree_ctx_partial_adv_t* ctx, uint8_t lseed[MQOM2_PARAM_SEED_SIZE]);


#endif /* __GGM_TREE_INCR_X1_ADV_H__ */
