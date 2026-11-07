#ifndef __GGM_TREE_DEFAULT_H__
#define __GGM_TREE_DEFAULT_H__

/* MQOM2 parameters */
#include "mqom2_parameters.h"
/* Encryption primitive */
#include "enc.h"
/* Common helpers */
#include "common.h"

/* Deal with namespacing */
#define GGMTree_Expand MQOM_NAMESPACE(GGMTree_Expand)
#define GGMTree_Open MQOM_NAMESPACE(GGMTree_Open)
#define GGMTree_PartiallyExpand MQOM_NAMESPACE(GGMTree_PartiallyExpand)
#define GGMTree_ExpandPath MQOM_NAMESPACE(GGMTree_ExpandPath)

int GGMTree_Expand(const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e, uint8_t node[MQOM2_PARAM_FULL_TREE_SIZE + 1][MQOM2_PARAM_SEED_SIZE], uint8_t lseed[MQOM2_PARAM_NB_EVALS][MQOM2_PARAM_SEED_SIZE]);

int GGMTree_Open(const uint8_t node[MQOM2_PARAM_FULL_TREE_SIZE + 1][MQOM2_PARAM_SEED_SIZE], uint32_t i_star, uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE]);

int GGMTree_PartiallyExpand(const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star, uint8_t lseed[MQOM2_PARAM_NB_EVALS][MQOM2_PARAM_SEED_SIZE]);

int GGMTree_ExpandPath(const uint8_t salt[MQOM2_PARAM_SALT_SIZE], const uint8_t rseed[MQOM2_PARAM_SEED_SIZE], const uint8_t delta[MQOM2_PARAM_SEED_SIZE], uint32_t e, uint32_t i_star, uint8_t path[MQOM2_PARAM_NB_EVALS_LOG][MQOM2_PARAM_SEED_SIZE], uint8_t lseed[MQOM2_PARAM_SEED_SIZE]);

#endif /* __GGM_TREE_DEFAULT_H__ */
