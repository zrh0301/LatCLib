#ifndef PQC_MERKLE_H
#define PQC_MERKLE_H

#include <stddef.h>
#include <stdint.h>

/* Leaves are numbered from left to right. Parent height is measured upward
 * from the deepest leaves. children contains two node_bytes-byte children.
 * Hash callbacks may write out over either input child. */
typedef void (*pqc_merkle_leaf_fn)(uint8_t *out, uint32_t leaf, void *ctx);
typedef void (*pqc_merkle_hash_fn)(uint8_t *out, const uint8_t *children,
                                  uint32_t height, uint32_t index, void *ctx);

typedef struct {
    uint32_t height;
    uint32_t index;
} pqc_merkle_node;

/* Build a complete tree of 2^height leaves, with height < 32. An optional
 * authentication path selects leaf_index, or UINT32_MAX to skip the path.
 * offset is a multiple of 2^height and offsets both leaf and parent indices.
 * Return -1 for invalid dimensions or allocation failure. On failure output
 * is unchanged. root and auth_path must be disjoint. */
int pqc_merkle_treehash(uint8_t *root, uint8_t *auth_path, size_t node_bytes,
                        uint32_t height, uint32_t leaf_index, uint32_t offset,
                        pqc_merkle_leaf_fn leaf, void *leaf_ctx,
                        pqc_merkle_hash_fn hash, void *hash_ctx);

/* Recompute a complete-tree root from one leaf and its authentication path.
 * Exact aliasing root == leaf is allowed. Other buffers must be disjoint.
 * The offset and height conventions match pqc_merkle_treehash. */
int pqc_merkle_root_from_path(uint8_t *root, const uint8_t *leaf,
                             const uint8_t *auth_path, size_t node_bytes,
                             uint32_t height, uint32_t leaf_index, uint32_t offset,
                             pqc_merkle_hash_fn hash, void *hash_ctx);

/* Multiproofs also support trees with a shortened bottom level. With L leaves
 * and B the largest power of two strictly below L, the first 2*(L-B) leaves
 * sit at height 0. The remaining leaves sit at height 1. A power-of-two L
 * gives a complete tree. L == 1 is a single leaf at height 0.
 *
 * selected contains count strictly increasing indices below leaf_count.
 * Authentication nodes are ordered by increasing height, then index.
 * auth_capacity must hold count*ceil(log2(leaf_count)) nodes. On failure
 * auth_count is unchanged. Output may contain partial entries. */
int pqc_merkle_auth_nodes(pqc_merkle_node *out, size_t *auth_count,
                          size_t auth_capacity, const uint32_t *selected,
                          size_t count, uint32_t leaf_count);

/* Build the tree and return the root and values at the supplied authentication
 * positions. auth_values has auth_count*node_bytes bytes. leaf is called once
 * per leaf in index order. Return -1 for invalid dimensions, positions or
 * allocation failure. On failure outputs are unchanged. */
int pqc_merkle_build(uint8_t *root, uint8_t *auth_values, size_t node_bytes,
                     uint32_t leaf_count, const pqc_merkle_node *auth,
                     size_t auth_count, pqc_merkle_leaf_fn leaf, void *leaf_ctx,
                     pqc_merkle_hash_fn hash, void *hash_ctx);

/* Rebuild a root from already-hashed selected leaves and their multiproof.
 * leaves has count*node_bytes bytes, in selected-index order. Authentication
 * values follow pqc_merkle_auth_nodes order. Missing or excess values fail.
 * Return -1 for invalid inputs or allocation failure. Output is unchanged on
 * failure. All input/output buffers must be disjoint. */
int pqc_merkle_root_from_proof(uint8_t *root, const uint32_t *selected,
                              const uint8_t *leaves, size_t count,
                              const uint8_t *auth_values, size_t auth_count,
                              size_t node_bytes, uint32_t leaf_count,
                              pqc_merkle_hash_fn hash, void *hash_ctx);

#endif
