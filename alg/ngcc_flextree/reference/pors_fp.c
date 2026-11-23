#include "pors_fp.h"
#include "address.h"
#include "hash.h"
#include "params/pors_precomputed-flextree-160f.h"

#include "ngcc/auxfunc.h"
#include "tree/merkle.h"
#include "util/bytes.h"

#include <string.h>

#define PORS_AUTH_NODE_CAPACITY ((size_t)SPX_PORS_FP_K * PORS_TREE_HEIGHT_CONST)

typedef struct {
    const spx_ctx *ctx;
    uint32_t secret_addr[8];
    uint32_t leaf_addr[8];
    const uint32_t *selected;
    uint8_t *secrets;
    size_t next_selected;
} pors_leaf_context;

/* Bind the generic tree builder to PORS secret and leaf address domains. */
static void generate_leaf(uint8_t *out, uint32_t index, void *data)
{
    pors_leaf_context *state = data;
    uint8_t secret[SPX_N];

    set_tree_index(state->secret_addr, index);
    prf_addr(secret, state->ctx, state->secret_addr);
    if (state->next_selected < SPX_PORS_FP_K &&
        state->selected[state->next_selected] == index) {
        memcpy(state->secrets + state->next_selected * SPX_N, secret, SPX_N);
        ++state->next_selected;
    }
    set_tree_index(state->leaf_addr, index);
    thash(out, secret, 1, state->ctx, state->leaf_addr);
}

static int message_to_indices(uint32_t indices[SPX_PORS_FP_K],
                               const uint8_t *randomizer, const uint8_t *message,
                               uint32_t counter)
{
    uint8_t input[SPX_N + COUNTER_SIZE + SPX_PORS_FP_MSG_BYTES];
    uint8_t digest[PORS_HPORS_INPUT_BYTES];
    pqc_u256 rank;

    memcpy(input, randomizer, SPX_N);
    pqc_store_be(input + SPX_N, COUNTER_SIZE, counter);
    memcpy(input + SPX_N + COUNTER_SIZE, message, SPX_PORS_FP_MSG_BYTES);
    if (pseudoXOF(sizeof(digest) * 8, input, sizeof(input) * 8, digest) != 0)
        return -1;
    pqc_u256_reduce_bytes_table(&rank, digest, sizeof(digest),
                                pors_hpors_byte_mod_table,
                                &pors_total_combination);
    return pqc_combination_unrank(indices, SPX_PORS_FP_T, SPX_PORS_FP_K,
                                  &rank, pors_binom_tables);
}

void pors_fp_sign(unsigned char *sig, unsigned char *pk,
                  const unsigned char *message, const spx_ctx *ctx,
                  const uint32_t fors_addr[8], const unsigned char *randomizer)
{
    pqc_merkle_node auth_nodes[PORS_AUTH_NODE_CAPACITY];
    uint32_t selected[SPX_PORS_FP_K];
    size_t auth_count = 0;
    uint32_t counter = 0;

    for (;;) {
        if (counter == UINT32_MAX ||
            message_to_indices(selected, randomizer, message, counter) != 0 ||
            pqc_merkle_auth_nodes(auth_nodes, &auth_count, PORS_AUTH_NODE_CAPACITY,
                                  selected, SPX_PORS_FP_K, SPX_PORS_FP_T) != 0)
            goto failed;
        if (auth_count <= SPX_PORS_FP_MAX_AUTH_NODES)
            break;
        ++counter;
    }

    pqc_store_be(sig, COUNTER_SIZE, counter);
    uint8_t *secrets = sig + COUNTER_SIZE;
    uint8_t *auth_values = secrets + SPX_PORS_FP_K * SPX_N;
    pors_leaf_context leaves = {.ctx = ctx, .selected = selected, .secrets = secrets};
    uint32_t node_addr[8] = {0};
    flextree_tree_context tree = {ctx, node_addr};

    copy_keypair_addr(leaves.secret_addr, fors_addr);
    copy_keypair_addr(leaves.leaf_addr, fors_addr);
    copy_keypair_addr(node_addr, fors_addr);
    set_type(leaves.secret_addr, SPX_ADDR_TYPE_FORSPRF);
    set_type(leaves.leaf_addr, SPX_ADDR_TYPE_FORSTREE);
    set_type(node_addr, SPX_ADDR_TYPE_FORSTREE);
    if (pqc_merkle_build(pk, auth_values, SPX_N, SPX_PORS_FP_T, auth_nodes, auth_count,
                          generate_leaf, &leaves, flextree_hash_tree_node, &tree) != 0 ||
        leaves.next_selected != SPX_PORS_FP_K)
        goto failed;
    memset(auth_values + auth_count * SPX_N, 0,
           (SPX_PORS_FP_MAX_AUTH_NODES - auth_count) * SPX_N);
    return;

failed:
    memset(sig, 0, SPX_PORS_FP_BYTES);
    memset(pk, 0, SPX_N);
}

void pors_fp_pk_from_sig(unsigned char *pk, const unsigned char *sig,
                         const unsigned char *message, const spx_ctx *ctx,
                         const uint32_t fors_addr[8], const unsigned char *randomizer)
{
    pqc_merkle_node auth_nodes[PORS_AUTH_NODE_CAPACITY];
    uint32_t selected[SPX_PORS_FP_K];
    uint8_t leaf_values[SPX_PORS_FP_K * SPX_N];
    const uint32_t counter = (uint32_t)pqc_load_be(sig, COUNTER_SIZE);
    size_t auth_count = 0;

    if (counter == UINT32_MAX ||
        message_to_indices(selected, randomizer, message, counter) != 0 ||
        pqc_merkle_auth_nodes(auth_nodes, &auth_count, PORS_AUTH_NODE_CAPACITY,
                              selected, SPX_PORS_FP_K, SPX_PORS_FP_T) != 0 ||
        auth_count > SPX_PORS_FP_MAX_AUTH_NODES)
        goto failed;

    const uint8_t *secrets = sig + COUNTER_SIZE;
    const uint8_t *auth_values = secrets + SPX_PORS_FP_K * SPX_N;
    /* The signature reserves a fixed authentication region. Its unused suffix
     * must match the zero padding produced by the signer. */
    for (size_t i = auth_count * SPX_N; i < SPX_PORS_FP_MAX_AUTH_NODES * SPX_N; ++i)
        if (auth_values[i] != 0)
            goto failed;

    uint32_t leaf_addr[8] = {0}, node_addr[8] = {0};
    flextree_tree_context tree = {ctx, node_addr};
    copy_keypair_addr(leaf_addr, fors_addr);
    copy_keypair_addr(node_addr, fors_addr);
    set_type(leaf_addr, SPX_ADDR_TYPE_FORSTREE);
    set_type(node_addr, SPX_ADDR_TYPE_FORSTREE);
    for (size_t i = 0; i < SPX_PORS_FP_K; ++i) {
        set_tree_index(leaf_addr, selected[i]);
        thash(leaf_values + i * SPX_N, secrets + i * SPX_N, 1, ctx, leaf_addr);
    }
    if (pqc_merkle_root_from_proof(pk, selected, leaf_values, SPX_PORS_FP_K,
                                   auth_values, auth_count, SPX_N, SPX_PORS_FP_T,
                                   flextree_hash_tree_node, &tree) == 0)
        return;
failed:
    memset(pk, 0, SPX_N);
}
