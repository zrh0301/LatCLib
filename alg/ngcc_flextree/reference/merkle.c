#include "tree/merkle.h"
#include "hash.h"
#include "util/bytes.h"
#include <stdint.h>
#include <string.h>

#include "wots.h"
#include "merkle.h"
#include "address.h"
#include "params.h"


/*
 * This generates a Merkle signature (WOTS signature followed by the Merkle
 * authentication path).  This is in this file because most of the complexity
 * is involved with the WOTS signature; the Merkle authentication path logic
 * uses the shared binary-tree traversal
 */ 
void merkle_sign(uint8_t *sig, unsigned char *root,
                 const spx_ctx *ctx,
                 uint32_t wots_addr[8], uint32_t tree_addr[8],
                 uint32_t idx_leaf, uint32_t *counter_out, uint32_t merkle_tree_height)
{
#define MAX_HASH_TRIALS_WOTS (1 << (20))
    unsigned char *auth_path = sig + SPX_WOTS_BYTES;
    wots_leaf_context info = {.ctx = ctx};
    unsigned steps[ SPX_WOTS_LEN ];

    /*Initial paramaters for custom thash & counter search*/
    unsigned char digest[SPX_N];
    uint32_t counter = 0;
    int csum;
    uint32_t to_sign = ~0;

    /*Initialize parameters for actual sign*/
    set_type(&tree_addr[0], SPX_ADDR_TYPE_HASHTREE);
    set_type(&info.pk_addr[0], SPX_ADDR_TYPE_WOTSPK);
    copy_subtree_addr(&info.leaf_addr[0], wots_addr);
    copy_subtree_addr(&info.pk_addr[0], wots_addr);

    /* Code for counter search */
    *counter_out = 0;
    if (idx_leaf != to_sign) {
        /*Set thash address for custom hash*/
        uint32_t *pk_addr = info.pk_addr;
        set_keypair_addr(pk_addr, idx_leaf);
        set_type(pk_addr, SPX_ADDR_TYPE_COMPRESS_WOTS);

        /*Search for correct counter */
        while (1) {
            counter++;
            if (counter > MAX_HASH_TRIALS_WOTS)
                return;
            pqc_store_be(((unsigned char *) (pk_addr))+(SPX_OFFSET_COUNTER) , COUNTER_SIZE, counter);
            thash(digest, root, 1, ctx, pk_addr);
            if (wots_digest_zero_bits_ok(digest))
            {
                csum = chain_lengths(steps, digest);
                if (csum == WANTED_CHECKSUM) 
                {
                    *counter_out = counter;
                    break;
                }
            }

        }

        /*Restore initial parameters for tree hash*/
        set_type(pk_addr, SPX_ADDR_TYPE_WOTSPK);
        pqc_store_be(((unsigned char *) (pk_addr))+(SPX_OFFSET_COUNTER) , COUNTER_SIZE, 0);
    }
    /* In this case we only try to generate the PK so no need to find the counter */
    else
    {
        chain_lengths(steps, root);
    }
    info.wots_sig = sig;
    info.wots_steps = steps;

    info.wots_sign_leaf = idx_leaf;

    flextree_tree_context tree = {ctx, tree_addr};
    if (pqc_merkle_treehash(root, auth_path, SPX_N, merkle_tree_height, idx_leaf, 0,
                            wots_gen_leaf, &info, flextree_hash_tree_node, &tree) != 0)
        memset(root, 0, SPX_N);
}

/* Compute root node of the top-most subtree. */
void merkle_gen_root(unsigned char *root, const spx_ctx *ctx)
{
    /* We do not need the auth path in key generation, but it simplifies the
       code to have just one treehash routine that computes both root and path
       in one function. */
    unsigned char auth_path[SPX_TREE_HEIGHT * SPX_N + SPX_WOTS_BYTES];
    uint32_t top_tree_addr[8] = {0};
    uint32_t wots_addr[8] = {0};
    uint32_t counter;

    set_layer_addr(top_tree_addr, SPX_D - 1);
    set_layer_addr(wots_addr, SPX_D - 1);

    merkle_sign(auth_path, root, ctx,
                wots_addr, top_tree_addr,
                ~0 /* ~0 means "don't bother generating an auth path */,
                &counter,SPX_TREE_HEIGHT);
}
