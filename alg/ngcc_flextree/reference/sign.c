#include "tree/merkle.h"
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "api.h"
#include "params.h"
#include "wots.h"
#include "pors_fp.h"
#include "hash.h"
#include "address.h"
#include "ngcc/drng.h"
#include "merkle.h"
#include "util/ct.h"
#include "util/bytes.h"

extern DRNG_ctx drng_algorithm;

/*
 * Returns the length of a secret key, in bytes
 */
unsigned long long crypto_sign_secretkeybytes(void)
{
    return CRYPTO_SECRETKEYBYTES;
}

/*
 * Returns the length of a public key, in bytes
 */
unsigned long long crypto_sign_publickeybytes(void)
{
    return CRYPTO_PUBLICKEYBYTES;
}

/*
 * Returns the length of a signature, in bytes
 */
unsigned long long crypto_sign_bytes(void)
{
    return CRYPTO_BYTES;
}

/*
 * Returns the length of the seed required to generate a key pair, in bytes
 */
unsigned long long crypto_sign_seedbytes(void)
{
    return CRYPTO_SEEDBYTES;
}

/*
 * Generates an SPX key pair given a seed of length
 * Format sk: [SK_SEED || SK_PRF || PUB_SEED || root]
 * Format pk: [PUB_SEED || root]
 */
int crypto_sign_seed_keypair(unsigned char *pk, unsigned char *sk,
                             const unsigned char *seed)
{
    spx_ctx ctx;

    /* Initialize SK_SEED, SK_PRF and PUB_SEED from seed. */
    memcpy(sk, seed, CRYPTO_SEEDBYTES);

    memcpy(pk, sk + 2*SPX_N, SPX_N);

    memcpy(ctx.pub_seed, pk, SPX_N);
    memcpy(ctx.sk_seed, sk, SPX_N);

    /* Compute root node of the top-most subtree. */
    merkle_gen_root(sk + 3*SPX_N, &ctx);

    memcpy(pk + SPX_N, sk + 3*SPX_N, SPX_N);

    return 0;
}

/*
 * Generates an SPX key pair.
 * Format sk: [SK_SEED || SK_PRF || PUB_SEED || root]
 * Format pk: [PUB_SEED || root]
 */
int crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
  unsigned char seed[CRYPTO_SEEDBYTES];
  if (get_random_number(&drng_algorithm, seed, sizeof(seed) * 8) != 0)
    return -1;
  return crypto_sign_seed_keypair(pk, sk, seed);
}

/**
 * Returns an array containing a detached signature.
 */
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen, const uint8_t *sk)
{
    spx_ctx ctx;

    const unsigned char *sk_prf = sk + SPX_N;
    const unsigned char *pk = sk + 2*SPX_N;

    unsigned char optrand[SPX_N];
    unsigned char mhash[SPX_PORS_FP_MSG_BYTES];
    unsigned char root[SPX_N];
    unsigned long long i;
    uint64_t tree;
    uint32_t idx_leaf;
    uint32_t wots_addr[8] = {0};
    uint32_t tree_addr[8] = {0};
    uint32_t counter;
    const unsigned char *R;
    memcpy(ctx.sk_seed, sk, SPX_N);
    memcpy(ctx.pub_seed, pk, SPX_N);

    set_type(wots_addr, SPX_ADDR_TYPE_WOTS);
    set_type(tree_addr, SPX_ADDR_TYPE_HASHTREE);

    /* Optionally, signing can be made non-deterministic using optrand.
       This can help counter side-channel attacks that would benefit from
       getting a large number of traces when the signer uses the same nodes. */
    if (get_random_number(&drng_algorithm, optrand, sizeof(optrand) * 8) != 0)
        return -1;
    /* Compute the digest randomization value. */
    if (gen_message_random(sig, sk_prf, optrand, m, mlen, &ctx) != 0)
        return -1;
    R = sig;

    /* Derive the message digest and leaf index from R, PK and M. */
    if (hash_message(mhash, &tree, &idx_leaf, sig, pk, m, mlen, &ctx) == -1)
        return -1;
    sig += SPX_N;

    set_tree_addr(wots_addr, tree);
    set_keypair_addr(wots_addr, idx_leaf);

    /* Sign the message hash using PORS+FP. */
    pors_fp_sign(sig, root, mhash, &ctx, wots_addr, R);

    sig += SPX_PORS_FP_BYTES;
    
    uint32_t merkle_tree_heights[SPX_D];
    uint32_t num_k = SPX_FULL_HEIGHT - SPX_TREE_HEIGHT*SPX_D;
    for (uint32_t j = 0; j < num_k; j++) {
        merkle_tree_heights[j] = SPX_TREE_HEIGHT+1;
    }
    for (uint32_t j = num_k; j < SPX_D; j++) {
        merkle_tree_heights[j] = SPX_TREE_HEIGHT;
    }

    for (i = 0; i < SPX_D; i++) {
        set_layer_addr(tree_addr, i);
        set_tree_addr(tree_addr, tree);

        copy_subtree_addr(wots_addr, tree_addr);
        set_keypair_addr(wots_addr, idx_leaf);

        merkle_sign(sig, root, &ctx, wots_addr, tree_addr, idx_leaf, &counter,merkle_tree_heights[i]);
        if (counter == 0)
            return -1;
        pqc_store_be(sig + SPX_WOTS_BYTES + merkle_tree_heights[i] * SPX_N,
                     COUNTER_SIZE, counter);
        sig += SPX_WOTS_BYTES + merkle_tree_heights[i] * SPX_N + COUNTER_SIZE;
        /* Update the indices for the next layer. */
        if (i + 1 < SPX_D) {
            idx_leaf = (tree & ((1 << merkle_tree_heights[i + 1]) - 1));
            tree = tree >> merkle_tree_heights[i + 1];
        }
    }

    *siglen = SPX_BYTES;
    return 0;
}

/**
 * Verifies a detached signature and message under a given public key.
 */
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen, const uint8_t *pk)
{
    spx_ctx ctx;
    const unsigned char *pub_root = pk + SPX_N;
    unsigned char mhash[SPX_PORS_FP_MSG_BYTES];
    unsigned char wots_pk[SPX_WOTS_BYTES];
    unsigned char root[SPX_N];
    unsigned char leaf[SPX_N];
    unsigned int i;
    uint64_t tree;
    uint32_t idx_leaf;
    uint32_t wots_addr[8] = {0};
    uint32_t tree_addr[8] = {0};
    uint32_t wots_pk_addr[8] = {0};
    uint32_t counter;
    const unsigned char *R = sig;

    if (siglen != SPX_BYTES) {
        return -1;
    }
    
    memcpy(ctx.pub_seed, pk, SPX_N);

    set_type(wots_addr, SPX_ADDR_TYPE_WOTS);
    set_type(tree_addr, SPX_ADDR_TYPE_HASHTREE);
    set_type(wots_pk_addr, SPX_ADDR_TYPE_WOTSPK);

    /* Derive the message digest and leaf index from R || PK || M. */
    if (hash_message(mhash, &tree, &idx_leaf, sig, pk, m, mlen, &ctx) == -1)
        return -1;
    sig += SPX_N;

    /* Layer correctly defaults to 0, so no need to set_layer_addr */
    set_tree_addr(wots_addr, tree);
    set_keypair_addr(wots_addr, idx_leaf);
    pors_fp_pk_from_sig(root, sig, mhash, &ctx, wots_addr, R);
    sig += SPX_PORS_FP_BYTES;
    
    uint32_t merkle_tree_heights[SPX_D];
    uint32_t num_k = SPX_FULL_HEIGHT - SPX_TREE_HEIGHT*SPX_D;
    for (uint32_t j = 0; j < num_k; j++) {
        merkle_tree_heights[j] = SPX_TREE_HEIGHT+1;
    }
    for (uint32_t j = num_k; j < SPX_D; j++) {
        merkle_tree_heights[j] = SPX_TREE_HEIGHT;
    }

    /* For each subtree.. */
    for (i = 0; i < SPX_D; i++) {
        set_layer_addr(tree_addr, i);
        set_tree_addr(tree_addr, tree);

        copy_subtree_addr(wots_addr, tree_addr);
        set_keypair_addr(wots_addr, idx_leaf);

        copy_keypair_addr(wots_pk_addr, wots_addr);

        /* The WOTS public key is only correct if the signature was correct. */
        /* Initially, root is the PORS+FP pk, but on subsequent iterations it is
           the root of the subtree below the currently processed subtree. */
        counter = (uint32_t)pqc_load_be(
            sig + SPX_WOTS_BYTES + merkle_tree_heights[i] * SPX_N, COUNTER_SIZE);
        if (counter == 0)
            return -1;
        wots_pk_from_sig(wots_pk, sig, root, &ctx, wots_addr, counter);
        sig += SPX_WOTS_BYTES;

        /* Compute the leaf node using the WOTS public key. */
        thash(leaf, wots_pk, SPX_WOTS_LEN, &ctx, wots_pk_addr);
        /* Compute the root node of this subtree. */
        flextree_tree_context subtree = {&ctx, tree_addr};
        if (pqc_merkle_root_from_path(root, leaf, sig, SPX_N, merkle_tree_heights[i],
                                      idx_leaf, 0, flextree_hash_tree_node, &subtree) != 0)
            return -1;
        sig += merkle_tree_heights[i] * SPX_N  + COUNTER_SIZE;

        /* Update the indices for the next layer. */
        if (i + 1 < SPX_D) {
            idx_leaf = (tree & ((1 << merkle_tree_heights[i + 1]) - 1));
            tree = tree >> merkle_tree_heights[i + 1];
        }
    }

    /* Check if the root node equals the root node in the public key. */
    if (!pqc_ct_equal(root, pub_root, SPX_N)) {
        return -1;
    }

    return 0;
}


/**
 * Returns an array containing the signature followed by the message.
 */
int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk)
{
    size_t siglen;

    if (crypto_sign_signature(sm, &siglen, m, (size_t)mlen, sk) != 0)
        return -1;

    memmove(sm + SPX_BYTES, m, mlen);
    *smlen = siglen + mlen;

    return 0;
}

/**
 * Verifies a given signature-message pair under a given public key.
 */
int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk)
{
    /* The API caller does not necessarily know what size a signature should be
       but SPHINCS+ signatures are always exactly SPX_BYTES. */
    if (smlen < SPX_BYTES) {
        memset(m, 0, smlen);
        *mlen = 0;
        return -1;
    }

    *mlen = smlen - SPX_BYTES;

    if (crypto_sign_verify(sm, SPX_BYTES, sm + SPX_BYTES, (size_t)*mlen, pk)) {
        memset(m, 0, smlen);
        *mlen = 0;
        return -1;
    }

    /* If verification was successful, move the message to the right place. */
    memmove(m, sm + SPX_BYTES, *mlen);

    return 0;
}
