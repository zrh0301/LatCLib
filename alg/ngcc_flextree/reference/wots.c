#include "util/bytes.h"
#include <stdint.h>
#include <string.h>

#include "hash.h"
#include "wots.h"
#include "address.h"
#include "params.h"

const unsigned int wots_w[]    = SPX_WOTS_W_ARRAY;
const unsigned int wots_wlog[] = SPX_WOTS_LOGW_ARRAY;

/**
 * Computes the chaining function.
 * out and in have to be n-byte arrays.
 *
 * Interprets in as start-th value of the chain.
 * addr has to contain the address of the chain.
 */
static void gen_chain(unsigned char *out, const unsigned char *in,
                      unsigned int start, unsigned int steps,
                      const spx_ctx *ctx, uint32_t addr[8],uint32_t chain_idx)
{
    uint32_t i;

    /* Initialize out with the value at position 'start'. */
    memcpy(out, in, SPX_N);

    /* Iterate 'steps' calls to the hash function. */
    for (i = start; i < (start+steps) && i < wots_w[chain_idx]; i++) {
        set_hash_addr(addr, i);
        thash(out, out, 1, ctx, addr);
    }
}

/**
 * Takes a message and derives the matching chain lengths.
 */
unsigned int chain_lengths(unsigned int *lengths, const unsigned char *msg)
{
    uint32_t bitbuf = 0;
    unsigned int bits_in_buf = 0;
    unsigned int csum = 0;
    const unsigned char *in = msg;

    for (unsigned int i = 0; i < SPX_WOTS_LEN; i++) {
        while (bits_in_buf < wots_wlog[i]) {
            bitbuf = (bitbuf << 8) | *in++;
            bits_in_buf += 8;
        }
        bits_in_buf -= wots_wlog[i];
        lengths[i] = (bitbuf >> bits_in_buf) & (wots_w[i] - 1u);
        /* sign-11-8: the specification fixes the digit sum, not its complement. */
        csum += lengths[i];
    }

    return csum;
}

/*
 * sign-11-7: the reference verifier tested the high WOTS_ZERO_BITS bits of the
 * final digest byte, a window that overlaps the bits chain_lengths consumes on
 * seven parameter sets. Builds with PQC_NGCC_FLEXTREE_STRICT check exactly the
 * trailing bits left over by the digit decoder. The default build keeps the
 * submitted predicate so that the published KAT stays reproducible.
 */
int wots_digest_zero_bits_ok(const unsigned char *digest)
{
#ifdef PQC_NGCC_FLEXTREE_STRICT
    unsigned int consumed = 0;
    unsigned int bytes;
    unsigned int i;

    for (i = 0; i < SPX_WOTS_LEN; i++)
        consumed += wots_wlog[i];
    if (consumed > 8u * SPX_N)
        return 0;

    bytes = consumed / 8u;
    if ((consumed & 7u) != 0u) {
        unsigned int remain = 8u - (consumed & 7u);
        if ((digest[bytes] & (unsigned char)(0xFFu >> (8u - remain))) != 0)
            return 0;
        bytes++;
    }
    for (; bytes < SPX_N; bytes++)
        if (digest[bytes] != 0)
            return 0;
    return 1;
#else
    const uint32_t mask = (~0U << (8 - WOTS_ZERO_BITS));
    return (digest[SPX_N - 1] & mask) == 0;
#endif
}

/**
 * Takes a WOTS signature and an n-byte message, computes a WOTS public key.
 *
 * Writes the computed public key to 'pk'.
 */
void wots_pk_from_sig(unsigned char *pk,
                      const unsigned char *sig, const unsigned char *msg,
                      const spx_ctx *ctx, uint32_t addr[8], uint32_t counter)
{
    unsigned int lengths[SPX_WOTS_LEN];
    uint32_t i;

    /*Initial parameters for validation of checksum*/
    int csum;
    unsigned char digest[SPX_N];

    /*Set thash address for custom hash to type 6 & PK format*/
    uint32_t wots_pk_addr[8] = {0};
    set_type(wots_pk_addr, SPX_ADDR_TYPE_COMPRESS_WOTS);
    copy_keypair_addr(wots_pk_addr, addr);

    /*Set padding*/
    pqc_store_be(((unsigned char *) (wots_pk_addr))+(SPX_OFFSET_COUNTER) , COUNTER_SIZE, counter);
    /*Calculate checksum*/
    thash(digest, msg, 1, ctx, wots_pk_addr);
    

    csum = chain_lengths(lengths, digest);

    /*Validate Checksum*/
    if ((csum != WANTED_CHECKSUM) || !wots_digest_zero_bits_ok(digest)){   
        /* The whole WOTS public key is read by the caller, so clear all of it
         * (the submitted code cleared only SPX_PK_BYTES, leaving most of the
         * buffer uninitialised). */
        memset(pk,0,SPX_WOTS_PK_BYTES);
    }
    else
    {
    for (i = 0; i < SPX_WOTS_LEN; i++) {
        set_chain_addr(addr, i);
        gen_chain(pk + i*SPX_N, sig + i*SPX_N,
                  lengths[i], wots_w[i] - 1 - lengths[i], ctx, addr,i);
        }
    }
}

/*
 * This generates a WOTS public key
 * It also generates the WOTS signature if leaf_info indicates
 * that we're signing with this WOTS key
 */
void wots_gen_leaf(unsigned char *dest, uint32_t leaf_idx, void *v_info) {
    wots_leaf_context *info = v_info;
    const spx_ctx *ctx = info->ctx;
    uint32_t *leaf_addr = info->leaf_addr;
    uint32_t *pk_addr = info->pk_addr;
    unsigned int i, k;
    unsigned char pk_buffer[ SPX_WOTS_BYTES ];
    unsigned char *buffer;
    uint32_t wots_k_mask;

    if (leaf_idx == info->wots_sign_leaf) {
        /* We're traversing the leaf that's signing; generate the WOTS */
        /* signature */
        wots_k_mask = 0;
    } else {
        /* Nope, we're just generating pk's; turn off the signature logic */
        wots_k_mask = ~0;
    }

    set_keypair_addr( leaf_addr, leaf_idx );
    set_keypair_addr( pk_addr, leaf_idx );

    for (i = 0, buffer = pk_buffer; i < SPX_WOTS_LEN; i++, buffer += SPX_N) {
        uint32_t wots_k = info->wots_steps[i] | wots_k_mask; /* Set wots_k to */
            /* the step if we're generating a signature, ~0 if we're not */

        /* Start with the secret seed */
        set_chain_addr(leaf_addr, i);
        set_hash_addr(leaf_addr, 0);
        set_type(leaf_addr, SPX_ADDR_TYPE_WOTSPRF);
 
        prf_addr(buffer, ctx, leaf_addr);

        set_type(leaf_addr, SPX_ADDR_TYPE_WOTS);

        /* Iterate down the WOTS chain */
        for (k=0;; k++) {
            /* Check if this is the value that needs to be saved as a */
            /* part of the WOTS signature */
            if (k == wots_k) {
                memcpy( info->wots_sig + i * SPX_N, buffer, SPX_N );
            }

            /* Check if we hit the top of the chain */
            if (k == wots_w[i] - 1) break;

            /* Iterate one step on the chain */
            set_hash_addr(leaf_addr, k);

            thash(buffer, buffer, 1, ctx, leaf_addr);
        }
    }

    /* Do the final thash to generate the public keys */
    thash(dest, pk_buffer, SPX_WOTS_LEN, ctx, pk_addr);
}
