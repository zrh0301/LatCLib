#ifndef SPX_WOTS_H
#define SPX_WOTS_H

#include <stdint.h>

#include "params.h"


extern const unsigned int wots_w[];
extern const unsigned int wots_wlog[];

/**
 * Takes a WOTS signature and an n-byte message, computes a WOTS public key.
 *
 * Writes the computed public key to 'pk'.
 */
void wots_pk_from_sig(unsigned char *pk,
                      const unsigned char *sig, const unsigned char *msg,
                      const spx_ctx *ctx, uint32_t addr[8], uint32_t counter);

/*
 * Compute the chain lengths needed for a given message hash
 */
unsigned int chain_lengths(unsigned int *lengths, const unsigned char *msg);

/*
 * Check the OTS digest bits that chain_lengths leaves unused (sign-11-7).
 */
int wots_digest_zero_bits_ok(const unsigned char *digest);


typedef struct {
    const spx_ctx *ctx;
    unsigned char *wots_sig;
    uint32_t wots_sign_leaf;
    uint32_t *wots_steps;
    uint32_t leaf_addr[8];
    uint32_t pk_addr[8];
} wots_leaf_context;

void wots_gen_leaf(unsigned char *dest, uint32_t leaf_idx, void *data);

#endif
