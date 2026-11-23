#ifndef SPX_HASH_H
#define SPX_HASH_H

#include <stdint.h>
#include "params.h"


void prf_addr(unsigned char *out, const spx_ctx *ctx,
              const uint32_t addr[8]);

int gen_message_random(unsigned char *R, const unsigned char *sk_prf,
                       const unsigned char *optrand,
                       const unsigned char *m, unsigned long long mlen,
                       const spx_ctx *ctx);

int hash_message(unsigned char *digest, uint64_t *tree, uint32_t *leaf_idx,
                  const unsigned char *R, const unsigned char *pk,
                  const unsigned char *m, unsigned long long mlen,
                  const spx_ctx *ctx);


void thash(unsigned char *out, const unsigned char *in, unsigned int inblocks,
           const spx_ctx *ctx, uint32_t addr[8]);

/* Protocol address binding for the shared binary-tree hash callback. */
typedef struct {
    const spx_ctx *ctx;
    uint32_t *addr;
} flextree_tree_context;

void flextree_hash_tree_node(uint8_t *out, const uint8_t *children,
                             uint32_t height, uint32_t index, void *data);

#endif
