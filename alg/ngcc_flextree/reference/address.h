#ifndef SPX_ADDRESS_H
#define SPX_ADDRESS_H

#include <stdint.h>

#define SPX_OFFSET_LAYER     3   /* The byte used to specify the Merkle tree layer */
#define SPX_OFFSET_TREE      8   /* The start of the 8 byte field used to specify the tree */
#define SPX_OFFSET_TYPE      19  /* The byte used to specify the hash type (reason) */
#define SPX_OFFSET_KP_ADDR2  22  /* The high byte used to specify the key pair (which one-time signature) */
#define SPX_OFFSET_KP_ADDR1  23  /* The low byte used to specify the key pair */
#define SPX_OFFSET_CHAIN_ADDR 27  /* The byte used to specify the chain address (which Winternitz chain) */
#define SPX_OFFSET_HASH_ADDR 31  /* The byte used to specify the hash address (where in the Winternitz chain) */
#define SPX_OFFSET_TREE_HGT  27  /* The byte used to specify the height of this node in the FORS or Merkle tree */
#define SPX_OFFSET_TREE_INDEX 28 /* The start of the 4 byte field used to specify the node in the FORS or Merkle tree */
#define SPX_OFFSET_COUNTER 24 /* The start of the 4 byte field used to specify the counter for the WOTS signature */


/* The hash types that are passed to set_type */
#define SPX_ADDR_TYPE_WOTS 0 
#define SPX_ADDR_TYPE_WOTSPK 1
#define SPX_ADDR_TYPE_HASHTREE 2
#define SPX_ADDR_TYPE_FORSTREE 3
#define SPX_ADDR_TYPE_WOTSPRF 4
#define SPX_ADDR_TYPE_FORSPRF 5
#define SPX_ADDR_TYPE_COMPRESS_WOTS 6

void set_layer_addr(uint32_t addr[8], uint32_t layer);

void set_tree_addr(uint32_t addr[8], uint64_t tree);

void set_type(uint32_t addr[8], uint32_t type);

/* Copies the layer and tree part of one address into the other */
void copy_subtree_addr(uint32_t out[8], const uint32_t in[8]);

/* These functions are used for WOTS and FORS addresses. */

void set_keypair_addr(uint32_t addr[8], uint32_t keypair);

void set_chain_addr(uint32_t addr[8], uint32_t chain);

void set_hash_addr(uint32_t addr[8], uint32_t hash);

void copy_keypair_addr(uint32_t out[8], const uint32_t in[8]);

/* These functions are used for all hash tree addresses (including FORS). */

void set_tree_height(uint32_t addr[8], uint32_t tree_height);

void set_tree_index(uint32_t addr[8], uint32_t tree_index);

#endif
