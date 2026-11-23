#ifndef SPX_PARAMS_H
#define SPX_PARAMS_H

#include <stdint.h>

/* Hash output length in bytes. */
#define SPX_N 20
/* Height of the hypertree. */
#define SPX_FULL_HEIGHT 67
/* Number of subtree layer. */
#define SPX_D 17
/* PORS+FP dimensions. */
#define SPX_PORS_FP_T 4721
#define SPX_PORS_FP_K 27
/* Winternitz parameter, */

#define SPX_WOTS_W_ARRAY {8, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16}
#define SPX_WOTS_LOGW_ARRAY {3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}

#define SPX_WOTS_LEN 40
#define WOTS_ZERO_BITS 1
#define WANTED_CHECKSUM  296  /*(SUM_W - SPX_WOTS_LEN) / 2)*/

/* For clarity */
#define SPX_ADDR_BYTES 32

#define SPX_WOTS_BYTES (SPX_WOTS_LEN * SPX_N)
#define SPX_WOTS_PK_BYTES SPX_WOTS_BYTES

/* Subtree size. */
#define SPX_TREE_HEIGHT (SPX_FULL_HEIGHT / SPX_D)

#if SPX_TREE_HEIGHT * SPX_D != SPX_FULL_HEIGHT
    #define SPX_BOTTOM_TREE_HEIGHT  SPX_TREE_HEIGHT + 1
#else
    #define SPX_BOTTOM_TREE_HEIGHT  SPX_TREE_HEIGHT
#endif

/* PORS+FP parameters. */
#define SPX_PORS_FP_MSG_BYTES SPX_N
#define SPX_PORS_FP_MAX_AUTH_NODES 155
#define SPX_PORS_FP_BYTES ((COUNTER_SIZE + (SPX_PORS_FP_K + SPX_PORS_FP_MAX_AUTH_NODES) * SPX_N))
#define SPX_PORS_FP_PK_BYTES SPX_N

/* Resulting SPX sizes. */
#define SPX_BYTES (SPX_N + SPX_PORS_FP_BYTES + SPX_D * SPX_WOTS_BYTES +\
                   SPX_FULL_HEIGHT * SPX_N + (SPX_D * COUNTER_SIZE))
#define SPX_PK_BYTES (2 * SPX_N)
#define SPX_SK_BYTES (2 * SPX_N + SPX_PK_BYTES)


/* Signature search counter width. */
#define COUNTER_SIZE 4

typedef struct {
    uint8_t pub_seed[SPX_N];
    uint8_t sk_seed[SPX_N];
} spx_ctx;

#endif
