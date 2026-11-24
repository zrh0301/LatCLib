#pragma once
#include <string.h>
#include <inttypes.h>
#include "tsuov_params.h"
#include "ngcc/auxfunc.h"

#define TSUOV_MGF_BLOCK_BYTES 32
#define TSUOV_MGF_POOLSIZE    TSUOV_MGF_BLOCK_BYTES

TYPEDEF_STRUCT (MGF_CTX,
  uint8_t  in     [TSUOV_SEED_LEN + 16] ;
  size_t   in_len                     ;
  uint8_t  pool     [TSUOV_MGF_BLOCK_BYTES] ;
  uint32_t pool_bytes                 ;
  uint32_t counter                    ;
) ;

extern MGF_CTX_s * MGF_init  (const uint8_t * seed, const size_t n0, MGF_CTX ctx) ;
extern MGF_CTX_s * MGF_update(const uint8_t * seed, const size_t n0, MGF_CTX ctx) ;
extern uint8_t   * MGF_yield (MGF_CTX ctx, uint8_t * dest, const size_t n1) ;
extern void        MGF_final (MGF_CTX ctx) ;
extern void        MGF_CTX_copy(MGF_CTX src, MGF_CTX dst) ;
