#include <string.h>
#include <stdlib.h>
#include "mgf.h"

static void mgf_refresh(MGF_CTX ctx)
{
  unsigned int c = ctx->counter++ ;
  ctx->in[ctx->in_len + 0] = (uint8_t)(c >> 24) ;
  ctx->in[ctx->in_len + 1] = (uint8_t)(c >> 16) ;
  ctx->in[ctx->in_len + 2] = (uint8_t)(c >>  8) ;
  ctx->in[ctx->in_len + 3] = (uint8_t)(c      ) ;
  sm3hash(256, ctx->in, (unsigned long long)(ctx->in_len + 4) * 8, ctx->pool) ;
  ctx->pool_bytes = TSUOV_MGF_BLOCK_BYTES ;
}

MGF_CTX_s * MGF_init(const uint8_t * seed, const size_t n0, MGF_CTX ctx)
{
  ctx->counter    = 0 ;
  ctx->pool_bytes = 0 ;
  memcpy(ctx->in, seed, n0) ;
  ctx->in_len = n0 ;
  return ctx ;
}

MGF_CTX_s * MGF_update(const uint8_t * seed, const size_t n0, MGF_CTX ctx)
{
  memcpy(ctx->in + ctx->in_len, seed, n0) ;
  ctx->in_len += n0 ;
  ctx->counter    = 0 ;
  ctx->pool_bytes = 0 ;
  return ctx ;
}

uint8_t * MGF_yield(MGF_CTX ctx, uint8_t * dest, const size_t n1)
{
  for(size_t i = 0; i < n1; i++){
    if(ctx->pool_bytes == 0) mgf_refresh(ctx) ;
    dest[i] = ctx->pool[TSUOV_MGF_BLOCK_BYTES - ctx->pool_bytes] ;
    ctx->pool_bytes-- ;
  }
  return dest ;
}

void MGF_final(MGF_CTX ctx)
{
  memset(ctx->pool, 0, sizeof(ctx->pool)) ;
  memset(ctx->in, 0, sizeof(ctx->in)) ;
}

void MGF_CTX_copy(MGF_CTX src, MGF_CTX dst)
{
  memcpy(dst, src, sizeof(MGF_CTX)) ;
}
