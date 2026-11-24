#pragma once

#include "tsuov_params.h"
#include "stdint.h"
#include <stddef.h>

#include "gf/gf31.h"
#include "gf/gf31_matrix.h"
#include "util/bits.h"

typedef gf31_t Fq ;

/* 域运算与 L=2 塔算子都来自库内 inc/gf/gf31.h 与 inc/gf/gf31_matrix.h。
 * 这里只保留参考实现使用的名字，方便算法层继续按原符号书写。 */
inline static int Fq_reduction(int Z){ return (int) gf31_reduce(Z) ; }
inline static Fq Fq_add(Fq X, Fq Y){ return gf31_add(X, Y) ; }
inline static Fq Fq_sub(Fq X, Fq Y){ return gf31_sub(X, Y) ; }
inline static Fq Fq_mul(Fq X, Fq Y){ return gf31_mul(X, Y) ; }
inline static Fq Fq_inv(Fq X){ return gf31_inv(X) ; }

#define phi_t2_q31_planar gf31_t2_phi_planar

/* Cursor-based wrappers over the shared LSB-first bit packer. The destination
 * must be zeroed before the first store (util/bits.h ORs the bits in). */
inline static void store_bits(
  int             x,
  const int       num_bits,
  uint8_t       * pool,
  size_t        * pool_bits
){
  pqc_store_bits_le(pool, *pool_bits, (uint32_t) x, (unsigned int) num_bits) ;
  *pool_bits += (size_t) num_bits ;
}

inline static int restore_bits(
  const uint8_t * pool,
  size_t        * pool_bits,
  const int       num_bits
){
  uint32_t x = pqc_load_bits_le(pool, *pool_bits, (unsigned int) num_bits) ;
  *pool_bits += (size_t) num_bits ;
  return (int) x ;
}

inline static void store_Fq(
  int             x,
  uint8_t       * pool,
  size_t        * pool_bits
){
  store_bits(x, TSUOV_ceil_log_2_q, pool, pool_bits) ;
}

inline static int restore_Fq(
  const uint8_t * pool,
  size_t * pool_bits,
  Fq * out
){
  int x = restore_bits(pool, pool_bits, TSUOV_ceil_log_2_q) ;
  /* The unused 5-bit code 31 must be rejected: otherwise it decodes to zero
   * during arithmetic and creates signature/public-key aliases (sign-31-6). */
  if((unsigned int) x >= TSUOV_q)
    return -1 ;
  *out = (Fq) x ;
  return 0 ;
}
