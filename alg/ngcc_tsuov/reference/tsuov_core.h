#pragma once

#include <stdlib.h>
#include <string.h>
#include "tsuov_params.h"
#include "Fql.h"
#include "mgf.h"
#include "ngcc/auxfunc.h"

#define TSUOV_PRG_CTX  MGF_CTX_s
#define PRG_init       PRG_init_sm3
#define PRG_yield      PRG_yield_sm3
#define PRG_final      PRG_final_sm3
#define PRG_copy       PRG_copy_sm3
#define RejSampPRG     RejSampPRG_sm3
#define TSUOV_PRG2_CTX MGF_CTX_s
#define PRG2_init      PRG_init_sm3
#define PRG2_yield     PRG_yield_sm3
#define PRG2_final     PRG_final_sm3
#define PRG2_copy      PRG_copy_sm3

#define SHIFT 5
#define BYTES2BLOCKS(BYTES)          (((BYTES)>>(SHIFT))+(((BYTES)&((1<<SHIFT)-1))?1:0))
#define BLOCK_ALIGNED_BYTES(BYTES)   (BYTES2BLOCKS(BYTES)<<(SHIFT))

typedef Fq         vector_V        [BLOCK_ALIGNED_BYTES(TSUOV_V)] ;
typedef Fq         vector_O        [BLOCK_ALIGNED_BYTES(TSUOV_O)] ;
typedef Fq         vector_m1       [BLOCK_ALIGNED_BYTES(TSUOV_m1)];
typedef Fq         vector_N        [BLOCK_ALIGNED_BYTES(TSUOV_N)] ;

typedef vector_V   VECTOR_V        [TSUOV_L]                             ;
typedef vector_O   VECTOR_O        [TSUOV_L]                             ;
typedef vector_N   VECTOR_N        [TSUOV_L]                             ;

typedef VECTOR_V   MATRIX_VxV      [TSUOV_V]                             ;
typedef VECTOR_V   MATRIX_OxV      [TSUOV_O]                             ;
typedef VECTOR_N   MATRIX_NxN      [TSUOV_N]                             ;
typedef VECTOR_N   MATRIX_OxN      [TSUOV_O]                             ;

typedef Fq         TSUOV_P3        [TSUOV_m1][TSUOV_O][TSUOV_L][TSUOV_O] ;

typedef VECTOR_O   MATRIX_kxO      [TSUOV_k];
typedef VECTOR_N   MATRIX_kxN      [TSUOV_k];
typedef VECTOR_V   MATRIX_kxV      [TSUOV_k]                             ;

typedef uint8_t TSUOV_SEED [TSUOV_SEED_LEN] ;
typedef uint8_t TSUOV_SALT [TSUOV_SALT_LEN] ;

TYPEDEF_STRUCT(TSUOV_SIGNATURE,
  TSUOV_SALT r  ;
  MATRIX_kxN s  ;
) ;

static inline void VECTOR_V_CLEAR_TAIL(VECTOR_V C){
  for(int k=0;k<TSUOV_L;k++) memset(C[k]+TSUOV_V, 0, BLOCK_ALIGNED_BYTES(TSUOV_V)-TSUOV_V) ;
}

MGF_CTX_s * PRG_init_sm3(const uint8_t seed[TSUOV_SEED_LEN]);
void PRG_yield_sm3(MGF_CTX_s * ctx, int length, uint8_t * dst);
void PRG_final_sm3(MGF_CTX_s * ctx);
MGF_CTX_s * PRG_copy_sm3(MGF_CTX_s * src);
int RejSampPRG_sm3(MGF_CTX_s * ctx, const uint64_t index,
                   const unsigned int length, const unsigned int tau, uint8_t * dst);

/* 0 on success, -1 when an internal allocation or the placeholder XOF fails. */
int TSUOV_KeyGen(const TSUOV_SEED seed_sk, const TSUOV_SEED seed_pk, TSUOV_P3 P3);
int TSUOV_Sign(const TSUOV_SEED seed_sk, const TSUOV_SEED seed_pk,
               const TSUOV_SEED seed_v, const TSUOV_SEED seed_r, const TSUOV_SEED seed_sol,
               const uint8_t message[], const size_t message_length, TSUOV_SIGNATURE sig);
/* 1 valid, 0 invalid, -1 on internal failure. */
int TSUOV_Verify(const TSUOV_SEED seed_pk, const TSUOV_P3 P3,
                 const uint8_t message[], const size_t message_length,
                 const TSUOV_SIGNATURE sig);

void store_TSUOV_P3(const TSUOV_P3 P3, uint8_t * pool, size_t * pool_bits);
/* Reject 5-bit codes equal to q (sign-31-6) and non-zero padding bits. */
int restore_TSUOV_P3(const uint8_t * pool, size_t * pool_bits, TSUOV_P3 P3);

inline static void store_TSUOV_SEED(const TSUOV_SEED seed, uint8_t * pool, size_t * pool_bits){
  size_t index = (*pool_bits >> 3) ;
  memcpy(pool+index, seed, TSUOV_SEED_LEN) ;
  *pool_bits += (TSUOV_SEED_LEN<<3) ;
}

inline static void restore_TSUOV_SEED(const uint8_t * pool, size_t * pool_bits, TSUOV_SEED seed){
  size_t index = (*pool_bits >> 3) ;
  memcpy(seed, pool+index, TSUOV_SEED_LEN) ;
  *pool_bits += (TSUOV_SEED_LEN<<3) ;
}

inline static void store_TSUOV_SALT(const TSUOV_SALT salt, uint8_t * pool, size_t * pool_bits){
  size_t index = (*pool_bits >> 3) ;
  memcpy(pool+index, salt, TSUOV_SALT_LEN) ;
  *pool_bits += (TSUOV_SALT_LEN<<3) ;
}

inline static void restore_TSUOV_SALT(const uint8_t * pool, size_t * pool_bits, TSUOV_SALT salt){
  size_t index = (*pool_bits >> 3) ;
  memcpy(salt, pool+index, TSUOV_SALT_LEN) ;
  *pool_bits += (TSUOV_SALT_LEN<<3) ;
}

inline static void store_TSUOV_SIGNATURE(const TSUOV_SIGNATURE sig, uint8_t * pool, size_t * pool_bits){
  store_TSUOV_SALT(sig->r, pool, pool_bits) ;
  for(size_t k=0; k<TSUOV_k; k++)
    for(size_t l=0; l<TSUOV_L; l++)
      for(size_t i=0; i<TSUOV_N; i++)
        store_Fq(sig->s[k][l][i], pool, pool_bits) ;
}

inline static int restore_TSUOV_SIGNATURE(const uint8_t * pool, size_t * pool_bits, TSUOV_SIGNATURE sig){
  restore_TSUOV_SALT(pool, pool_bits, sig->r) ;
  for(size_t k=0; k<TSUOV_k; k++)
    for(size_t l=0; l<TSUOV_L; l++)
      for(size_t i=0; i<TSUOV_N; i++)
        if (restore_Fq(pool, pool_bits, &sig->s[k][l][i]) != 0)
          return -1 ;
  for(size_t i = 0; i < TSUOV_k; i++)
    for(size_t j = 0; j < TSUOV_L; j++)
      memset(&(sig->s[i][j][TSUOV_N]), 0, (BLOCK_ALIGNED_BYTES(TSUOV_N) - TSUOV_N) * sizeof(uint8_t));
  return 0 ;
}
