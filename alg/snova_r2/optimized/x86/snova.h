#ifndef SNOVA_H
#define SNOVA_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifndef SNOVA_q
#include "snova_params.h"
#endif

#if ((SNOVA_v * SNOVA_l * SNOVA_o) & 1) != 0

#error "Not supported"
#endif

#ifndef SNOVA_m1
#define SNOVA_m1 ((SNOVA_o * SNOVA_r + SNOVA_l - 1) / SNOVA_l)
#endif

#if SNOVA_q != 16
#define SYMMETRIC
#endif

#ifndef SNOVA_r
#define SNOVA_r SNOVA_l
#endif

#ifndef FIXED_ABQ
#define FIXED_ABQ ((SNOVA_q != 16) || (SNOVA_l < 4))
#endif

#ifdef AESCTR
#define PKX_NAME _aes_
#else
#define PKX_NAME _
#endif

#if SNOVA_l == SNOVA_r
#define PARAM_JOIN_(o, p, a, b, c, d, f) _snova_##a##_##b##_##c##_##d##p##o##_##f
#define PARAM_JOIN(o, p, a, b, c, d, f) PARAM_JOIN_(o, p, a, b, c, d, f)
#define SNOVA_NAMESPACE(f) PARAM_JOIN(SNOVA_OPT, PKX_NAME, SNOVA_v, SNOVA_o, SNOVA_q, SNOVA_l, f)
#else
#define PARAM_JOIN_(o, p, a, b, c, d, r, f) _snova_##a##_##b##_##c##_##d##x##r##p##o##_##f
#define PARAM_JOIN(o, p, a, b, c, d, r, f) PARAM_JOIN_(o, p, a, b, c, d, r, f)
#define SNOVA_NAMESPACE(f) PARAM_JOIN(SNOVA_OPT, PKX_NAME, SNOVA_v, SNOVA_o, SNOVA_q, SNOVA_l, SNOVA_r, f)
#endif

#define SEED_LENGTH_PUBLIC 16
#define SEED_LENGTH_PRIVATE 32
#define SEED_LENGTH (SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE)

#define BYTES_SALT 16
#define BYTES_DIGEST 64

#define SNOVA_n (SNOVA_v + SNOVA_o)
#define SNOVA_l2 (SNOVA_l * SNOVA_l)
#define SNOVA_r2 (SNOVA_r * SNOVA_r)
#define SNOVA_lr (SNOVA_l * SNOVA_r)

#ifdef SYMMETRIC
#define NUMGF_PK (SNOVA_m1 * SNOVA_o * SNOVA_l * (SNOVA_o * SNOVA_l + 1) / 2)
#else
#define NUMGF_PK (SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2)
#endif

#define NUMGF_SIGNATURE (SNOVA_n * SNOVA_lr)

#if SNOVA_q == 11
#define Q_A 0
#define Q_B 3
#define Q_C 6
#define PACK_GF 16
#define PACK_BYTES 7

#elif SNOVA_q == 13
#define Q_A 2
#define Q_B 11
#define Q_C 3
#define PACK_GF 15
#define PACK_BYTES 7

#elif SNOVA_q == 16
#define PACK_GF 2
#define PACK_BYTES 1

#elif SNOVA_q == 17
#define Q_A 1
#define Q_B 11
#define Q_C 10
#define PACK_GF 15
#define PACK_BYTES 8

#elif SNOVA_q == 19
#define Q_A 1
#define Q_B 3
#define Q_C 15
#define PACK_GF 15
#define PACK_BYTES 8

#elif SNOVA_q == 23
#define Q_A 1
#define Q_B 11
#define Q_C 22
#define PACK_GF 7
#define PACK_BYTES 4

#elif SNOVA_q == 29
#define Q_A 3
#define Q_B 12
#define Q_C 11
#define PACK_GF 13
#define PACK_BYTES 8

#elif SNOVA_q == 31
#define Q_A 2
#define Q_B 5
#define Q_C 8
#define PACK_GF 8
#define PACK_BYTES 5

#else
#error "Parameters not supported"
#endif

#define BYTES_GF(x) ((PACK_BYTES * (x) + PACK_GF - 1) / PACK_GF)

#define BYTES_PK (BYTES_GF(NUMGF_PK) + SEED_LENGTH_PUBLIC)
#define BYTES_SIGNATURE (BYTES_GF(NUMGF_SIGNATURE) + BYTES_SALT)

#define GF16_HASH (SNOVA_o * SNOVA_l * SNOVA_r)
#define BYTES_HASH (BYTES_GF(GF16_HASH))

#ifndef SNOVA_alpha
#define SNOVA_alpha (SNOVA_l * SNOVA_r + SNOVA_r)
#endif

#ifdef SYMMETRIC
#if FIXED_ABQ
#define NUM_GEN_PUB_GF (SNOVA_m1 * (SNOVA_v * (SNOVA_v + 1) / 2 + SNOVA_v * SNOVA_o) * SNOVA_l2)
#else
#define NUM_GEN_PUB_GF                                                         \
    (SNOVA_m1 * (SNOVA_v * (SNOVA_v + 1) / 2 + SNOVA_v * SNOVA_o) * SNOVA_l2 + \
     (SNOVA_o * SNOVA_alpha) * (SNOVA_r2 + SNOVA_lr) + 2 * SNOVA_o * SNOVA_alpha * SNOVA_l)
#endif
#else
#define NUM_GEN_PUB_GF                                                                                                     \
    (SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2 + (SNOVA_o * SNOVA_alpha) * (SNOVA_r2 + SNOVA_lr) + \
     2 * SNOVA_o * SNOVA_alpha * SNOVA_l)
#endif
#define NUM_PUB_GF                                                                                                       \
    (SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2 + SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr) + \
     2 * SNOVA_o * SNOVA_alpha * SNOVA_l)

#if SNOVA_q != 16
#define NUM_GEN_PUB_BYTES (NUM_GEN_PUB_GF)
#else
#define NUM_GEN_PUB_BYTES ((NUM_GEN_PUB_GF + 1) / 2)
#endif
#define NUM_GEN_SEC_BYTES (BYTES_GF(SNOVA_v * SNOVA_lr))

#define i_prime(mi, alpha) ((alpha + mi) % SNOVA_m1)

typedef struct {
	uint16_t P11[SNOVA_m1 * SNOVA_n * SNOVA_n * SNOVA_l2];
	uint16_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];
	uint16_t F21[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2];
#ifndef SYMMETRIC
	uint16_t F12[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2];
#endif
	uint16_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
	uint16_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
	uint16_t Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
	uint16_t Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2];
	uint16_t q1[SNOVA_o * SNOVA_alpha * SNOVA_l];
	uint16_t q2[SNOVA_o * SNOVA_alpha * SNOVA_l];
	uint8_t sk_seed[SEED_LENGTH];
} expanded_SK;

typedef struct {
	uint16_t P[SNOVA_m1 * SNOVA_n * SNOVA_n * SNOVA_l2];
	uint8_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
	uint8_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
	uint8_t q1[SNOVA_o * SNOVA_alpha * SNOVA_l];
	uint8_t q2[SNOVA_o * SNOVA_alpha * SNOVA_l];
	uint8_t pk_seed[SEED_LENGTH_PUBLIC];
} expanded_PK;

int SNOVA_NAMESPACE(genkeys)(uint8_t* pk, uint8_t* sk, const uint8_t* seed);
int SNOVA_NAMESPACE(sk_expand)(expanded_SK* skx, const uint8_t* sk);
int SNOVA_NAMESPACE(sign)(const expanded_SK* skx, uint8_t* sig, const uint8_t* digest, const size_t len_digest,
                          const uint8_t *salt);
int SNOVA_NAMESPACE(pk_expand)(expanded_PK* pkx, const uint8_t* pk);
int SNOVA_NAMESPACE(verify)(const expanded_PK* pkx, const uint8_t* sig, const uint8_t* digest, const size_t len_digest);

#endif
