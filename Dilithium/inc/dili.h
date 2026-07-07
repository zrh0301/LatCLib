#ifndef DILI_H
#define DILI_H

#define DILI_N 256
#define DILI_Q 8380417
#define DILI_QBITS 23
#define DILI_D 13

#define DILI_K 8  // Max for Dilithium5
#define DILI_L 7  // Max for Dilithium5

#define DILI_BETA 78
#define DILI_OMEGA 80
#define DILI_TAU 60
#define DILI_GAMMA1 (1<<19)
#define DILI_GAMMA2 ((DILI_Q-1)/32)
#define DILI_ALPHA (2*DILI_GAMMA2)

#define DILI_POLYT1_PACKEDBYTES 320                       // t1: 10-bit x 256 / 8
#define DILI_POLYZ_PACKEDBYTES(g1) ((g1)==(1<<17) ? 576 : 640)  // z: 18-bit(576) / 20-bit(640)

#include <stdint.h>
#include "poly.h"
#include "keccak.h"

typedef struct {
    uint8_t rho[32];
    poly t1[DILI_K];
} dili_pk;

typedef struct {
    uint8_t rho[32];
    uint8_t K[32];
    uint8_t tr[64];
    poly s1[DILI_L];
    poly s2[DILI_K];
    poly t0[DILI_K];
} dili_sk;

typedef struct {
    dili_pk pk;
    dili_sk sk;
} dili_keypair_t;

typedef struct {
    int mode;             // 0 for Dilithium2, 1 for Dilithium3, 2 for Dilithium5
    int k, l;             // 矩阵维度
    int eta;              // 密钥系数范围
    int tau;              // 挑战多项式中 ±1 的个数
    int beta;             // 拒绝采样界
    int omega;            // hint 的最大权重
    int gamma1;           // 掩码 y 的范围
    int gamma2;           // 舍入参数
    int w1_packed_bytes;  // 每个 w1 多项式打包后的字节数 (D2: 192, D3/D5: 128)
    int c_tilde_bytes;    // c~ 挑战哈希字节数 = lambda/4 (ML-DSA-44/65/87: 32/48/64)
} dili_param;

void dili_keygen(dili_keypair_t* keypair, dili_param param);
int dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk, dili_param param);
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk, dili_param param);

// Packing/Unpacking functions
void dili_pack_pk(uint8_t *pk_bytes, const uint8_t *rho, const poly *t1, int k);
void dili_unpack_pk(uint8_t *rho, poly *t1, const uint8_t *pk_bytes, int k);

#endif // DILI_H