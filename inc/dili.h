//参数
#define DILI_N 256
#define DILI_Q 8380417
#define DILI_QBITS 23
#define DILI_D 13
#define DILI_BETA 78
#define DILITHIUM2

#ifdef DILITHIUM2
#define DILI_OMEGA 80
#define DILI_TAU 39
#define DILI_GAMMA1 (1<<17)
#define DILI_GAMMA2 ((DILI_Q-1)/88)
#define DILI_K 4
#define DILI_L 4
#define DILI_ALPHA (2*DILI_GAMMA2)
#endif

#include <stdint.h>
#include "poly.h"
#include "keccak.h"

typedef struct {
    uint8_t rho[32];
    poly t1[DILI_N];
} dili_pk;

typedef struct {
    uint8_t rho[32];
    uint8_t K[32];
    uint8_t tr[64];
    poly s1[DILI_N];
    poly s2[DILI_N];
    poly t0[DILI_N];
} dili_sk;

typedef struct {
    dili_pk pk;
    dili_sk sk;
} dili_keypair_t;

typedef struct {
    int mode; // 0 for Dilithium2, 1 for Dilithium3, 2 for Dilithium5
    int k, l;
    int eta;
} dili_param;

void dili_keygen(dili_keypair_t* keypair, dili_param param);
int dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk, dili_param param);
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk, dili_param param);

// Packing/Unpacking functions
void dili_pack_pk(uint8_t *pk_bytes, const uint8_t *rho, const poly *t1, int k);
void dili_unpack_pk(uint8_t *rho, poly *t1, const uint8_t *pk_bytes, int k);