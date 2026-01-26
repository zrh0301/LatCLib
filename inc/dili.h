//参数
#define DILI_N 256
#define DILI_Q 8380417
#define DILI_QBITS 23
#define DILI_D 14
#define DILI_BETA 392
#define DILI_OMEGA 80
#include <stdint.h>

typedef struct {
    uint8_t seed_A[32];
    uint8_t t1[DILI_N];
} dili_pk;

typedef struct {
    uint8_t seed_A[32];
    uint8_t s1[DILI_N];
    uint8_t s2[DILI_N];
    uint8_t t0[DILI_N];
} dili_sk;

typedef struct {
    dili_pk pk;
    dili_sk sk;
} dili_keypair_t;

void dili_keygen(dili_keypair_t* keypair);
void dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk);
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk);