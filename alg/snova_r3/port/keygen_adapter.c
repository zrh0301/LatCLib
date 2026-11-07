/* A separate namespace preserves the native expanded-SK/PK ABI unchanged. */
#if defined(__GFNI__)
#ifndef R3_ROWS
#define R3_ROWS 4
#endif
#ifndef R3_CHAINS
#define R3_CHAINS 2
#endif
#ifndef R3_KEY_ORDER
#if SNOVA_l == 5 && SNOVA_o <= 5
#define R3_KEY_ORDER 1
#else
#define R3_KEY_ORDER 0
#endif
#endif
#endif
#undef SNOVA_OPT
#define SNOVA_OPT 11
#undef R3_SIGN
#define R3_SIGN 0
#undef R3_VERIFY
#define R3_VERIFY 0
#include "../ported/snova_opt_16.c"
int r3_port_genkeys(uint8_t *pk,uint8_t *sk,const uint8_t *seed) {
    return SNOVA_NAMESPACE(genkeys)(pk,sk,seed);
}
int r3_port_sk_expand(void *skx,const uint8_t *sk) {
    return SNOVA_NAMESPACE(sk_expand)((expanded_SK*)skx,sk);
}
