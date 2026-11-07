#include "api.h"
#include "pqc_signature.h"
#include <stdint.h>
static int sign_adapter(uint8_t *out,size_t *len,const uint8_t *m,size_t n,const uint8_t *sk){
 if(!out||!len||!m||!sk||n>SIZE_MAX-(size_t)CRYPTO_BYTES)return -1;
 unsigned long long l=0;int rc=crypto_sign(out,&l,m,(unsigned long long)n,sk);if(l>SIZE_MAX)return -1;*len=(size_t)l;return rc;
}
static int open_adapter(uint8_t *out,size_t *len,const uint8_t *sm,size_t n,const uint8_t *pk){
 if(!out||!len||!sm||!pk||n<(size_t)CRYPTO_BYTES)return -1;
 unsigned long long l=0;int rc=crypto_sign_open(out,&l,sm,(unsigned long long)n,pk);if(l>SIZE_MAX)return -1;*len=(size_t)l;return rc;
}
static const pqc_signature algorithm={"mqom:" CRYPTO_ALGNAME,CRYPTO_PUBLICKEYBYTES,CRYPTO_SECRETKEYBYTES,CRYPTO_BYTES,crypto_sign_keypair,sign_adapter,open_adapter};
const pqc_signature *pqc_mqom_algorithm(void){return &algorithm;}
