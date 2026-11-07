#include "mayo.h"
#include "pqc_signature.h"
#include <stdint.h>
static int keypair(uint8_t *pk,uint8_t *sk){return mayo_keypair(NULL,pk,sk);}
static int sign_adapter(uint8_t *out,size_t *len,const uint8_t *m,size_t n,const uint8_t *sk){if(n>SIZE_MAX-PARAM_NAME(sig_bytes))return -1;return mayo_sign(NULL,out,len,m,n,sk);}
static int open_adapter(uint8_t *out,size_t *len,const uint8_t *sm,size_t n,const uint8_t *pk){if(n<PARAM_NAME(sig_bytes))return -1;return mayo_open(NULL,out,len,sm,n,pk);}
static const pqc_signature algorithm={PARAM_NAME(name),PARAM_NAME(cpk_bytes),PARAM_NAME(csk_bytes),PARAM_NAME(sig_bytes),keypair,sign_adapter,open_adapter};
const pqc_signature *pqc_mayo_algorithm(void){return &algorithm;}
