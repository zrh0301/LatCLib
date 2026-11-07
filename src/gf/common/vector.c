#include "gf16.h"
#include "gf16_vector.h"
#include "gf256.h"
#include "gf2.h"
#include "gf_backend.h"
uint8_t gf256_mul(uint8_t a,uint8_t b) {
 uint8_t p=0; for(unsigned i=0;i<8;i++){p^=a & (uint8_t)(0u-(b&1u)); a=(uint8_t)((a<<1)^(0x1bu & (0u-(a>>7))));b>>=1;}return p;
}
uint8_t gf256_inv(uint8_t a) { uint8_t p=1;for(int i=7;i>=0;i--){p=gf256_mul(p,p);if((254u>>i)&1u)p=gf256_mul(p,a);}return p; }
void gf2_vector_xor(uint64_t *o,const uint64_t *a,const uint64_t *b,size_t n){for(size_t i=0;i<n;i++)o[i]=a[i]^b[i];}
uint8_t gf2_vector_dot(const uint64_t *a,const uint64_t *b,size_t n){uint64_t p=0;for(size_t i=0;i<n;i++)p^=a[i]&b[i];p^=p>>32;p^=p>>16;p^=p>>8;p^=p>>4;p^=p>>2;p^=p>>1;return p&1;}
