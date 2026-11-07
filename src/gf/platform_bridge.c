#include "gfext.h"
#include <stdatomic.h>
#include "../../inc/gf16.h"
#include "../../inc/gf256.h"
/* All public wrappers initialize imported read-only tables before use. */
static atomic_int initialized=0;
static void init(void){int expected=0;if(atomic_compare_exchange_strong_explicit(&initialized,&expected,1,memory_order_acq_rel,memory_order_acquire)){pqc_kernel_gf16_tables_init();pqc_kernel_gf256_tower_init();atomic_store_explicit(&initialized,2,memory_order_release);}else{while(atomic_load_explicit(&initialized,memory_order_acquire)!=2){}}}
const char *pqc_gf_backend_name(void){
#ifdef GF16_ENABLE_GFNI
return "x86-avx2-gfni";
#else
return GF16_BACKEND_NAME;
#endif
}
void gf16_vector_mul(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t n){
 init(); size_t i=0;
 for(;i+GF16_VEC_BYTES<=n;i+=GF16_VEC_BYTES)
  pqc_kernel_gf16_vec_store(o+i,pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec_load(a+i),pqc_kernel_gf16_vec_load(b+i)));
 for(;i<n;i++)o[i]=gf16_mul(a[i],b[i]);
}
/* Use variable multiplication for secret coefficients: no secret-indexed table. */
void gf16_vector_madd(uint8_t *o,const uint8_t *a,uint8_t k,size_t n){init();size_t i=0;for(;i+GF16_VEC_BYTES<=n;i+=GF16_VEC_BYTES){pqc_kernel_gf16_vec x=pqc_kernel_gf16_vec_load(a+i),y=pqc_kernel_gf16_vec_set1(k&15);pqc_kernel_gf16_vec_store(o+i,pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(o+i),pqc_kernel_gf16_vec_mul_bitplane(x,y)));}for(;i<n;i++){uint8_t x=a[i]&15,p=0,b=k&15;for(int j=0;j<4;j++){p^=x&(uint8_t)(0u-(b&1));x=(uint8_t)(((x<<1)&15)^(3u&(0u-(x>>3))));b>>=1;}o[i]^=p;}}
void gf16_packed_mul_const(uint8_t *o,const uint8_t *a,uint8_t k,size_t n){size_t i=0;init();for(;i+GF16_VEC_BYTES<=n;i+=GF16_VEC_BYTES){pqc_kernel_gf16_vec v=pqc_kernel_gf16_vec_load(a+i),b=pqc_kernel_gf16_vec_set1(k&15);pqc_kernel_gf16_vec lo=pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec_mask_u8(v,15),b),hi=pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec_srli4(v),b);pqc_kernel_gf16_vec_store(o+i,pqc_kernel_gf16_vec_or(lo,pqc_kernel_gf16_vec_slli4(hi)));}for(;i<n;i++){uint8_t v=a[i],lo=v&15,hi=v>>4,p=0,q=0,b=k&15;for(int j=0;j<4;j++){p^=lo&(uint8_t)(0u-(b&1));q^=hi&(uint8_t)(0u-(b&1));lo=(uint8_t)(((lo<<1)&15)^(3u&(0u-(lo>>3))));hi=(uint8_t)(((hi<<1)&15)^(3u&(0u-(hi>>3))));b>>=1;}o[i]=p|(q<<4);}}
void gf256_vector_mul(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t n){init();
#ifdef GF16_ENABLE_GFNI
pqc_kernel_gf256_mul_gfni_array(o,a,b,n);
#else
pqc_kernel_gf256_mul_array(o,a,b,n);
#endif
}
void gf256_vector_madd(uint8_t *o,const uint8_t *a,uint8_t k,size_t n){init();for(size_t i=0;i<n;i++)o[i]^=gf256_mul(a[i],k);}

void gf256_tower_vector_mul(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t n){init();pqc_kernel_gf256_tower_mul_array(o,a,b,n);}
