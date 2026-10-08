#include "test.h"
#include "polyvec.h"
#include "keccak.h"
#include <stdio.h>
#include <string.h>

#define ROUNDS 10000
static void report(const char *name, uint64_t cycles) {
    printf("%-28s %12llu cycles/op\n", name, (unsigned long long)(cycles / ROUNDS));
}
#define MEASURE(name, ...) do { uint64_t _start=test_ticks(); for(int i=0;i<ROUNDS;i++){ __VA_ARGS__; } uint64_t b=test_ticks(); report(name,b-_start); } while(0)

static int rng_cb(void *ctx, uint8_t *out, size_t n) { (void)ctx; memset(out, 0x5a, n); return 0; }
int main(void) {
    fp32_ctx f; ring_ctx r; ntt_plan p; poly a,b,c,v[4],m[4*4];
    uint8_t bytes[168], digest[64]; KECCAK1600_CTX h;
    fp32_init(&f, 8380417); ring_init_negacyclic(&r,&f,256); ntt_plan_init_256(&p,&r,8);
    for(int i=0;i<256;i++){a[i]=i;b[i]=255-i;}
    for(int i=0;i<4;i++){memcpy(v[i],a,sizeof(a));}
    printf("Operator performance (CPU cycles)\n");
    MEASURE("fp32 add", fp32_add(a[0],b[0],&f));
    MEASURE("fp32 sub", fp32_sub(a[0],b[0],&f));
    MEASURE("fp32 mul", fp32_mul(a[0],b[0],&f));
    MEASURE("fp32 montgomery reduce", fp32_mont_reduce((int64_t)a[0]*b[0],&f));
    MEASURE("fp32 inverse", fp32_inv(&c[0],a[1],&f));
    MEASURE("poly add (256)", POLY_ADD(c,a,b,256));
    MEASURE("poly sub (256)", POLY_SUB(c,a,b,256));
    MEASURE("poly reduce (256)", poly_reduce_mod(c,256,f.q));
    MEASURE("poly NTT (256)", memcpy(c,a,sizeof(a)), poly_ntt(c,&p));
    MEASURE("poly inverse NTT (256)", memcpy(c,a,sizeof(a)), poly_invntt(c,&p));
    MEASURE("poly NTT pointwise", poly_mul_ntt(c,a,b,&p));
    MEASURE("polyvec add (4x256)", polyvec_add(v,v,v,4));
    MEASURE("polyvec NTT (4)", polyvec_ntt(v,4,&p));
    MEASURE("polyvec inner product (4)", polyvec_pointwise_acc(c,(const poly *)v,(const poly *)v,4,&p));
    MEASURE("polyvec matrix product (4x4)", polyvec_matrix_pointwise(v,m,4,4,(const poly *)v,&p));
    MEASURE("CBD sampling", sample_CBD(c,256,bytes,sizeof(bytes),&p));
    MEASURE("random bytes (168)", pqc_randombytes(bytes,sizeof(bytes),rng_cb,NULL));
    MEASURE("SHA3-256 (168 bytes)", sha3_init(&h,0x06,256), sha3_update(&h,bytes,sizeof(bytes)), sha3_final(digest,&h));
    return 0;
}
