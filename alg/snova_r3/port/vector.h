#ifndef R3_PORT_VECTOR_H
#define R3_PORT_VECTOR_H
#include <string.h>
#ifndef R3_METHOD
#define R3_METHOD 0
#endif
#if defined(__AVX2__)
#include <immintrin.h>
typedef __m256i pv;
#define PV_BYTES 32
#define px _mm256_xor_si256
#define pa _mm256_and_si256
#define po _mm256_or_si256
#define ps _mm256_set1_epi8
#define pz() _mm256_setzero_si256()
#define pe _mm256_cmpeq_epi8
static inline pv pl(const void *p) { return _mm256_loadu_si256((const __m256i*)p); }
static inline void pst(void *p,pv v) { _mm256_storeu_si256((__m256i*)p,v); }
static inline pv shr4(pv v) { return pa(_mm256_srli_epi16(v,4),ps(15)); }
static inline pv shl4(pv v) { return pa(_mm256_slli_epi16(v,4),ps(240)); }
static inline pv shuf(pv t,pv i) { return _mm256_shuffle_epi8(t,i); }
static inline pv psub(pv a,pv b) { return _mm256_sub_epi8(a,b); }
static inline pv pdouble(pv b) { return _mm256_add_epi8(b,b); }
#elif defined(__aarch64__)
#include <arm_neon.h>
typedef uint8x16_t pv;
#define PV_BYTES 16
#define px veorq_u8
#define pa vandq_u8
#define po vorrq_u8
#define ps vdupq_n_u8
#define pz() vdupq_n_u8(0)
#define pe vceqq_u8
#define psub vsubq_u8
static inline pv pl(const void *p) { return vld1q_u8((const uint8_t*)p); }
static inline void pst(void *p,pv v) { vst1q_u8((uint8_t*)p,v); }
static inline pv shr4(pv v) { return vshrq_n_u8(v,4); }
static inline pv shl4(pv v) { return vshlq_n_u8(v,4); }
static inline pv shuf(pv t,pv i) { return vqtbl1q_u8(t,i); }
static inline pv pdouble(pv b) { return vaddq_u8(b,b); }
#else
#error "The full port requires AVX2 or AArch64 NEON"
#endif
#define PAD(x) (((x)+PV_BYTES-1)/PV_BYTES*PV_BYTES)

static inline uint8_t pmul(uint8_t a,uint8_t b) {
    uint16_t t = (uint16_t)(gf16_expand(a) * b);
    unsigned p=(t&15)^((t>>3)&30)^((t>>6)&60)^((t>>9)&120);
    return (uint8_t)((p^(p>>4)^((p>>4)<<1))&15);
}
static inline pv pxtime(pv b) {
    return px(pa(pdouble(b),ps(15)),pa(pe(pa(b,ps(8)),ps(8)),ps(3)));
}
static inline pv ppacked_xtime(pv b) {
#if defined(__AVX2__)
    pv carry=pa(_mm256_srli_epi16(b,3),ps(17));
#else
    pv carry=pa(vshrq_n_u8(b,3),ps(17));
#endif
    return px(pa(pdouble(b),ps(238)),px(carry,pdouble(carry)));
}
static inline pv pbit(pv a,pv b) {
    pv b2=pxtime(b),b4=pxtime(b2),b8=pxtime(b4);
    pv t1=pa(b,pe(pa(a,ps(1)),ps(1)));
    pv t2=pa(b2,pe(pa(a,ps(2)),ps(2)));
    pv t4=pa(b4,pe(pa(a,ps(4)),ps(4)));
    pv t8=pa(b8,pe(pa(a,ps(8)),ps(8)));
    return px(px(t1,t2),px(t4,t8));
}
static inline pv pmulv(pv a,pv b) {
#if defined(__GFNI__) && R3_METHOD == 1
    pv c=_mm256_gf2p8mul_epi8(a,b);
    return _mm256_gf2p8affine_epi64_epi8(c,_mm256_set1_epi64x(0x113264c800000000ULL),0);
#elif R3_METHOD == 3
    static const uint8_t logs[32]={0,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12,0,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};
    static const uint8_t exps[32]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9,1,1,2,4,8,3,6,12,11,5,10,7,14,15,13,9,1};
    pv x=shuf(pl(logs),a),y=shuf(pl(logs),b),s=psub(x,psub(pz(),y));
#if defined(__AVX2__)
    pv carry=_mm256_cmpgt_epi8(s,ps(14));
#else
    pv carry=vcgtq_u8(s,ps(14));
#endif
    s=psub(s,pa(carry,ps(15)));
    pv nonzero=px(po(pe(a,pz()),pe(b,pz())),ps(255));
    return pa(shuf(pl(exps),s),nonzero);
#else
    return pbit(a,b);
#endif
}

#if defined(__GFNI__)
static inline uint64_t paff(uint8_t c) {
    return ((0ULL-(uint64_t)(c&1))&0x0102040810204080ULL)^
           ((0ULL-(uint64_t)((c>>1)&1))&0x0809020480902040ULL)^
           ((0ULL-(uint64_t)((c>>2)&1))&0x040c090240c09020ULL)^
           ((0ULL-(uint64_t)((c>>3)&1))&0x02060c092060c090ULL);
}
#endif
/* B and C stay packed throughout this contraction. A supplies scalar coefficients. */
static inline pv ptab(uint8_t c);
#include "packed_a.h"
static void pmm_packed(uint8_t *c,int cs,const uint8_t *a,int as,
                       const uint8_t *b,int bs,int rows,int inner,int packed_cols,int public_coeff) {
#if R3_METHOD == 8
    (void)public_coeff;
    uint8_t ap[((rows+1)/2)*PAD(inner)];
    ppack_a(ap,PAD(inner),a,as,rows,inner);
    pmm_ab_split(c,cs,ap,PAD(inner),b,bs,rows,inner,packed_cols);
    SNOVA_SECURE_CLEAR(ap);
#elif R3_METHOD == 7
    (void)public_coeff;
    for(int i=0;i<rows;i+=4)for(int j=0;j<PAD(packed_cols);j+=PV_BYTES) {
        pv acc[4]={pz(),pz(),pz(),pz()};
        for(int k=0;k<inner;k++) {
            pv b1=pl(b+k*bs+j),b2=ppacked_xtime(b1),b4=ppacked_xtime(b2),b8=ppacked_xtime(b4);
            for(int r=0;r<4 && i+r<rows;r++) {
                unsigned s=a[(i+r)*as+k];
                pv t1=pa(b1,ps((uint8_t)-(int)(s&1))),t2=pa(b2,ps((uint8_t)-(int)((s>>1)&1)));
                pv t4=pa(b4,ps((uint8_t)-(int)((s>>2)&1))),t8=pa(b8,ps((uint8_t)-(int)((s>>3)&1)));
                acc[r]=px(acc[r],px(px(t1,t2),px(t4,t8)));
            }
        }
        for(int r=0;r<4 && i+r<rows;r++)pst(c+(i+r)*cs+j,acc[r]);
    }
#else
#if defined(__GFNI__)
    uint64_t public_matrices[16];
    for(int k=0;k<16;k++)public_matrices[k]=paff(k);
#else
    (void)public_coeff;
#endif
    for(int i=0;i<rows;i++)for(int j=0;j<PAD(packed_cols);j+=PV_BYTES) {
        pv acc=pz();
        for(int k=0;k<inner;k++) {
#if defined(__GFNI__)
            pv mat=_mm256_set1_epi64x(public_coeff?public_matrices[a[i*as+k]]:paff(a[i*as+k]));
            acc=px(acc,_mm256_gf2p8affine_epi64_epi8(pl(b+k*bs+j),mat,0));
#else
            pv tab=ptab(a[i*as+k]),v=pl(b+k*bs+j);
            acc=px(acc,po(shuf(tab,pa(v,ps(15))),shl4(shuf(tab,shr4(v)))));
#endif
        }
        pst(c+i*cs+j,acc);
    }
#endif
}
/* Tables are built in registers, including for secret coefficients. */
static inline pv ptab(uint8_t c) {
    static const uint8_t ids[32]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                                0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
#if defined(__GFNI__)
    pv v=_mm256_gf2p8mul_epi8(pl(ids),ps(c));
    return _mm256_gf2p8affine_epi64_epi8(v,_mm256_set1_epi64x(0x113264c800000000ULL),0);
#else
    return pbit(pl(ids),ps(c));
#endif
}
static inline pv ppairtab(uint8_t a,uint8_t b) { return po(ptab(a),shl4(ptab(b))); }
static inline pv pconst(pv b,uint8_t a) {
#if R3_METHOD == 0
    return shuf(ptab(a),b);
#else
    return pmulv(ps(a),b);
#endif
}
/* Flat row-major matrices; strides include zero padding. C is overwritten.
 * MTK2 computes two output rows in each byte until the final row store.
 * The GFNI/bit-plane alternatives consume the same padded raw matrices. */
static void pmm_base(uint8_t *c,int cs,const uint8_t *a,int as,
                const uint8_t *b,int bs,int rows,int inner,int cols) {
#if R3_METHOD == 5 || R3_METHOD == 7 || R3_METHOD == 8
    int packed_cols=(cols+1)/2,bs2=PAD(packed_cols);
    uint8_t bp[inner*bs2],cp[rows*bs2];
    memset(bp,0,sizeof(bp));
    for(int k=0;k<inner;k++)for(int j=0;j<cols;j++)bp[k*bs2+j/2]|=b[k*bs+j]<<((j&1)*4);
    pmm_packed(cp,bs2,a,as,bp,bs2,rows,inner,packed_cols,0);
    memset(c,0,(size_t)rows*cs);
    for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)c[i*cs+j]=(cp[i*bs2+j/2]>>((j&1)*4))&15;
    SNOVA_SECURE_CLEAR(bp);SNOVA_SECURE_CLEAR(cp);
#else
    int stride=PAD(cols);
    for(int i=0;i<rows;i+=2) {
        for(int j=0;j<stride;j+=PV_BYTES) {
            pv s0=pz(),s1=pz();
            for(int k=0;k<inner;k++) {
                pv v=pl(b+k*bs+j);
#if R3_METHOD == 0
                pv t=ppairtab(a[i*as+k],i+1<rows?a[(i+1)*as+k]:0);
                s0=px(s0,shuf(t,v));
#else
#if defined(__GFNI__) && R3_METHOD == 1
                s0=px(s0,_mm256_gf2p8mul_epi8(ps(a[i*as+k]),v));
                if(i+1<rows)s1=px(s1,_mm256_gf2p8mul_epi8(ps(a[(i+1)*as+k]),v));
#else
                s0=px(s0,pmulv(ps(a[i*as+k]),v));
                if(i+1<rows)s1=px(s1,pmulv(ps(a[(i+1)*as+k]),v));
#endif
#endif
            }
#if R3_METHOD == 0
            s1=shr4(s0);s0=pa(s0,ps(15));
#elif defined(__GFNI__) && R3_METHOD == 1
            pv reduce=_mm256_set1_epi64x(0x113264c800000000ULL);
            s0=_mm256_gf2p8affine_epi64_epi8(s0,reduce,0);
            s1=_mm256_gf2p8affine_epi64_epi8(s1,reduce,0);
#endif
            pst(c+i*cs+j,s0);
            if(i+1<rows)pst(c+(i+1)*cs+j,s1);
        }
    }
#endif
}
#include "tiled.h"
#include "expanded.h"
#include "row_tile.h"
static inline void pmm(uint8_t *c,int cs,const uint8_t *a,int as,
                       const uint8_t *b,int bs,int rows,int inner,int cols) {
#if defined(__GFNI__) && R3_METHOD == 1 && defined(R3_WIDE_TILE) && R3_WIDE_TILE
    if(cols>PV_BYTES){pmm_tile(c,cs,a,as,b,bs,rows,inner,cols);return;}
#endif
#if defined(__GFNI__) && R3_METHOD == 1 && (R3_ROWS != 2 || R3_CHAINS != 1)
    pmm_row_tile(c,cs,a,as,b,bs,rows,inner,cols);
#elif R3_METHOD == 4
    pmm_expanded(c,cs,a,as,b,bs,rows,inner,cols);
#elif R3_METHOD == 6
    pmm_sameb4(c,cs,a,as,b,bs,rows,inner,cols);
#elif R3_TILE > 1 && R3_METHOD != 5 && R3_METHOD != 7 && R3_METHOD != 8
    pmm_tile(c,cs,a,as,b,bs,rows,inner,cols);
#else
    pmm_base(c,cs,a,as,b,bs,rows,inner,cols);
#endif
}

/* Only call this with public coefficients (verification). The table contains
 * the 256 fixed field maps, not key-dependent precomputation. */
static void pmm_public(uint8_t *c,int cs,const uint8_t *a,int as,
                       const uint8_t *b,int bs,int rows,int inner,int cols) {
#if R3_METHOD == 0
    static pv tables[256];static int ready;
    if(!ready){for(int k=0;k<256;k++)tables[k]=ppairtab(k&15,k>>4);ready=1;}
    for(int i=0;i<rows;i+=2)for(int j=0;j<PAD(cols);j+=PV_BYTES) {
        pv acc=pz();
        for(int k=0;k<inner;k++) {
            unsigned pair=a[i*as+k]|((i+1<rows?a[(i+1)*as+k]:0)<<4);
            acc=px(acc,shuf(tables[pair],pl(b+k*bs+j)));
        }
        pst(c+i*cs+j,pa(acc,ps(15)));
        if(i+1<rows)pst(c+(i+1)*cs+j,shr4(acc));
    }
#else
    pmm(c,cs,a,as,b,bs,rows,inner,cols);
#endif
}
static void ptrans(uint8_t *out,int os,const uint8_t *in,int is,int rows,int cols) {
    memset(out,0,(size_t)cols*os);
    for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)out[j*os+i]=in[i*is+j];
}
static void pblock16(uint8_t *out,int stride,const uint16_t *in,int rows,int cols) {
    memset(out,0,(size_t)rows*SNOVA_l*stride);
    for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        for(int a=0;a<SNOVA_l;a++)for(int b=0;b<SNOVA_l;b++)
            out[(i*SNOVA_l+a)*stride+j*SNOVA_l+b]=(uint8_t)in[(i*cols+j)*SNOVA_l2+a*SNOVA_l+b];
}
static void pblock8(uint8_t *out,int stride,const uint8_t *in,int rows,int cols) {
    memset(out,0,(size_t)rows*SNOVA_l*stride);
    for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        for(int a=0;a<SNOVA_l;a++)for(int b=0;b<SNOVA_l;b++)
            out[(i*SNOVA_l+a)*stride+j*SNOVA_l+b]=in[(i*cols+j)*SNOVA_l2+a*SNOVA_l+b];
}
static void punblock16(uint16_t *out,const uint8_t *in,int stride,int rows,int cols,const uint8_t *add) {
    for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        for(int a=0;a<SNOVA_l;a++)for(int b=0;b<SNOVA_l;b++) {
            int ix=(i*cols+j)*SNOVA_l2+a*SNOVA_l+b;
            out[ix]=in[(i*SNOVA_l+a)*stride+j*SNOVA_l+b]^(add?add[ix]:0);
        }
}
#endif
