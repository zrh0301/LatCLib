/* SPDX-License-Identifier: MIT
 * Direct quadratic contraction: MTK2 for public P coefficients, then
 * packed dual-output bit-plane products with a shared variable right operand.
 * Only field-constant tables persist. All signature/key intermediates are local.
 */
#ifndef ROUND3_NEON_H
#define ROUND3_NEON_H
#if defined(__aarch64__)
#include <arm_neon.h>
static inline uint8x16_t round3_neon_xtime(uint8x16_t a) {
    return veorq_u8(vandq_u8(vshlq_n_u8(a,1),vdupq_n_u8(0xee)),
        vmulq_u8(vandq_u8(vshrq_n_u8(a,3),vdupq_n_u8(0x11)),vdupq_n_u8(3)));
}
static inline uint8x16_t round3_neon_select(uint8x16_t b1,uint8x16_t b2,
                                          uint8x16_t b4,uint8x16_t b8,uint8_t a) {
    uint8_t m[4];
    for(int k=0;k<4;k++)m[k]=(uint8_t)(((0u-((a>>k)&1u))&15u)|((0u-((a>>(k+4))&1u))&240u));
    return veorq_u8(veorq_u8(vandq_u8(b1,vdupq_n_u8(m[0])),vandq_u8(b2,vdupq_n_u8(m[1]))),
                    veorq_u8(vandq_u8(b4,vdupq_n_u8(m[2])),vandq_u8(b8,vdupq_n_u8(m[3]))));
}
static inline void round3_neon_contract(const uint16_t *P,const uint8_t *W,
                                        uint16_t *out,int n,int sign_layout) {
    enum {NB=SNOVA_lr32/16, RP=(SNOVA_l+1)/2};
    /* Public-only coefficient tables are rebuilt per call, included in timing. */
    uint8_t tabs[256][16];
    for(int p=0;p<256;p++)for(int x=0;x<16;x++) {
        uint8_t lo=0,hi=0,a=p&15,b=p>>4;
        for(int k=0;k<4;k++) {
            unsigned mask=0u-((x>>k)&1u);lo^=a&mask;hi^=b&mask;
            a=(uint8_t)(((a<<1)&15)^((0u-(a>>3))&3));
            b=(uint8_t)(((b<<1)&15)^((0u-(b>>3))&3));
        }
        tabs[p][x]=(uint8_t)(lo|(hi<<4));
    }
    for(int mi=0;mi<SNOVA_m1;mi++) {
        uint8x16_t total[RP][SNOVA_r][NB];
        uint8x16_t right[RP][NB];
        memset(total,0,sizeof(total));
        for(int ni=0;ni<n;ni++) {
            memset(right,0,sizeof(right));
            for(int nj=0;nj<n;nj++)for(int k=0;k<SNOVA_l;k++) {
                const uint8_t *w=W+(sign_layout?(k*n+nj):(nj*SNOVA_l+k))*SNOVA_lr32;
                const uint16_t *p=P+((mi*n+ni)*n+nj)*SNOVA_l2;
                for(int pair=0;pair<RP;pair++) {
                    unsigned lo=p[(2*pair)*SNOVA_l+k];
                    unsigned hi=2*pair+1<SNOVA_l?p[(2*pair+1)*SNOVA_l+k]:0;
                    uint8x16_t table=vld1q_u8(tabs[lo|(hi<<4)]);
                    for(int b=0;b<NB;b++)right[pair][b]=veorq_u8(right[pair][b],vqtbl1q_u8(table,vld1q_u8(w+16*b)));
                }
            }
            for(int k=0;k<SNOVA_l;k++)for(int b=0;b<NB;b++) {
                uint8x16_t raw=right[k/2][b];
                raw=(k&1)?vshrq_n_u8(raw,4):vandq_u8(raw,vdupq_n_u8(15));
                uint8x16_t b1=vorrq_u8(raw,vshlq_n_u8(raw,4));
                uint8x16_t b2=round3_neon_xtime(b1),b4=round3_neon_xtime(b2),b8=round3_neon_xtime(b4);
                const uint8_t *w=W+(sign_layout?(k*n+ni):(ni*SNOVA_l+k))*SNOVA_lr32;
                for(int pair=0;pair<RP;pair++)for(int i=0;i<SNOVA_r;i++) {
                    uint8_t a=w[2*pair*SNOVA_r+i];
                    if(2*pair+1<SNOVA_l)a|=w[(2*pair+1)*SNOVA_r+i]<<4;
                    total[pair][i][b]=veorq_u8(total[pair][i][b],round3_neon_select(b1,b2,b4,b8,a));
                }
            }
        }
        for(int a=0;a<SNOVA_l;a++)for(int i=0;i<SNOVA_r;i++)for(int b=0;b<NB;b++) {
            uint8x16_t v=total[a/2][i][b];
            v=(a&1)?vshrq_n_u8(v,4):vandq_u8(v,vdupq_n_u8(15));
            uint16_t *d=out+((mi*SNOVA_l+a)*SNOVA_r+i)*SNOVA_lr32+16*b;
            vst1q_u16(d,vmovl_u8(vget_low_u8(v)));vst1q_u16(d+8,vmovl_u8(vget_high_u8(v)));
        }
        snova_secure_clear(total,sizeof(total));
        snova_secure_clear(right,sizeof(right));
    }
}
#endif
#endif
