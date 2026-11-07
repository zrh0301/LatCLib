/* SPDX-License-Identifier: MIT */
#ifndef R3_ARM_STAGE_H
#define R3_ARM_STAGE_H
#if defined(__aarch64__)
#include <arm_neon.h>
#include "round3_arm_tables.h"
#ifndef R3_ARM_CONTRACT
#define R3_ARM_CONTRACT 1
#endif
#ifndef R3_ARM_LINEAR
#define R3_ARM_LINEAR 1
#endif
#ifndef R3_ARM_SECONDARY
#define R3_ARM_SECONDARY 1
#endif
#define R3_ARM_NB (SNOVA_lr32/16)
#define R3_ARM_RP ((SNOVA_l+1)/2)
#define R3_ARM_RR ((SNOVA_r2+15)/16*16)

/* Secret pair -> register table: no secret-indexed memory read. */
static inline uint8x16_t r3_arm_secret_table(uint8_t pair) {
    const uint8x16_t ids={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    uint8x16_t b1=vdupq_n_u8(pair),b2=round3_neon_xtime(b1);
    uint8x16_t b4=round3_neon_xtime(b2),b8=round3_neon_xtime(b4);
    uint8x16_t t1=vandq_u8(b1,vtstq_u8(ids,vdupq_n_u8(1)));
    uint8x16_t t2=vandq_u8(b2,vtstq_u8(ids,vdupq_n_u8(2)));
    uint8x16_t t4=vandq_u8(b4,vtstq_u8(ids,vdupq_n_u8(4)));
    uint8x16_t t8=vandq_u8(b8,vtstq_u8(ids,vdupq_n_u8(8)));
    return veorq_u8(veorq_u8(t1,t2),veorq_u8(t4,t8));
}

/* W^T P W: pair the output rows; consume each P-row contraction immediately.
 * P is public even in signing. W is secret in signing and public in verify.
 * Output preserves the raw-u16 ABI. Odd l has a zero high output. */
static void r3_arm_contract(const uint16_t *P,const uint8_t *W,
                           uint16_t *out,int n,int sign_layout) {
    uint8x16_t inputs[R3_ARM_NB];
    for(int mi=0;mi<SNOVA_m1;mi++) {
        uint8x16_t total[R3_ARM_RP][SNOVA_r][R3_ARM_NB];
        uint8x16_t right[R3_ARM_RP][R3_ARM_NB];
        memset(total,0,sizeof(total));
        for(int ni=0;ni<n;ni++) {
            memset(right,0,sizeof(right));
            for(int nj=0;nj<n;nj++)for(int k=0;k<SNOVA_l;k++) {
                const uint8_t *w=W+(sign_layout?k*n+nj:nj*SNOVA_l+k)*SNOVA_lr32;
                const uint16_t *p=P+((mi*n+ni)*n+nj)*SNOVA_l2;
                for(int b=0;b<R3_ARM_NB;b++)inputs[b]=vld1q_u8(w+16*b);
                for(int pair=0;pair<R3_ARM_RP;pair++) {
                    unsigned a=p[2*pair*SNOVA_l+k];
                    if(2*pair+1<SNOVA_l)a|=p[(2*pair+1)*SNOVA_l+k]<<4;
                    uint8x16_t tab=vld1q_u8(r3_arm_pairs[a]);
                    for(int b=0;b<R3_ARM_NB;b++)right[pair][b]=veorq_u8(right[pair][b],vqtbl1q_u8(tab,inputs[b]));
                }
            }
            for(int k=0;k<SNOVA_l;k++) {
                for(int b=0;b<R3_ARM_NB;b++)inputs[b]=(k&1)?vshrq_n_u8(right[k/2][b],4):vandq_u8(right[k/2][b],vdupq_n_u8(15));
                const uint8_t *w=W+(sign_layout?k*n+ni:ni*SNOVA_l+k)*SNOVA_lr32;
                for(int pair=0;pair<R3_ARM_RP;pair++)for(int i=0;i<SNOVA_r;i++) {
                    uint8_t a=w[2*pair*SNOVA_r+i];
                    if(2*pair+1<SNOVA_l)a|=w[(2*pair+1)*SNOVA_r+i]<<4;
                    uint8x16_t tab=sign_layout?r3_arm_secret_table(a):vld1q_u8(r3_arm_pairs[a]);
                    for(int b=0;b<R3_ARM_NB;b++)total[pair][i][b]=veorq_u8(total[pair][i][b],vqtbl1q_u8(tab,inputs[b]));
                }
            }
        }
        for(int a=0;a<SNOVA_l;a++)for(int i=0;i<SNOVA_r;i++)for(int b=0;b<R3_ARM_NB;b++) {
            uint8x16_t v=total[a/2][i][b];
            v=(a&1)?vshrq_n_u8(v,4):vandq_u8(v,vdupq_n_u8(15));
            uint16_t *d=out+((mi*SNOVA_l+a)*SNOVA_r+i)*SNOVA_lr32+16*b;
            vst1q_u16(d,vmovl_u8(vget_low_u8(v)));
            vst1q_u16(d+8,vmovl_u8(vget_high_u8(v)));
        }
        snova_secure_clear(total,sizeof(total));snova_secure_clear(right,sizeof(right));
    }
    snova_secure_clear(inputs,sizeof(inputs));
}

/* Two secret F folds share W loads. Pair F output coordinates in nibbles;
 * write raw u16 only at the native sign-linear boundary, with no PK/SK cache. */
static void r3_arm_linear(uint16_t *out12,uint16_t *out21,
                         const uint16_t *F12,const uint16_t *F21,const uint8_t *W) {
    uint8x16_t inputs[R3_ARM_NB];
    for(int mi=0;mi<SNOVA_m1;mi++)for(int o=0;o<SNOVA_o;o++) {
        uint8x16_t s12[R3_ARM_RP][R3_ARM_NB],s21[R3_ARM_RP][R3_ARM_NB];
        memset(s12,0,sizeof(s12));memset(s21,0,sizeof(s21));
        for(int n=0;n<SNOVA_v;n++)for(int k=0;k<SNOVA_l;k++) {
            const uint16_t *f12=F12+((mi*SNOVA_v+n)*SNOVA_o+o)*SNOVA_l2;
            const uint16_t *f21=F21+((mi*SNOVA_o+o)*SNOVA_v+n)*SNOVA_l2;
            for(int b=0;b<R3_ARM_NB;b++)inputs[b]=vld1q_u8(W+(k*SNOVA_v+n)*SNOVA_lr32+16*b);
            for(int p=0;p<R3_ARM_RP;p++) {
                uint8_t a=f12[k*SNOVA_l+2*p],c=f21[2*p*SNOVA_l+k];
                if(2*p+1<SNOVA_l){a|=f12[k*SNOVA_l+2*p+1]<<4;c|=f21[(2*p+1)*SNOVA_l+k]<<4;}
                uint8x16_t ta=r3_arm_secret_table(a),tc=r3_arm_secret_table(c);
                for(int b=0;b<R3_ARM_NB;b++) {
                    s12[p][b]=veorq_u8(s12[p][b],vqtbl1q_u8(ta,inputs[b]));
                    s21[p][b]=veorq_u8(s21[p][b],vqtbl1q_u8(tc,inputs[b]));
                }
            }
        }
        for(int i=0;i<SNOVA_l;i++) {
            uint8_t row21[SNOVA_lr32];
            for(int b=0;b<R3_ARM_NB;b++) {
                uint8x16_t a=s12[i/2][b],c=s21[i/2][b];
                a=(i&1)?vshrq_n_u8(a,4):vandq_u8(a,vdupq_n_u8(15));
                c=(i&1)?vshrq_n_u8(c,4):vandq_u8(c,vdupq_n_u8(15));
                uint16_t *d=out12+((mi*SNOVA_l+i)*SNOVA_o+o)*SNOVA_lr32+16*b;
                vst1q_u16(d,vmovl_u8(vget_low_u8(a)));vst1q_u16(d+8,vmovl_u8(vget_high_u8(a)));
                vst1q_u8(row21+16*b,c);
            }
            for(int ab=0;ab<SNOVA_l;ab++)for(int j=0;j<SNOVA_r;j++)
                out21[((mi*SNOVA_l+ab)*SNOVA_o+o)*SNOVA_lr+i*SNOVA_r+j]=row21[ab*SNOVA_r+j];
            snova_secure_clear(row21,sizeof(row21));
        }
        snova_secure_clear(s12,sizeof(s12));snova_secure_clear(s21,sizeof(s21));
    }
    snova_secure_clear(inputs,sizeof(inputs));
}

/* Fuse public q1/q2 coefficients before contracting U. U is converted once
 * per call. A/B use the original expanded arithmetic; outputs are raw u16. */
static void r3_arm_secondary(uint16_t *out,const uint16_t *U,
                             const uint16_t *A,const uint16_t *B,
                             const uint16_t *q1,const uint16_t *q2) {
    uint8_t u[SNOVA_m1*SNOVA_l2][R3_ARM_RR];
    memset(u,0,sizeof(u));
    for(int z=0;z<SNOVA_m1*SNOVA_l2;z++)for(int t=0;t<SNOVA_r2;t++)u[z][t]=(uint8_t)U[z*SNOVA_r2+t];
    memset(out,0,SNOVA_o*SNOVA_lr*sizeof(*out));
    for(int mi=0;mi<SNOVA_o;mi++)for(int al=0;al<SNOVA_alpha;al++) {
        int ix=mi*SNOVA_alpha+al,mp=i_prime(mi,al);
        uint8_t g[R3_ARM_RR];uint16_t gb[SNOVA_lr]={0};
        uint8x16_t acc[R3_ARM_RR/16];memset(acc,0,sizeof(acc));
        for(int a=0;a<SNOVA_l;a++)for(int b=0;b<SNOVA_l;b++) {
            unsigned coeff=gf16_compress(gf16_compress(q1[ix*SNOVA_l+a])*q2[ix*SNOVA_l+b]);
            uint8x16_t tab=vld1q_u8(r3_arm_pairs[coeff]);
            for(int t=0;t<R3_ARM_RR/16;t++)acc[t]=veorq_u8(acc[t],vqtbl1q_u8(tab,vld1q_u8(u[(mp*SNOVA_l+a)*SNOVA_l+b]+t*16)));
        }
        for(int t=0;t<R3_ARM_RR/16;t++)vst1q_u8(g+16*t,acc[t]);
        for(int i=0;i<SNOVA_r;i++)for(int j=0;j<SNOVA_l;j++) {
            uint16_t x=0;
            for(int k=0;k<SNOVA_r;k++)x^=g[i*SNOVA_r+k]*B[ix*SNOVA_lr+k*SNOVA_l+j];
            gb[i*SNOVA_l+j]=gf16_compress(x);
        }
        for(int i=0;i<SNOVA_r;i++)for(int j=0;j<SNOVA_l;j++) {
            uint16_t x=0;
            for(int k=0;k<SNOVA_r;k++)x^=A[ix*SNOVA_r2+i*SNOVA_r+k]*gb[k*SNOVA_l+j];
            out[mi*SNOVA_lr+i*SNOVA_l+j]^=gf16_compress(x);
        }
        snova_secure_clear(acc,sizeof(acc));snova_secure_clear(g,sizeof(g));snova_secure_clear(gb,sizeof(gb));
    }
    snova_secure_clear(u,sizeof(u));
}
#endif
#endif
