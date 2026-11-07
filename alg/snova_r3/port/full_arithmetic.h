#include "vector.h"
#if !FIXED_ABQ
#error "The full port currently targets the fixed-ABQ Round3 submission parameters"
#endif
#ifndef R3_KEY
#define R3_KEY 1
#endif
#ifndef R3_SIGN
#define R3_SIGN 1
#endif
#ifndef R3_VERIFY
#define R3_VERIFY 1
#endif
#ifndef R3_FVV
#define R3_FVV 1
#endif
#ifndef R3_KEY_ORDER
#define R3_KEY_ORDER 0
#endif
#define PL SNOVA_l
#define PR SNOVA_r
#define PV SNOVA_v
#define PO SNOVA_o
#define PN SNOVA_n
#define PM SNOVA_m1
#define PAL SNOVA_alpha
#define VL (PV*PL)
#define OL (PO*PL)
#define NL (PN*PL)
#define EQ (PO*PL*PR)
#define GS PAD(EQ+1)
#ifndef R3_PUBLIC_E
#define R3_PUBLIC_E 0
#endif

#if R3_PUBLIC_E
/* This is a linear map of FIXED PUBLIC ABQ, never a pk/sk/message cache.
 * Its one-time construction is analogous to upstream's fixed ABQ setup. */
static uint8_t pe_map[PM*PL*PL*PR*PR][PAD(EQ)];
static int pe_ready;
static void pe_init(void) {
    if(pe_ready)return;
    for(int m=0;m<PO;m++)for(int al=0;al<PAL;al++) {
        int ix=m*PAL+al,mp=i_prime(m,al);
        for(int c=0;c<PL;c++)for(int d=0;d<PL;d++) {
            uint8_t q=pmul(gf16_compress(fixedq1[ix*PL+c]),gf16_compress(fixedq2[ix*PL+d]));
            for(int k=0;k<PR;k++)for(int t=0;t<PR;t++)
                for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)
                    pe_map[((mp*PL+c)*PL+d)*PR*PR+k*PR+t][m*PR*PL+i*PL+j]^=
                        pmul(q,pmul(gf16_compress(fixedAm[ix*PR*PR+i*PR+k]),gf16_compress(fixedBm[ix*PR*PL+t*PL+j])));
        }
    }
    pe_ready=1;
}
#endif

/* All layouts below are call-local. Only upstream's fixed public ABQ is static. */
static void pt12(uint8_t t[VL][PAD(OL)],uint8_t tl[OL][PAD(VL)],const uint8_t *seed) {
    uint8_t raw[PV*PO*PL*PL];
    expand_T12(raw,seed);
    pblock8(&t[0][0],PAD(OL),raw,PV,PO);
    memset(tl,0,OL*PAD(VL));
    for(int v=0;v<PV;v++)for(int o=0;o<PO;o++)
        for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)
            tl[o*PL+i][v*PL+j]=raw[(v*PO+o)*PL*PL+i*PL+j];
    SNOVA_SECURE_CLEAR(raw);
}

int SNOVA_NAMESPACE(genkeys)(uint8_t *pk,uint8_t *sk,const uint8_t *seed) {
#if !R3_KEY
    return SNOVA_NAMESPACE(original_genkeys)(pk,sk,seed);
#else
    SNOVA_INIT
    uint8_t t[VL][PAD(OL)],tl[OL][PAD(VL)];
    uint8_t pub[NUM_PUB_GF],p[VL][PAD(VL)],f[VL][PAD(OL)];
    uint8_t p21[OL][PAD(VL)],p22[OL][PAD(OL)],tmp[OL][PAD(OL)];
    uint8_t raw[NUMGF_PK];
    pt12(t,tl,seed+SEED_LENGTH_PUBLIC);
    expand_public(pub,seed);
    for(int m=0;m<PM;m++) {
        pblock8(&p[0][0],PAD(VL),pub+m*PV*PV*PL*PL,PV,PV);
#if R3_KEY_ORDER == 1
        pmm(&p21[0][0],PAD(VL),&tl[0][0],PAD(VL),&p[0][0],PAD(VL),OL,VL,VL);
        const uint8_t *raw21=pub+PM*PV*PN*PL*PL+m*PO*PV*PL*PL;
        for(int o=0;o<PO;o++)for(int v=0;v<PV;v++)
            for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)
                p21[o*PL+i][v*PL+j]^=raw21[(o*PV+v)*PL*PL+i*PL+j];
        pblock8(&f[0][0],PAD(OL),pub+PM*PV*PV*PL*PL+m*PV*PO*PL*PL,PV,PO);
        pmm(&p22[0][0],PAD(OL),&p21[0][0],PAD(VL),&t[0][0],PAD(OL),OL,VL,OL);
        pmm(&tmp[0][0],PAD(OL),&tl[0][0],PAD(VL),&f[0][0],PAD(OL),OL,VL,OL);
#else
        pmm(&f[0][0],PAD(OL),&p[0][0],PAD(VL),&t[0][0],PAD(OL),VL,VL,OL);
        const uint8_t *p12=pub+PM*PV*PV*PL*PL+m*PV*PO*PL*PL;
        for(int v=0;v<PV;v++)for(int o=0;o<PO;o++)
            for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)
                f[v*PL+i][o*PL+j]^=p12[(v*PO+o)*PL*PL+i*PL+j];
        pblock8(&p21[0][0],PAD(VL),pub+PM*PV*PN*PL*PL+m*PO*PV*PL*PL,PO,PV);
        pmm(&p22[0][0],PAD(OL),&tl[0][0],PAD(VL),&f[0][0],PAD(OL),OL,VL,OL);
        pmm(&tmp[0][0],PAD(OL),&p21[0][0],PAD(VL),&t[0][0],PAD(OL),OL,VL,OL);
#endif
        for(int a=0;a<PO;a++)for(int b=0;b<PO;b++)
            for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)
                raw[((m*PO+a)*PO+b)*PL*PL+i*PL+j]=p22[a*PL+i][b*PL+j]^tmp[a*PL+i][b*PL+j];
    }
    memcpy(pk,seed,SEED_LENGTH_PUBLIC);compress_gf(pk+SEED_LENGTH_PUBLIC,raw,NUMGF_PK);
    memcpy(sk,seed,SEED_LENGTH_PUBLIC+SEED_LENGTH_PRIVATE);
#if HASH_PK
    shake256(sk+SEED_LENGTH_PUBLIC+SEED_LENGTH_PRIVATE,BYTES_PK_HASH,pk,BYTES_PK);
#endif
    SNOVA_SECURE_CLEAR(t);SNOVA_SECURE_CLEAR(tl);SNOVA_SECURE_CLEAR(f);
    SNOVA_SECURE_CLEAR(tmp);SNOVA_SECURE_CLEAR(p22);
#if R3_KEY_ORDER == 1
    SNOVA_SECURE_CLEAR(p21);
#endif
    return 0;
#endif
}

int SNOVA_NAMESPACE(sk_expand)(expanded_SK *skx,const uint8_t *sk) {
#if !R3_KEY
    return SNOVA_NAMESPACE(original_sk_expand)(skx,sk);
#else
    SNOVA_INIT
    memcpy(skx->pk_seed,sk,SEED_LENGTH_PUBLIC);
    memcpy(skx->sk_seed,sk+SEED_LENGTH_PUBLIC,SEED_LENGTH_PRIVATE);
#if HASH_PK
    memcpy(skx->pk_hash,sk+SEED_LENGTH_PUBLIC+SEED_LENGTH_PRIVATE,BYTES_PK_HASH);
#endif
    uint8_t t[VL][PAD(OL)],tl[OL][PAD(VL)],pub[NUM_PUB_GF];
    uint8_t p[VL][PAD(VL)],f12[VL][PAD(OL)],f21[OL][PAD(VL)];
    pt12(t,tl,skx->sk_seed);
    for(int v=0;v<PV;v++)for(int o=0;o<PO;o++)
        for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)
            skx->T12[(v*PO+o)*PL*PL+i*PL+j]=gf16_expand(t[v*PL+i][o*PL+j]);
    expand_public(pub,skx->pk_seed);
    for(int i=0;i<PM*PV*PV*PL*PL;i++)skx->P11[i]=pub[i];
    for(int m=0;m<PM;m++) {
        pblock8(&p[0][0],PAD(VL),pub+m*PV*PV*PL*PL,PV,PV);
        pmm(&f12[0][0],PAD(OL),&p[0][0],PAD(VL),&t[0][0],PAD(OL),VL,VL,OL);
        pmm(&f21[0][0],PAD(VL),&tl[0][0],PAD(VL),&p[0][0],PAD(VL),OL,VL,VL);
        punblock16(skx->F12+m*PV*PO*PL*PL,&f12[0][0],PAD(OL),PV,PO,
                   pub+PM*PV*PV*PL*PL+m*PV*PO*PL*PL);
        punblock16(skx->F21+m*PO*PV*PL*PL,&f21[0][0],PAD(VL),PO,PV,
                   pub+PM*PV*PN*PL*PL+m*PO*PV*PL*PL);
    }
    SNOVA_SECURE_CLEAR(t);SNOVA_SECURE_CLEAR(tl);
    SNOVA_SECURE_CLEAR(f12);SNOVA_SECURE_CLEAR(f21);
    return 0;
#endif
}

typedef struct {
    uint8_t a[PR][PAD(PR)],b[PR][PAD(PL)],q1t[PL][PAD(PL)],q2[PL][PAD(PL)],qp[PL*PL];
} pabq_raw;
static pabq_raw public_abq[PO*PAL];
static int public_abq_ready;
static void pabq_init(void) {
    if(public_abq_ready)return;
    for(int index=0;index<PO*PAL;index++) {
        pabq_raw *p=public_abq+index;
        for(int i=0;i<PR;i++)for(int j=0;j<PR;j++)p->a[i][j]=gf16_compress(fixedAm[index*PR*PR+i*PR+j]);
        for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)p->b[i][j]=gf16_compress(fixedBm[index*PR*PL+i*PL+j]);
        for(int i=0;i<PL;i++)for(int j=0;j<PL;j++) {
            p->q1t[i][j]=gf16_compress(fixedQ1[index*PL*PL+j*PL+i]);
            p->q2[i][j]=gf16_compress(fixedQ2[index*PL*PL+i*PL+j]);
            p->qp[i*PL+j]=pmul(gf16_compress(fixedq1[index*PL+i]),gf16_compress(fixedq2[index*PL+j]));
        }
    }
    public_abq_ready=1;
}
static void pabq(int index,uint8_t a[PR][PAD(PR)],uint8_t b[PR][PAD(PL)],
                 uint8_t q1t[PL][PAD(PL)],uint8_t q2[PL][PAD(PL)]) {
    pabq_init();
    memcpy(a,public_abq[index].a,PR*PAD(PR));memcpy(b,public_abq[index].b,PR*PAD(PL));
    memcpy(q1t,public_abq[index].q1t,PL*PAD(PL));memcpy(q2,public_abq[index].q2,PL*PAD(PL));
}

/* L = A X^T Q1^T, R = Q2 X B. Concatenate the variable blocks,
 * so subsequent contractions operate on contiguous rank rows. */
typedef struct {
    uint8_t xt[PR][PAD(NL)],ax[PR][PAD(NL)];
    uint8_t xx[PL][PAD(PN*PR)],qx[PL][PAD(PN*PR)];
    uint8_t ac[PL][PAD(PN*PR)],lc[PL][PAD(PN*PR)];
    uint8_t rc[PR][PAD(NL)],rr[PL][PAD(NL)];
} plr_scratch;
static void plr(uint8_t *left,int ls,uint8_t *right,int rs,const uint8_t *x,int count,plr_scratch *w,
                const uint8_t a[PR][PAD(PR)],const uint8_t b[PR][PAD(PL)],
                const uint8_t q1t[PL][PAD(PL)],const uint8_t q2[PL][PAD(PL)]) {
    uint8_t (*xt)[PAD(NL)]=w->xt,(*ax)[PAD(NL)]=w->ax;
    uint8_t (*xx)[PAD(PN*PR)]=w->xx,(*qx)[PAD(PN*PR)]=w->qx;
    memset(xt,0,sizeof(w->xt));memset(xx,0,sizeof(w->xx));
    for(int n=0;n<count;n++)for(int i=0;i<PL;i++)for(int j=0;j<PR;j++) {
        xt[j][n*PL+i]=x[n*PL*PR+i*PR+j];xx[i][n*PR+j]=x[n*PL*PR+i*PR+j];
    }
    pmm(&ax[0][0],PAD(NL),&a[0][0],PAD(PR),&xt[0][0],PAD(NL),PR,PR,count*PL);
    pmm(&qx[0][0],PAD(PN*PR),&q2[0][0],PAD(PL),&xx[0][0],PAD(PN*PR),PL,PL,count*PR);
    memset(w->ac,0,sizeof(w->ac));memset(w->rc,0,sizeof(w->rc));
    uint8_t q1[PL][PAD(PL)],bt[PL][PAD(PR)];
    ptrans(&q1[0][0],PAD(PL),&q1t[0][0],PAD(PL),PL,PL);
    ptrans(&bt[0][0],PAD(PR),&b[0][0],PAD(PL),PR,PL);
    for(int n=0;n<count;n++)for(int i=0;i<PR;i++)for(int k=0;k<PL;k++)w->ac[k][i*count+n]=ax[i][n*PL+k];
    for(int n=0;n<count;n++)for(int i=0;i<PL;i++)for(int k=0;k<PR;k++)w->rc[k][i*count+n]=qx[i][n*PR+k];
    pmm(&w->lc[0][0],sizeof(w->lc[0]),&q1[0][0],PAD(PL),&w->ac[0][0],sizeof(w->ac[0]),PL,PL,PR*count);
    pmm(&w->rr[0][0],sizeof(w->rr[0]),&bt[0][0],PAD(PR),&w->rc[0][0],sizeof(w->rc[0]),PL,PR,PL*count);
    memset(left,0,PR*ls);memset(right,0,count*PL*rs);
    for(int n=0;n<count;n++)for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)left[i*ls+n*PL+j]=w->lc[j][i*count+n];
    for(int n=0;n<count;n++)for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)right[(n*PL+i)*rs+j]=w->rr[j][i*count+n];
}

/* Direct Round2 evaluation: W_i=[S^0 X_i,...,S^(l-1) X_i],
 * U_m=W^T P_m W. All P blocks use the Round3 m1 index, not o. */
static void pdirect(uint8_t out[EQ],const uint16_t *p,const uint8_t *x,int count,int public_inputs) {
    uint8_t w[NL][PAD(PL*PR)]={0},wt[PL*PR][PAD(NL)];
    uint8_t pf[NL][PAD(NL)],right[NL][PAD(PL*PR)];
    uint8_t u[PM][PL*PR][PAD(PL*PR)];
    for(int n=0;n<count;n++)for(int i=0;i<PL;i++)for(int a=0;a<PL;a++)
        for(int j=0;j<PR;j++)for(int k=0;k<PL;k++)
            w[n*PL+i][a*PR+j]^=pmul((uint8_t)gf16_compress(gf_Sx[a*PL*PL+i*PL+k]),x[n*PL*PR+k*PR+j]);
    ptrans(&wt[0][0],PAD(NL),&w[0][0],PAD(PL*PR),count*PL,PL*PR);
#if R3_METHOD == 5 || R3_METHOD == 7 || R3_METHOD == 8
    uint8_t wp[NL][PAD((PL*PR+1)/2)]={0},rp[NL][PAD((PL*PR+1)/2)],upack[PL*PR][PAD((PL*PR+1)/2)];
    for(int i=0;i<count*PL;i++)for(int j=0;j<PL*PR;j++)wp[i][j/2]|=w[i][j]<<((j&1)*4);
#if R3_METHOD == 8
    uint8_t wtp[(PL*PR+1)/2][PAD(NL)],pp[(NL+1)/2][PAD(NL)];
    ppack_a(&wtp[0][0],PAD(NL),&wt[0][0],PAD(NL),PL*PR,count*PL);
#endif
#endif
    for(int m=0;m<PM;m++) {
        pblock16(&pf[0][0],PAD(NL),p+m*count*count*PL*PL,count,count);
#if R3_METHOD == 5 || R3_METHOD == 7 || R3_METHOD == 8
#if R3_METHOD == 8
        ppack_a(&pp[0][0],PAD(NL),&pf[0][0],PAD(NL),count*PL,count*PL);
        pmm_ab_split(&rp[0][0],sizeof(rp[0]),&pp[0][0],PAD(NL),&wp[0][0],sizeof(wp[0]),count*PL,count*PL,(PL*PR+1)/2);
        pmm_ab_split(&upack[0][0],sizeof(upack[0]),&wtp[0][0],PAD(NL),&rp[0][0],sizeof(rp[0]),PL*PR,count*PL,(PL*PR+1)/2);
#else
        pmm_packed(&rp[0][0],sizeof(rp[0]),&pf[0][0],PAD(NL),&wp[0][0],sizeof(wp[0]),count*PL,count*PL,(PL*PR+1)/2,1);
        pmm_packed(&upack[0][0],sizeof(upack[0]),&wt[0][0],PAD(NL),&rp[0][0],sizeof(rp[0]),PL*PR,count*PL,(PL*PR+1)/2,public_inputs);
#endif
        memset(u[m],0,sizeof(u[m]));
        for(int i=0;i<PL*PR;i++)for(int j=0;j<PL*PR;j++)u[m][i][j]=(upack[i][j/2]>>((j&1)*4))&15;
#else
        pmm_public(&right[0][0],PAD(PL*PR),&pf[0][0],PAD(NL),&w[0][0],PAD(PL*PR),count*PL,count*PL,PL*PR);
        if(public_inputs)pmm_public(&u[m][0][0],PAD(PL*PR),&wt[0][0],PAD(NL),&right[0][0],PAD(PL*PR),PL*PR,count*PL,PL*PR);
        else pmm(&u[m][0][0],PAD(PL*PR),&wt[0][0],PAD(NL),&right[0][0],PAD(PL*PR),PL*PR,count*PL,PL*PR);
#endif
    }
    uint8_t up[PM][PL*PL][PAD(PR*PR)]={0};
    for(int m=0;m<PM;m++)for(int c=0;c<PL;c++)for(int d=0;d<PL;d++)
        for(int i=0;i<PR;i++)for(int j=0;j<PR;j++)up[m][c*PL+d][i*PR+j]=u[m][c*PR+i][d*PR+j];
    memset(out,0,EQ);
#if R3_PUBLIC_E
    pe_init();
    uint8_t uf[PAD(PM*PL*PL*PR*PR)]={0},ef[PAD(EQ)];
    for(int m=0;m<PM;m++)for(int cd=0;cd<PL*PL;cd++)for(int ij=0;ij<PR*PR;ij++)
        uf[(m*PL*PL+cd)*PR*PR+ij]=up[m][cd][ij];
    pmm(ef,PAD(EQ),uf,sizeof(uf),&pe_map[0][0],PAD(EQ),1,PM*PL*PL*PR*PR,EQ);
    memcpy(out,ef,EQ);SNOVA_SECURE_CLEAR(uf);SNOVA_SECURE_CLEAR(ef);
#else
    for(int m=0;m<PO;m++)for(int al=0;al<PAL;al++) {
        int ix=m*PAL+al,mp=i_prime(m,al);
        uint8_t g[PR][PAD(PR)]={0},gb[PR][PAD(PL)],ag[PR][PAD(PL)];
        uint8_t a[PR][PAD(PR)],b[PR][PAD(PL)],q1t[PL][PAD(PL)],q2[PL][PAD(PL)];
        pabq(ix,a,b,q1t,q2);
        uint8_t packedg[PAD(PR*PR)];
        for(int k=0;k<PAD(PR*PR);k+=PV_BYTES) {
            pv acc=pz();
            for(int c=0;c<PL;c++)for(int d=0;d<PL;d++) {
                uint8_t coeff=public_abq[ix].qp[c*PL+d];
                acc=px(acc,pconst(pl(up[mp][c*PL+d]+k),coeff));
            }
            pst(packedg+k,acc);
        }
        for(int i=0;i<PR;i++)for(int j=0;j<PR;j++)g[i][j]=packedg[i*PR+j];
        pmm(&gb[0][0],PAD(PL),&g[0][0],PAD(PR),&b[0][0],PAD(PL),PR,PR,PL);
        pmm(&ag[0][0],PAD(PL),&a[0][0],PAD(PR),&gb[0][0],PAD(PL),PR,PR,PL);
        for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)out[m*PR*PL+i*PL+j]^=ag[i][j];
        SNOVA_SECURE_CLEAR(g);SNOVA_SECURE_CLEAR(gb);SNOVA_SECURE_CLEAR(ag);
        SNOVA_SECURE_CLEAR(packedg);
    }
#endif
    SNOVA_SECURE_CLEAR(w);SNOVA_SECURE_CLEAR(wt);SNOVA_SECURE_CLEAR(right);SNOVA_SECURE_CLEAR(u);
    SNOVA_SECURE_CLEAR(up);
#if R3_METHOD == 5 || R3_METHOD == 7 || R3_METHOD == 8
    SNOVA_SECURE_CLEAR(wp);SNOVA_SECURE_CLEAR(rp);SNOVA_SECURE_CLEAR(upack);
#if R3_METHOD == 8
    SNOVA_SECURE_CLEAR(wtp);
#endif
#endif
}

static void pstructured(uint8_t out[EQ],const uint16_t *p,const uint8_t *x,int count) {
    uint8_t pf[PM][NL][PAD(NL)];
    plr_scratch work;
    uint8_t left[PR][PAD(NL)],right[NL][PAD(PL)],lp[PR][PAD(NL)],h[PR][PAD(PL)];
    for(int m=0;m<PM;m++)pblock16(&pf[m][0][0],PAD(NL),p+m*count*count*PL*PL,count,count);
    memset(out,0,EQ);
    for(int m=0;m<PO;m++)for(int al=0;al<PAL;al++) {
        uint8_t a[PR][PAD(PR)],b[PR][PAD(PL)],q1t[PL][PAD(PL)],q2[PL][PAD(PL)];
        pabq(m*PAL+al,a,b,q1t,q2);
        plr(&left[0][0],PAD(NL),&right[0][0],PAD(PL),x,count,&work,a,b,q1t,q2);
        pmm(&lp[0][0],PAD(NL),&left[0][0],PAD(NL),&pf[i_prime(m,al)][0][0],PAD(NL),PR,count*PL,count*PL);
        pmm(&h[0][0],PAD(PL),&lp[0][0],PAD(NL),&right[0][0],PAD(PL),PR,count*PL,PL);
        for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)out[m*PL*PR+i*PL+j]^=h[i][j];
    }
    SNOVA_SECURE_CLEAR(left);SNOVA_SECURE_CLEAR(right);SNOVA_SECURE_CLEAR(lp);SNOVA_SECURE_CLEAR(h);
    SNOVA_SECURE_CLEARI(work);
}

static int psolve(uint8_t g[EQ][GS],uint8_t solution[EQ]) {
    int singular=0;
    for(int i=0;i<EQ;i++) {
        for(int j=i+1;j<EQ;j++) {
            pv mask=ps((uint8_t)(ct_is_not_zero(g[i][i])-1));
            for(int k=(i/PV_BYTES)*PV_BYTES;k<GS;k+=PV_BYTES)
                pst(g[i]+k,px(pl(g[i]+k),pa(mask,pl(g[j]+k))));
        }
        uint8_t x=g[i][i],x2=pmul(x,x),x4=pmul(x2,x2),x8=pmul(x4,x4);
        uint8_t inv=pmul(pmul(x2,x4),x8);
        singular|=1-ct_is_not_zero(x);
        pv tab=ptab(inv);
        for(int k=(i/PV_BYTES)*PV_BYTES;k<GS;k+=PV_BYTES)pst(g[i]+k,shuf(tab,pl(g[i]+k)));
        for(int j=i+1;j<EQ;j++) {
            tab=ptab(g[j][i]);
            for(int k=(i/PV_BYTES)*PV_BYTES;k<GS;k+=PV_BYTES)
                pst(g[j]+k,px(pl(g[j]+k),shuf(tab,pl(g[i]+k))));
        }
    }
    memset(solution,0,EQ);
    for(int i=EQ-1;i>=0;i--) {
        uint8_t x=g[i][EQ];
        for(int j=i+1;j<EQ;j++)x^=pmul(g[i][j],solution[j]);
        solution[i]=x;
    }
    return singular;
}

int SNOVA_NAMESPACE(sign)(const expanded_SK *skx,uint8_t *sig,const uint8_t *digest,size_t len_digest,const uint8_t *salt) {
#if !R3_SIGN
    return SNOVA_NAMESPACE(original_sign)(skx,sig,digest,len_digest,salt);
#else
    SNOVA_INIT
    uint8_t target[GF16_HASH],hb[BYTES_HASH],x[PN*PL*PR];
#if HASH_PK
    hash_combined(hb,digest,len_digest,skx->pk_hash,salt);
#else
    hash_combined(hb,digest,len_digest,skx->pk_seed,salt);
#endif
    convert_bytes_to_GF(target,hb,GF16_HASH);
    uint8_t f12[PM][VL][PAD(OL)],f21[PM][OL][PAD(VL)];
#if R3_FVV == 2
    uint8_t f11[PM][VL][PAD(VL)],lp[PR][PAD(VL)],quad[PR][PAD(PL)];
#endif
    for(int m=0;m<PM;m++) {
        pblock16(&f12[m][0][0],PAD(OL),skx->F12+m*PV*PO*PL*PL,PV,PO);
        pblock16(&f21[m][0][0],PAD(VL),skx->F21+m*PO*PV*PL*PL,PO,PV);
#if R3_FVV == 2
        pblock16(&f11[m][0][0],PAD(VL),skx->P11+m*PV*PV*PL*PL,PV,PV);
#endif
    }
    int result=-1;
    plr_scratch work;
    uint8_t left[PR][PAD(VL)],right[VL][PAD(PL)],h12[PR][PAD(OL)],h21[OL][PAD(PL)];
    uint8_t c12[PR][PAD(OL)],c21[OL][PAD(PL)],c21t[PL][PAD(OL)];
    uint8_t hc12[PL][PAD(PO*PR)],cc12[PL][PAD(PO*PR)],hc21[PL][PAD(OL)],cc21[PL][PAD(OL)];
    for(int attempt=1;attempt<255;attempt++) {
        uint8_t attempt_byte=(uint8_t)attempt,vinegar[NUM_GEN_SEC_BYTES],fv[EQ],solution[EQ];
        uint8_t g[EQ][GS]={0},gacc[PO][PR][PL][PR][PAD(OL)]={0};shake_t state;
        shake256_init(&state);shake_absorb(&state,skx->sk_seed,SEED_LENGTH_PRIVATE);
        shake_absorb(&state,digest,BYTES_DIGEST);shake_absorb(&state,salt,BYTES_SALT);
        shake_absorb(&state,&attempt_byte,1);shake_finalize(&state);
        shake_squeeze(vinegar,NUM_GEN_SEC_BYTES,&state);
        convert_bytes_to_GF(x,vinegar,PV*PL*PR);
#if R3_FVV == 2
        memset(fv,0,sizeof(fv));
#else
        pdirect(fv,skx->P11,x,PV,0);
#endif
        for(int m=0;m<PO;m++)for(int al=0;al<PAL;al++) {
            int mp=i_prime(m,al);
            uint8_t a[PR][PAD(PR)],b[PR][PAD(PL)],q1t[PL][PAD(PL)],q2[PL][PAD(PL)];
            memset(c12,0,sizeof(c12));memset(c21,0,sizeof(c21));
            pabq(m*PAL+al,a,b,q1t,q2);
            plr(&left[0][0],PAD(VL),&right[0][0],PAD(PL),x,PV,&work,a,b,q1t,q2);
#if R3_FVV == 2
            pmm(&lp[0][0],PAD(VL),&left[0][0],PAD(VL),&f11[mp][0][0],PAD(VL),PR,VL,VL);
            pmm(&quad[0][0],PAD(PL),&lp[0][0],PAD(VL),&right[0][0],PAD(PL),PR,VL,PL);
            for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)fv[m*PR*PL+i*PL+j]^=quad[i][j];
#endif
            pmm(&h12[0][0],PAD(OL),&left[0][0],PAD(VL),&f12[mp][0][0],PAD(OL),PR,VL,OL);
            pmm(&h21[0][0],PAD(PL),&f21[mp][0][0],PAD(VL),&right[0][0],PAD(PL),OL,VL,PL);
            uint8_t q2t[PL][PAD(PL)];
            ptrans(&q2t[0][0],PAD(PL),&q2[0][0],PAD(PL),PL,PL);
            memset(hc12,0,sizeof(hc12));memset(hc21,0,sizeof(hc21));
            for(int o=0;o<PO;o++)for(int i=0;i<PR;i++)for(int k=0;k<PL;k++)hc12[k][o*PR+i]=h12[i][o*PL+k];
            for(int o=0;o<PO;o++)for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)hc21[i][o*PL+j]=h21[o*PL+i][j];
            pmm(&cc12[0][0],sizeof(cc12[0]),&q2t[0][0],PAD(PL),&hc12[0][0],sizeof(hc12[0]),PL,PL,PO*PR);
            pmm(&cc21[0][0],sizeof(cc21[0]),&q1t[0][0],PAD(PL),&hc21[0][0],sizeof(hc21[0]),PL,PL,OL);
            memset(c21t,0,sizeof(c21t));
            for(int o=0;o<PO;o++)for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)c12[i][o*PL+j]=cc12[j][o*PR+i];
            for(int o=0;o<PO;o++)for(int i=0;i<PL;i++)for(int j=0;j<PL;j++)c21t[j][o*PL+i]=cc21[i][o*PL+j];
            for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)for(int u=0;u<PR;u++)
                for(int k=0;k<PAD(OL);k+=PV_BYTES) {
                    pv acc=pl(gacc[m][i][j][u]+k);
                    acc=px(acc,pconst(pl(c12[i]+k),b[u][j]));
                    acc=px(acc,pconst(pl(c21t[j]+k),a[i][u]));
                    pst(gacc[m][i][j][u]+k,acc);
                }
        }
        for(int m=0;m<PO;m++)for(int i=0;i<PR;i++)for(int j=0;j<PL;j++)
            for(int o=0;o<PO;o++)for(int k=0;k<PL;k++)for(int u=0;u<PR;u++)
                g[m*PL*PR+i*PL+j][o*PL*PR+k*PR+u]=gacc[m][i][j][u][o*PL+k];
        for(int i=0;i<EQ;i++)g[i][EQ]=target[i]^fv[i];
        int singular=psolve(g,solution);
        if(!singular) {
            memcpy(x+PV*PL*PR,solution,EQ);
            for(int v=0;v<PV;v++)for(int o=0;o<PO;o++)for(int i=0;i<PL;i++)for(int j=0;j<PR;j++)
                for(int k=0;k<PL;k++)x[v*PL*PR+i*PR+j]^=pmul(gf16_compress(skx->T12[(v*PO+o)*PL*PL+i*PL+k]),solution[o*PL*PR+k*PR+j]);
            compress_gf(sig,x,PN*PL*PR);memcpy(sig+BYTES_SIGNATURE-BYTES_SALT,salt,BYTES_SALT);result=0;
        }
        SNOVA_SECURE_CLEAR(g);SNOVA_SECURE_CLEAR(vinegar);SNOVA_SECURE_CLEAR(fv);
        SNOVA_SECURE_CLEAR(gacc);
        SNOVA_SECURE_CLEAR(solution);SNOVA_SECURE_CLEARI(state);
        if(!singular)break;
    }
    if(result)memset(sig,0,BYTES_SIGNATURE);
    SNOVA_SECURE_CLEAR(x);SNOVA_SECURE_CLEAR(f12);SNOVA_SECURE_CLEAR(f21);
    SNOVA_SECURE_CLEAR(left);SNOVA_SECURE_CLEAR(right);SNOVA_SECURE_CLEAR(h12);SNOVA_SECURE_CLEAR(h21);
    SNOVA_SECURE_CLEAR(c12);SNOVA_SECURE_CLEAR(c21);SNOVA_SECURE_CLEAR(c21t);SNOVA_SECURE_CLEARI(work);
    SNOVA_SECURE_CLEAR(hc12);SNOVA_SECURE_CLEAR(cc12);SNOVA_SECURE_CLEAR(hc21);SNOVA_SECURE_CLEAR(cc21);
#if R3_FVV == 2
    SNOVA_SECURE_CLEAR(lp);SNOVA_SECURE_CLEAR(quad);
#endif
    return result;
#endif
}

int SNOVA_NAMESPACE(verify)(const expanded_PK *pkx,const uint8_t *sig,const uint8_t *digest,size_t len_digest) {
#if !R3_VERIFY
    return SNOVA_NAMESPACE(original_verify)(pkx,sig,digest,len_digest);
#else
    SNOVA_INIT
    uint8_t x[PN*PL*PR],out[EQ],target[GF16_HASH],hb[BYTES_HASH];
    convert_bytes_to_GF(x,sig,PN*PL*PR);
#if R3_VERIFY == 2
    pstructured(out,pkx->P,x,PN);
#else
    pdirect(out,pkx->P,x,PN,1);
#endif
#if HASH_PK
    hash_combined(hb,digest,len_digest,pkx->pk_hash,sig+BYTES_SIGNATURE-BYTES_SALT);
#else
    hash_combined(hb,digest,len_digest,pkx->pk_seed,sig+BYTES_SIGNATURE-BYTES_SALT);
#endif
    convert_bytes_to_GF(target,hb,GF16_HASH);
    unsigned diff=0;for(int i=0;i<GF16_HASH;i++)diff|=out[i]^target[i];
    return diff?-1:0;
#endif
}
