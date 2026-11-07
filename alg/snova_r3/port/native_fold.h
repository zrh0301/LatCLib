/* Raw bytes at the existing native core boundary; no protocol/API changes. */
#define gf16_expand(a) (((a)|((a)<<3)|((a)<<6)|((a)<<9))&0x1111)
#define gf16_compress(a) (((a)&15)^(((a)>>3)&14)^(((a)>>6)&14)^(((a)>>7)&1)^(((a)>>9)&14)^(((a)>>10)&3)^(((a)>>12)&14)^(((a)>>13)&7))
#ifndef R3_NATIVE_FOLD
#if defined(__GFNI__) && SNOVA_l == 5 && SNOVA_o >= 5
#define R3_NATIVE_FOLD 1
#else
#define R3_NATIVE_FOLD 0
#endif
#endif
#ifndef R3_NATIVE_KEY_FOLD
#define R3_NATIVE_KEY_FOLD 0
#endif
#ifndef R3_FULL_KEY
#define R3_FULL_KEY 1
#endif
#if R3_NATIVE_FOLD
#ifndef R3_ROWS
#define R3_ROWS 4
#endif
#ifndef R3_CHAINS
#define R3_CHAINS 2
#endif
#endif
#include "vector.h"
static void pwrite_blocks(uint8_t *out,const uint8_t *in,int stride,int rows,int cols,int add) {
    for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        for(int a=0;a<SNOVA_l;a++)for(int b=0;b<SNOVA_l;b++) {
            int ix=(i*cols+j)*SNOVA_l2+a*SNOVA_l+b;
            uint8_t value=in[(i*SNOVA_l+a)*stride+j*SNOVA_l+b];
            if(add)out[ix]^=value;else out[ix]=value;
        }
}
static void p_fold_b(uint8_t *out,const uint8_t *in,const uint8_t *t,int rows,int add) {
    const int vl=SNOVA_v*SNOVA_l,ol=SNOVA_o*SNOVA_l;
    uint8_t tf[vl][PAD(ol)],af[rows*SNOVA_l][PAD(vl)],cf[rows*SNOVA_l][PAD(ol)];
    pblock8(&tf[0][0],PAD(ol),t,SNOVA_v,SNOVA_o);
    for(int m=0;m<SNOVA_m1;m++) {
        pblock8(&af[0][0],PAD(vl),in+m*rows*SNOVA_v*SNOVA_l2,rows,SNOVA_v);
        pmm(&cf[0][0],PAD(ol),&af[0][0],PAD(vl),&tf[0][0],PAD(ol),rows*SNOVA_l,vl,ol);
        pwrite_blocks(out+m*rows*SNOVA_o*SNOVA_l2,&cf[0][0],PAD(ol),rows,SNOVA_o,add);
    }
    SNOVA_SECURE_CLEAR(tf);SNOVA_SECURE_CLEAR(af);SNOVA_SECURE_CLEAR(cf);
}
static void p_fold_a(uint8_t *out,const uint8_t *t,const uint8_t *in,int cols,int add) {
    const int vl=SNOVA_v*SNOVA_l,ol=SNOVA_o*SNOVA_l,cl=cols*SNOVA_l;
    uint8_t tf[ol][PAD(vl)],bf[vl][PAD(cl)],cf[ol][PAD(cl)];
    memset(tf,0,sizeof(tf));
    for(int v=0;v<SNOVA_v;v++)for(int o=0;o<SNOVA_o;o++)
        for(int i=0;i<SNOVA_l;i++)for(int j=0;j<SNOVA_l;j++)
            tf[o*SNOVA_l+i][v*SNOVA_l+j]=t[(v*SNOVA_o+o)*SNOVA_l2+i*SNOVA_l+j];
    for(int m=0;m<SNOVA_m1;m++) {
        pblock8(&bf[0][0],PAD(cl),in+m*SNOVA_v*cols*SNOVA_l2,SNOVA_v,cols);
        pmm(&cf[0][0],PAD(cl),&tf[0][0],PAD(vl),&bf[0][0],PAD(cl),ol,vl,cl);
        pwrite_blocks(out+m*SNOVA_o*cols*SNOVA_l2,&cf[0][0],PAD(cl),SNOVA_o,cols,add);
    }
    SNOVA_SECURE_CLEAR(tf);SNOVA_SECURE_CLEAR(bf);SNOVA_SECURE_CLEAR(cf);
}
#undef gf16_expand
#undef gf16_compress
