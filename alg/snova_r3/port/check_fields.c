#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define SNOVA_l 4
#define SNOVA_l2 16
#define gf16_expand(a) (((a)|((a)<<3)|((a)<<6)|((a)<<9))&0x1111)
#define gf16_compress(a) (((a)&15)^(((a)>>3)&14)^(((a)>>6)&14)^(((a)>>7)&1)^(((a)>>9)&14)^(((a)>>10)&3)^(((a)>>12)&14)^(((a)>>13)&7))
static void erase(void *p,size_t n) { volatile uint8_t *q=p;while(n--)*q++=0; }
#define SNOVA_SECURE_CLEAR(a) erase(a,sizeof(a))
#include "vector.h"
static uint8_t oracle(uint8_t a,uint8_t b) {
    unsigned r=0;
    for(int k=0;k<4;k++){if(b&1)r^=a;b>>=1;a=(a<<1)^((a&8)?19:0);}
    return (uint8_t)r;
}
#define CHECK(x) do{if(!(x)){fprintf(stderr,"field correctness failed at line %d method=%d\n",__LINE__,R3_METHOD);return 1;}}while(0)
int main(void) {
    uint8_t av[PV_BYTES],bv[PV_BYTES],cv[PV_BYTES];
    for(int a=0;a<16;a++) {
        for(int i=0;i<PV_BYTES;i++){av[i]=a;bv[i]=i&15;}
        pst(cv,pmulv(pl(av),pl(bv)));
        for(int i=0;i<PV_BYTES;i++){CHECK(cv[i]==oracle(a,bv[i]));CHECK(pmul(a,bv[i])==cv[i]);}
    }
    for(int pair=0;pair<256;pair++) {
        for(int i=0;i<PV_BYTES;i++)av[i]=i&15;
        pst(cv,shuf(ppairtab(pair&15,pair>>4),pl(av)));
        for(int i=0;i<PV_BYTES;i++)CHECK(cv[i]==(oracle(pair&15,av[i])|(oracle(pair>>4,av[i])<<4)));
    }
    uint8_t pb[256],pc[256];for(int i=0;i<256;i++)pb[i]=(uint8_t)i;
    for(uint8_t a=0;a<16;a++)for(int pub=0;pub<2;pub++) {
        pmm_packed(pc,256,&a,1,pb,256,1,1,256,pub);
        for(int i=0;i<256;i++)CHECK(pc[i]==(oracle(a,i&15)|(oracle(a,i>>4)<<4)));
    }
    const int cols_list[]={1,3,15,16,17,31,32,33,61};
    const int rows_list[]={1,2,3,5},inner_list[]={1,5,33};
    uint8_t a[5][64],b[33][64],guard[5*64+2];
    for(int i=0;i<5;i++)for(int j=0;j<64;j++)a[i][j]=(i*7+j*3+1)&15;
    int cases=0;
    for(unsigned rr=0;rr<sizeof(rows_list)/sizeof(*rows_list);rr++)
    for(unsigned kk=0;kk<sizeof(inner_list)/sizeof(*inner_list);kk++)
    for(unsigned cc=0;cc<sizeof(cols_list)/sizeof(*cols_list);cc++) {
        int rows=rows_list[rr],inner=inner_list[kk],cols=cols_list[cc],stride=PAD(cols);
        memset(b,0,sizeof(b));memset(guard,0xa5,sizeof(guard));
        for(int i=0;i<inner;i++)for(int j=0;j<cols;j++)b[i][j]=(i*5+j*9+2)&15;
        pmm(guard+1,stride,&a[0][0],64,&b[0][0],64,rows,inner,cols);
        CHECK(guard[0]==0xa5);CHECK(guard[rows*stride+1]==0xa5);
        for(int i=0;i<rows;i++)for(int j=0;j<cols;j++) {
            uint8_t expected=0;for(int k=0;k<inner;k++)expected^=oracle(a[i][k],b[k][j]);
            CHECK(guard[1+i*stride+j]==expected);
        }
        cases++;
    }
    printf("field_correctness=PASS method=%d matrix_shapes=%d counter=not-measured\n",R3_METHOD,cases);
    return 0;
}
