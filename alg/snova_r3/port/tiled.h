#ifndef R3_TILE
#define R3_TILE 4
#endif
/* Two output rows by R3_TILE vector columns. Independent accumulators
 * expose throughput without changing the order of field XOR reduction. */
static void pmm_tile(uint8_t *c,int cs,const uint8_t *a,int as,
                     const uint8_t *b,int bs,int rows,int inner,int cols) {
    const int stride=PAD(cols);
    for(int i=0;i<rows;i+=2)for(int j=0;j<stride;j+=R3_TILE*PV_BYTES) {
        pv s0[R3_TILE],s1[R3_TILE];
        for(int g=0;g<R3_TILE;g++){s0[g]=pz();s1[g]=pz();}
        for(int k=0;k<inner;k++) {
#if R3_METHOD == 0
            pv table=ppairtab(a[i*as+k],i+1<rows?a[(i+1)*as+k]:0);
#else
            pv a0=ps(a[i*as+k]),a1=ps(i+1<rows?a[(i+1)*as+k]:0);
#endif
            for(int g=0;g<R3_TILE;g++)if(j+g*PV_BYTES<stride) {
                pv v=pl(b+k*bs+j+g*PV_BYTES);
#if R3_METHOD == 0
                s0[g]=px(s0[g],shuf(table,v));
#elif defined(__GFNI__) && R3_METHOD == 1
                s0[g]=px(s0[g],_mm256_gf2p8mul_epi8(a0,v));
                s1[g]=px(s1[g],_mm256_gf2p8mul_epi8(a1,v));
#else
                s0[g]=px(s0[g],pmulv(a0,v));
                s1[g]=px(s1[g],pmulv(a1,v));
#endif
            }
        }
        for(int g=0;g<R3_TILE;g++)if(j+g*PV_BYTES<stride) {
#if R3_METHOD == 0
            s1[g]=shr4(s0[g]);s0[g]=pa(s0[g],ps(15));
#elif defined(__GFNI__) && R3_METHOD == 1
            pv reduce=_mm256_set1_epi64x(0x113264c800000000ULL);
            s0[g]=_mm256_gf2p8affine_epi64_epi8(s0[g],reduce,0);
            s1[g]=_mm256_gf2p8affine_epi64_epi8(s1[g],reduce,0);
#endif
            pst(c+i*cs+j+g*PV_BYTES,s0[g]);
            if(i+1<rows)pst(c+(i+1)*cs+j+g*PV_BYTES,s1[g]);
        }
    }
}
