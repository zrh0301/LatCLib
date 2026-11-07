#ifndef R3_ROWS
#define R3_ROWS 2
#endif
#ifndef R3_CHAINS
#define R3_CHAINS 1
#endif
#if defined(__GFNI__)
/* Several output rows share each B load. Independent inner chains are XORed
 * before the single polynomial reduction; all scratch stays in this call. */
static void pmm_row_tile(uint8_t *c,int cs,const uint8_t *a,int as,
                         const uint8_t *b,int bs,int rows,int inner,int cols) {
    for(int i=0;i<rows;i+=R3_ROWS)for(int j=0;j<PAD(cols);j+=PV_BYTES) {
        pv acc[R3_ROWS][R3_CHAINS];
        for(int r=0;r<R3_ROWS;r++)for(int t=0;t<R3_CHAINS;t++)acc[r][t]=pz();
        for(int k=0;k<inner;k+=R3_CHAINS) {
            for(int t=0;t<R3_CHAINS;t++)if(k+t<inner) {
                pv v=pl(b+(k+t)*bs+j);
                for(int r=0;r<R3_ROWS;r++)if(i+r<rows)
                    acc[r][t]=px(acc[r][t],_mm256_gf2p8mul_epi8(ps(a[(i+r)*as+k+t]),v));
            }
        }
        pv reduce=_mm256_set1_epi64x(0x113264c800000000ULL);
        for(int r=0;r<R3_ROWS;r++)if(i+r<rows) {
            pv sum=acc[r][0];for(int t=1;t<R3_CHAINS;t++)sum=px(sum,acc[r][t]);
            pst(c+(i+r)*cs+j,_mm256_gf2p8affine_epi64_epi8(sum,reduce,0));
        }
    }
}
#endif
