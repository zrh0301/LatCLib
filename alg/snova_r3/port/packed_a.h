/* Two A rows share each byte. Each row produces its own packed-B output;
 * this is a 2-by-2 outer-product tile, not nibblewise diagonal multiplication. */
static void ppack_a(uint8_t *out,int os,const uint8_t *a,int as,int rows,int inner) {
    memset(out,0,(size_t)((rows+1)/2)*os);
    for(int i=0;i<rows;i+=2)for(int k=0;k<inner;k++)
        out[(i/2)*os+k]=a[i*as+k]|((i+1<rows?a[(i+1)*as+k]:0)<<4);
}
static inline pv pselect_planes(pv a,pv b1,pv b2,pv b4,pv b8) {
    pv t1=pa(b1,pe(pa(a,ps(1)),ps(1)));
    pv t2=pa(b2,pe(pa(a,ps(2)),ps(2)));
    pv t4=pa(b4,pe(pa(a,ps(4)),ps(4)));
    pv t8=pa(b8,pe(pa(a,ps(8)),ps(8)));
    return px(px(t1,t2),px(t4,t8));
}
static void pmm_ab_split(uint8_t *c,int cs,const uint8_t *ap,int aps,
                         const uint8_t *bp,int bps,int rows,int inner,int packed_cols) {
    for(int i=0;i<rows;i+=4)for(int j=0;j<PAD(packed_cols);j+=PV_BYTES) {
        pv acc[4]={pz(),pz(),pz(),pz()};
        for(int k=0;k<inner;k++) {
            pv b1=pl(bp+k*bps+j),b2=ppacked_xtime(b1),b4=ppacked_xtime(b2),b8=ppacked_xtime(b4);
            for(int pair=0;pair<2;pair++)if(i+2*pair<rows) {
                pv packed=ps(ap[(i/2+pair)*aps+k]);
                pv lo=pa(packed,ps(15)),hi=shr4(packed);
                acc[2*pair]=px(acc[2*pair],pselect_planes(lo,b1,b2,b4,b8));
                acc[2*pair+1]=px(acc[2*pair+1],pselect_planes(hi,b1,b2,b4,b8));
            }
        }
        for(int r=0;r<4;r++)if(i+r<rows)pst(c+(i+r)*cs+j,acc[r]);
    }
}
