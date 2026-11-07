/* One-sided expanded-u16 multiplication, XOR accumulation, one compression
 * per output vector. Raw B and padded row layout match the other kernels. */
static void pmm_expanded(uint8_t *c,int cs,const uint8_t *a,int as,
                         const uint8_t *b,int bs,int rows,int inner,int cols) {
    for(int i=0;i<rows;i++)for(int j=0;j<PAD(cols);j+=PV_BYTES/2) {
#if defined(__AVX2__)
        __m256i acc=_mm256_setzero_si256();
        for(int k=0;k<inner;k++) {
            __m256i bv=_mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i*)(b+k*bs+j)));
            acc=_mm256_xor_si256(acc,_mm256_mullo_epi16(_mm256_set1_epi16(gf16_expand(a[i*as+k])),bv));
        }
        __m256i p=_mm256_xor_si256(_mm256_and_si256(acc,_mm256_set1_epi16(15)),
                   _mm256_and_si256(_mm256_srli_epi16(acc,3),_mm256_set1_epi16(30)));
        p=_mm256_xor_si256(p,_mm256_and_si256(_mm256_srli_epi16(acc,6),_mm256_set1_epi16(60)));
        p=_mm256_xor_si256(p,_mm256_and_si256(_mm256_srli_epi16(acc,9),_mm256_set1_epi16(120)));
        __m256i hi=_mm256_srli_epi16(p,4);
        p=_mm256_and_si256(_mm256_xor_si256(p,_mm256_xor_si256(hi,_mm256_slli_epi16(hi,1))),_mm256_set1_epi16(15));
        p=_mm256_permute4x64_epi64(_mm256_packus_epi16(p,_mm256_setzero_si256()),0xd8);
        _mm_storeu_si128((__m128i*)(c+i*cs+j),_mm256_castsi256_si128(p));
#else
        uint16x8_t acc=vdupq_n_u16(0);
        for(int k=0;k<inner;k++) {
            uint16x8_t bv=vmovl_u8(vld1_u8(b+k*bs+j));
            acc=veorq_u16(acc,vmulq_u16(vdupq_n_u16(gf16_expand(a[i*as+k])),bv));
        }
        uint16x8_t p=veorq_u16(vandq_u16(acc,vdupq_n_u16(15)),vandq_u16(vshrq_n_u16(acc,3),vdupq_n_u16(30)));
        p=veorq_u16(p,vandq_u16(vshrq_n_u16(acc,6),vdupq_n_u16(60)));
        p=veorq_u16(p,vandq_u16(vshrq_n_u16(acc,9),vdupq_n_u16(120)));
        uint16x8_t hi=vshrq_n_u16(p,4);
        p=vandq_u16(veorq_u16(p,veorq_u16(hi,vshlq_n_u16(hi,1))),vdupq_n_u16(15));
        vst1_u8(c+i*cs+j,vmovn_u16(p));
#endif
    }
}

/* Four independent left coefficients share B, 2B, 4B, 8B. */
static void pmm_sameb4(uint8_t *c,int cs,const uint8_t *a,int as,
                       const uint8_t *b,int bs,int rows,int inner,int cols) {
    for(int i=0;i<rows;i+=4)for(int j=0;j<PAD(cols);j+=PV_BYTES) {
        pv acc[4]={pz(),pz(),pz(),pz()};
        for(int k=0;k<inner;k++) {
            pv b1=pl(b+k*bs+j),b2=pxtime(b1),b4=pxtime(b2),b8=pxtime(b4);
            for(int r=0;r<4 && i+r<rows;r++) {
                unsigned s=a[(i+r)*as+k];
                pv t1=pa(b1,ps((uint8_t)-(int)(s&1)));
                pv t2=pa(b2,ps((uint8_t)-(int)((s>>1)&1)));
                pv t4=pa(b4,ps((uint8_t)-(int)((s>>2)&1)));
                pv t8=pa(b8,ps((uint8_t)-(int)((s>>3)&1)));
                acc[r]=px(acc[r],px(px(t1,t2),px(t4,t8)));
            }
        }
        for(int r=0;r<4 && i+r<rows;r++)pst(c+(i+r)*cs+j,acc[r]);
    }
}
