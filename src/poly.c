#include "../inc/types.h"
#include "../inc/rand.h"
#include "../inc/poly.h"

void poly_rand(lssi_t* a, int coefs){
    int tmp_size = MAX_POLY_COEFS * LQBYTS;
    uint8_t tmp[tmp_size];
    int i, j, k;
    rand_bytes(tmp, coefs * LQBYTS);

    k = 0;
    for(i=0; i<coefs; ++i){
        a[i] = 0;
        for(j=0; j<LQBYTS; ++j){
            a[i] = (a[i]<<8) | tmp[k++];
        }
    }
}

void poly_mod_barrett_lssi(lssi_t* c, const lssi_t* a, const lssi_t q, int qbits, int coefs){
    int k = qbits>>1;
    lsdi_t mu = (1ll<<k)/q;
    lssi_t q_half = (q-1)/2;
    int i;

    for(i=0;i<coefs;++i){
        lsdi_t val = a[i];
        lsdi_t quotient = (val * mu)>>k;

        lssi_t r = (lssi_t)(val - quotient*q);
        
        if(r > 2*q || r < -2*q) r %= q;
        
        if(r > q_half) r -= q;
        if(r < -q_half) r+= q;
        c[i] = r;
    }
}

void poly_mod_barrett_lsdi(lssi_t* c, const lsdi_t* a, const lssi_t q, int qbits, int coefs){
    int k = qbits>>1;
    lsdi_t mu = (lsdi_t)(((__uint128_t)1 << k )/q);
    lssi_t q_half = (q-1)/2;
    int i;
    for(i=0;i<coefs;++i){
        lsdi_t val = a[i];
        lsdi_t quotient = (lsdi_t)(((__int128_t)val * mu) >> k);
        lsdi_t r = val - quotient*q;

        if(r > 2*q || r < -2*q) r %= q;
        
        if(r>q_half) r -= q;
        if(r<-q_half) r+=q;

        c[i] = (lssi_t)r;
    }
}

int poly_mod_equal(const lssi_t *a, const lssi_t *b, const lssi_t q, int qbits, int coefs){
    poly c, d;
    poly_mod_barrett_lssi(c, a, q, qbits, coefs);
    poly_mod_barrett_lssi(d, b, q, qbits, coefs);
    return POLY_EQ(c, d, coefs);
}

int poly_mod_is_zero(const lssi_t *a, const lssi_t q, int qbits, int coefs){
    poly c;
    poly_mod_barrett_lssi(c, a, q, qbits, coefs);
    return POLY_IS_ZERO(c, coefs);
}

void poly_neg(lssi_t *c, const lssi_t *a, const lssi_t q, int qbits, int coefs){
    poly tmp;
    int i;
    for(i=0; i<coefs; ++i){
        tmp[i] = -a[i];
    }
    poly_mod_barrett_lssi(c, tmp, q, qbits, coefs);
}

void poly_uniform_eta(poly *a, lp_ctx *ctx, int ETA) {
    int i, j;
    uint8_t buf[MAX_POLY_COEFS * LQBYTS];
    rand_bytes(buf, sizeof(buf));
    sample_eta((int32_t*)a, MAX_POLY_COEFS, buf, sizeof(buf), ctx, ETA);
}

void poly_uniform(poly *a, lp_ctx *ctx) {
    int i, j;
    uint8_t buf[MAX_POLY_COEFS * LQBYTS];
    rand_bytes(buf, sizeof(buf));
    sample_uniform((int32_t*)a, MAX_POLY_COEFS, buf, sizeof(buf), ctx);
}

void poly_uniform_CBD(poly *a, lp_ctx *ctx) {
    int i, j;
    uint8_t buf[MAX_POLY_COEFS * LQBYTS];
    rand_bytes(buf, sizeof(buf));
    sample_CBD((int32_t*)a, MAX_POLY_COEFS, buf, sizeof(buf), ctx);
}

unsigned int sample_uniform(int32_t *a,unsigned int len,const uint8_t *buf,unsigned int buflen,lp_ctx *ctx){
    const uint32_t Q = ctx->q;
    unsigned int ctr = 0;
    unsigned int pos = 0;
    uint32_t t;

    while (ctr < len && pos + 3 <= buflen) {
        uint32_t b0 = buf[pos], b1 = buf[pos+1], b2 = buf[pos+2];
        pos += 3;
        t = b0 | (b1 << 8) | (b2 << 16);
        t &= 0x7FFFFF;
        if (t < Q) {
            a[ctr++] = (int32_t)t;
        }
    }
    return ctr;
}

unsigned int sample_eta(int32_t *a,unsigned int len,const uint8_t *buf,unsigned int buflen,lp_ctx *ctx,int ETA){
    unsigned int ctr = 0;
    unsigned int pos = 0;
    uint32_t t0, t1;
    while (ctr < len && pos < buflen) {
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;
        if (ETA == 2) {
            if (t0 < 15) {
                t0 = t0 - (205 * t0 >> 10) * 5;
                a[ctr++] = ETA - (int32_t)t0;
            }
            if (t1 < 15 && ctr < len) {
                t1 = t1 - (205 * t1 >> 10) * 5;
                a[ctr++] = ETA - (int32_t)t1;
            }
        } else if (ETA == 4) {
            if (t0 < 9) {
                a[ctr++] = ETA - (int32_t)t0;
            }
            if (t1 < 9 && ctr < len) {
                a[ctr++] = ETA - (int32_t)t1;
            }
        } else {
            uint32_t range = 2 * ETA + 1;
            if (t0 < range) {
                a[ctr++] = (int32_t)t0 - ETA;
            }
            if (t1 < range && ctr < len) {
                a[ctr++] = (int32_t)t1 - ETA;
            }
        }
    }
    return ctr;
}

unsigned int sample_gamma1(int32_t *a,unsigned int len,const uint8_t *buf,unsigned int buflen,lp_ctx *ctx){
    unsigned int ctr = 0;
    unsigned int pos = 0;
    const int32_t GAMMA1_17 = (1 << 17);
    const int32_t GAMMA1_19 = (1 << 19);
    // Use GAMMA1_17 for Dilithium2, GAMMA1_19 for Dilithium3/5
    int32_t GAMMA1 = GAMMA1_17;
    
    if (GAMMA1 == GAMMA1_17) {
        while (ctr < len && pos + 8 <= buflen) {
            a[ctr++] = ((buf[pos] | ((uint32_t)buf[pos+1] << 8) | ((uint32_t)buf[pos+2] << 16)) & 0x3FFFF) - GAMMA1;
            if (ctr >= len) break;
            a[ctr++] = (((buf[pos+2] >> 2) | ((uint32_t)buf[pos+3] << 6) | ((uint32_t)buf[pos+4] << 14)) & 0x3FFFF) - GAMMA1;
            if (ctr >= len) break;
            a[ctr++] = (((buf[pos+4] >> 4) | ((uint32_t)buf[pos+5] << 4) | ((uint32_t)buf[pos+6] << 12)) & 0x3FFFF) - GAMMA1;
            if (ctr >= len) break;
            a[ctr++] = (((buf[pos+6] >> 6) | ((uint32_t)buf[pos+7] << 2) | ((uint32_t)buf[pos+8] << 10)) & 0x3FFFF) - GAMMA1;
            pos += 9;
        }
    } else if (GAMMA1 == GAMMA1_19) {
        while (ctr < len && pos + 4 <= buflen) {
            a[ctr++] = ((buf[pos] | ((uint32_t)buf[pos+1] << 8) | ((uint32_t)buf[pos+2] << 16)) & 0xFFFFF) - GAMMA1;
            if (ctr >= len) break;
            a[ctr++] = (((buf[pos+2] >> 4) | ((uint32_t)buf[pos+3] << 4) | ((uint32_t)buf[pos+4] << 12)) & 0xFFFFF) - GAMMA1;
            pos += 5;
        }
    }
    
    return ctr;
}

unsigned int sample_gamma2(int32_t *a,unsigned int len,const uint8_t *buf,unsigned int buflen,lp_ctx *ctx){
    unsigned int ctr = 0;
    unsigned int pos = 0;
    const uint32_t Q = ctx->q;
    const int32_t GAMMA2_88 = (Q - 1) / 88;
    const int32_t GAMMA2_32 = (Q - 1) / 32;
    int32_t GAMMA2 = GAMMA2_32;
    if (GAMMA2 == GAMMA2_88) {
        while (ctr < len && pos + 2 <= buflen) {
            a[ctr++] = buf[pos] & 0x3F;
            if (ctr >= len) break;
            a[ctr++] = ((buf[pos] >> 6) | ((buf[pos+1] & 0x0F) << 2)) & 0x3F;
            if (ctr >= len) break;
            a[ctr++] = ((buf[pos+1] >> 4) | ((buf[pos+2] & 0x03) << 4)) & 0x3F;
            if (ctr >= len) break;
            a[ctr++] = (buf[pos+2] >> 2) & 0x3F;
            pos += 3;
        }
    } else if (GAMMA2 == GAMMA2_32) {
        while (ctr < len && pos < buflen) {
            a[ctr++] = buf[pos] & 0x0F;
            if (ctr >= len) break;
            a[ctr++] = buf[pos] >> 4;
            pos++;
        }
    }
    
    return ctr;
}

unsigned int sample_CBD(int32_t *a,unsigned int len,const uint8_t *buf,unsigned int buflen,lp_ctx *ctx){
    unsigned int ctr = 0;
    unsigned int pos = 0;
    uint32_t t, d;
    int16_t a_val, b_val;
    
    const int eta = 2;
    
    if (eta == 2) {
        while (ctr < len && pos < buflen) {
            t = buf[pos++];
            d = t & 0x55555555;
            d += (t >> 1) & 0x55555555;
            for (int i = 0; i < 2 && ctr < len; i++) {
                a_val = (d >> (4*i)) & 0x03;
                b_val = (d >> (4*i + 2)) & 0x03;
                a[ctr++] = a_val - b_val;
            }
        }
    } else if (eta == 3) {
        while (ctr < len && pos + 2 <= buflen) {
            uint32_t t24 = buf[pos] | ((uint32_t)buf[pos+1] << 8) | ((uint32_t)buf[pos+2] << 16);
            pos += 3;
            
            for (int i = 0; i < 4 && ctr < len; i++) {
                uint32_t bits = (t24 >> (6*i)) & 0x3F;
                a_val = (bits & 0x07);
                a_val += ((bits >> 3) & 0x07);
                int a_sum = ((bits & 1) + ((bits >> 1) & 1) + ((bits >> 2) & 1));
                int b_sum = (((bits >> 3) & 1) + ((bits >> 4) & 1) + ((bits >> 5) & 1));
                a[ctr++] = a_sum - b_sum;
            }
        }
    }
    
    return ctr;
}