#include "poly.h"
#include "dili_sample.h"
#include "rand.h"

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
