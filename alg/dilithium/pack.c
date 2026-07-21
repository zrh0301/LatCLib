#include "pack.h"
#include "dili.h"
#include <string.h>

// t1: 10-bit 打包, 每 4 系数 -> 5 字节 (标准布局)
void dili_polyt1_pack(uint8_t *r, const poly p, int n, int d) {
    unsigned int i;
    (void)d;
    for(i = 0; i < (unsigned)n/4; i++) {
        r[5*i+0] = (uint8_t)(p[4*i+0] >> 0);
        r[5*i+1] = (uint8_t)((p[4*i+0] >> 8) | (p[4*i+1] << 2));
        r[5*i+2] = (uint8_t)((p[4*i+1] >> 6) | (p[4*i+2] << 4));
        r[5*i+3] = (uint8_t)((p[4*i+2] >> 4) | (p[4*i+3] << 6));
        r[5*i+4] = (uint8_t)(p[4*i+3] >> 2);
    }
}

void dili_polyt1_unpack(poly p, const uint8_t *a, int n, int d) {
    unsigned int i;
    (void)d;
    for(i = 0; i < (unsigned)n/4; i++) {
        p[4*i+0] = ((a[5*i+0] >> 0) | ((uint32_t)a[5*i+1] << 8)) & 0x3FF;
        p[4*i+1] = ((a[5*i+1] >> 2) | ((uint32_t)a[5*i+2] << 6)) & 0x3FF;
        p[4*i+2] = ((a[5*i+2] >> 4) | ((uint32_t)a[5*i+3] << 4)) & 0x3FF;
        p[4*i+3] = ((a[5*i+3] >> 6) | ((uint32_t)a[5*i+4] << 2)) & 0x3FF;
    }
}

void dili_pack_pk(uint8_t *pk_bytes, const uint8_t *rho, const poly *t1, int k) {
    unsigned int i;
    
    for(i = 0; i < 32; i++)
        pk_bytes[i] = rho[i];
    pk_bytes += 32;
    
    for(i = 0; i < k; i++) {
        dili_polyt1_pack(pk_bytes, t1[i], DILI_N, DILI_D);
        pk_bytes += DILI_POLYT1_PACKEDBYTES;
    }
}

void dili_unpack_pk(uint8_t *rho, poly *t1, const uint8_t *pk_bytes, int k) {
    unsigned int i;

    for(i = 0; i < 32; i++)
        rho[i] = pk_bytes[i];
    pk_bytes += 32;

    for(i = 0; i < k; i++) {
        dili_polyt1_unpack(t1[i], pk_bytes, DILI_N, DILI_D);
        pk_bytes += DILI_POLYT1_PACKEDBYTES;
    }
}

// z: 标准打包, 存 t = gamma1 - z
// gamma1=2^17 -> 18-bit, 每 4 系数 9 字节; gamma1=2^19 -> 20-bit, 每 2 系数 5 字节
void pack_z(uint8_t* r, const poly z, int gamma1){
    int i;
    uint32_t t[4];
    if(gamma1 == (1 << 17)) {
        for(i = 0; i < DILI_N/4; i++){
            t[0] = (uint32_t)(gamma1 - z[4*i+0]);
            t[1] = (uint32_t)(gamma1 - z[4*i+1]);
            t[2] = (uint32_t)(gamma1 - z[4*i+2]);
            t[3] = (uint32_t)(gamma1 - z[4*i+3]);
            r[9*i+0] = (uint8_t)(t[0]);
            r[9*i+1] = (uint8_t)(t[0] >> 8);
            r[9*i+2] = (uint8_t)((t[0] >> 16) | (t[1] << 2));
            r[9*i+3] = (uint8_t)(t[1] >> 6);
            r[9*i+4] = (uint8_t)((t[1] >> 14) | (t[2] << 4));
            r[9*i+5] = (uint8_t)(t[2] >> 4);
            r[9*i+6] = (uint8_t)((t[2] >> 12) | (t[3] << 6));
            r[9*i+7] = (uint8_t)(t[3] >> 2);
            r[9*i+8] = (uint8_t)(t[3] >> 10);
        }
    } else {  // gamma1 = 2^19
        for(i = 0; i < DILI_N/2; i++){
            t[0] = (uint32_t)(gamma1 - z[2*i+0]);
            t[1] = (uint32_t)(gamma1 - z[2*i+1]);
            r[5*i+0] = (uint8_t)(t[0]);
            r[5*i+1] = (uint8_t)(t[0] >> 8);
            r[5*i+2] = (uint8_t)((t[0] >> 16) | (t[1] << 4));
            r[5*i+3] = (uint8_t)(t[1] >> 4);
            r[5*i+4] = (uint8_t)(t[1] >> 12);
        }
    }
}

void unpack_z(poly z, const uint8_t* r, int gamma1){
    int i;
    if(gamma1 == (1 << 17)) {
        for(i = 0; i < DILI_N/4; i++){
            z[4*i+0] = (r[9*i+0] | ((uint32_t)r[9*i+1] << 8) | ((uint32_t)r[9*i+2] << 16)) & 0x3FFFF;
            z[4*i+1] = ((r[9*i+2] >> 2) | ((uint32_t)r[9*i+3] << 6) | ((uint32_t)r[9*i+4] << 14)) & 0x3FFFF;
            z[4*i+2] = ((r[9*i+4] >> 4) | ((uint32_t)r[9*i+5] << 4) | ((uint32_t)r[9*i+6] << 12)) & 0x3FFFF;
            z[4*i+3] = ((r[9*i+6] >> 6) | ((uint32_t)r[9*i+7] << 2) | ((uint32_t)r[9*i+8] << 10)) & 0x3FFFF;
            z[4*i+0] = gamma1 - z[4*i+0];
            z[4*i+1] = gamma1 - z[4*i+1];
            z[4*i+2] = gamma1 - z[4*i+2];
            z[4*i+3] = gamma1 - z[4*i+3];
        }
    } else {  // gamma1 = 2^19
        for(i = 0; i < DILI_N/2; i++){
            z[2*i+0] = (r[5*i+0] | ((uint32_t)r[5*i+1] << 8) | ((uint32_t)r[5*i+2] << 16)) & 0xFFFFF;
            z[2*i+1] = ((r[5*i+2] >> 4) | ((uint32_t)r[5*i+3] << 4) | ((uint32_t)r[5*i+4] << 12)) & 0xFFFFF;
            z[2*i+0] = gamma1 - z[2*i+0];
            z[2*i+1] = gamma1 - z[2*i+1];
        }
    }
}

void pack_w1(uint8_t* buf, const poly w1, int gamma2){
    int i;
    
    if(gamma2 == ((DILI_Q-1)/88)) {  // Dilithium2
        for(i = 0; i < DILI_N / 4; i++) {
            uint32_t t0 = w1[4*i+0] & 0x3F;  // 6 bits
            uint32_t t1 = w1[4*i+1] & 0x3F;
            uint32_t t2 = w1[4*i+2] & 0x3F;
            uint32_t t3 = w1[4*i+3] & 0x3F;
            
            buf[3*i+0] = t0 | (t1 << 6);
            buf[3*i+1] = (t1 >> 2) | (t2 << 4);
            buf[3*i+2] = (t2 >> 4) | (t3 << 2);
        }
    } else {  // Dilithium3/5
        for(i = 0; i < DILI_N / 2; i++) {
            uint32_t t0 = w1[2*i+0] & 0x0F;  // 4 bits
            uint32_t t1 = w1[2*i+1] & 0x0F;
            buf[i] = t0 | (t1 << 4);
        }
    }
}

int pack_h(uint8_t* r, uint8_t h[][256], int k, int omega){
    int offset, i, j;
    offset = 0;
    memset(r, 0, omega+k);
    for(i=0; i<k; ++i){
        for(j=0; j<DILI_N; ++j){
            if(h[i][j]){
                if(offset >= omega) return 1;
                r[offset++] = (uint8_t)j;
            } 
        }
        r[omega+i] = (uint8_t)offset;
    }
    return 0;
}

int unpack_h(uint8_t h[][256], const uint8_t* r, int k, int omega){
    int prev_limit, limit, idx, i, j;
    prev_limit = 0;
    for(i=0; i<k; ++i){
        for(j=0; j<DILI_N; ++j){
            h[i][j] = 0;
        }
        limit = r[omega+i];
        if(limit < prev_limit || limit > omega) {
            return -1;
        }
        for(j=prev_limit; j<limit; ++j){
            idx = r[j];
            if(idx >= DILI_N) {
                return -1;
            }
            // 同一多项式内索引必须严格递增, 否则拒绝 (防止签名可锻性)
            if(j > prev_limit && idx <= r[j-1]) {
                return -1;
            }
            h[i][idx] = 1;
        }
        prev_limit = limit;
    }
    // 尾部填充字节必须全 0
    for(j = prev_limit; j < omega; ++j){
        if(r[j] != 0) {
            return -1;
        }
    }
    return 0;
}
