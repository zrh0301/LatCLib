#include "pack.h"
#include "dili.h"
#include <string.h>

void dili_polyt1_pack(uint8_t *r, const poly p, int n, int d) {
    unsigned int i;
    for(i = 0; i < n; i++) {
        uint32_t t = p[i];
        r[i*2] = t & 0xFF;
        r[i*2+1] = (t >> 8) & 0xFF;
    }
}

void dili_polyt1_unpack(poly p, const uint8_t *a, int n, int d) {
    unsigned int i;
    for(i = 0; i < n; i++) {
        p[i] = a[i*2] | ((uint32_t)a[i*2+1] << 8);
        p[i] &= (1 << (d-1)) - 1;
    }
}

void dili_pack_pk(uint8_t *pk_bytes, const uint8_t *rho, const poly *t1, int k) {
    unsigned int i;
    
    for(i = 0; i < 32; i++)
        pk_bytes[i] = rho[i];
    pk_bytes += 32;
    
    for(i = 0; i < k; i++) {
        dili_polyt1_pack(pk_bytes, t1[i], DILI_N, DILI_D);
        pk_bytes += 2 * DILI_N;
    }
}

void dili_unpack_pk(uint8_t *rho, poly *t1, const uint8_t *pk_bytes, int k) {
    unsigned int i;
    
    for(i = 0; i < 32; i++)
        rho[i] = pk_bytes[i];
    pk_bytes += 32;
    
    for(i = 0; i < k; i++) {
        dili_polyt1_unpack(t1[i], pk_bytes, DILI_N, DILI_D);
        pk_bytes += 2 * DILI_N;
    }
}

void pack_z(uint8_t* r, const poly z){
    int i;
    uint32_t t;
    for(i=0; i<DILI_N; ++i){
        t = (uint32_t)(DILI_GAMMA1 - z[i]);
        r[3*i+0] = t & 0xff;
        r[3*i+1] = (t>>8) & 0xff;
        r[3*i+2] = (t>>16) & 0xff;
    }
}

void unpack_z(poly z, const uint8_t* r){
    int i;
    uint32_t t;
    for(i=0; i<DILI_N; ++i){
        t = r[3*i+0];
        t |= (uint32_t)(r[3*i+1]<<8);
        t |= (uint32_t)(r[3*i+2]<<16);
        z[i] = DILI_GAMMA1 - (lssi_t)t;
    }
}

void pack_w1(uint8_t* buf, const poly w1){
    int i;
    for(i = 0; i < DILI_N / 4; i++) {
        uint32_t t0 = w1[4*i+0] & 0x3F;
        uint32_t t1 = w1[4*i+1] & 0x3F;
        uint32_t t2 = w1[4*i+2] & 0x3F;
        uint32_t t3 = w1[4*i+3] & 0x3F;
        
        buf[3*i+0] = t0 | (t1 << 6);
        buf[3*i+1] = (t1 >> 2) | (t2 << 4);
        buf[3*i+2] = (t2 >> 4) | (t3 << 2);
    }
}

int pack_h(uint8_t* r, uint8_t h[][256], int k){
    int offset, i, j;
    offset = 0;
    memset(r, 0, DILI_OMEGA+k);
    for(i=0; i<k; ++i){
        for(j=0; j<DILI_N; ++j){
            if(h[i][j]){
                if(offset >= DILI_OMEGA) return 1; // error
                r[offset++] = (uint8_t)j;
            } 
        }
        r[DILI_OMEGA+i] = (uint8_t)offset;
    }
    return 0; // success
}

int unpack_h(uint8_t h[][256], const uint8_t* r, int k){
    int prev_limit, limit, idx, i, j;
    prev_limit = 0;
    for(i=0; i<k; ++i){
        for(j=0; j<DILI_N; ++j){
            h[i][j] = 0;
        }
        limit = r[DILI_OMEGA+i];
        if(limit < prev_limit || limit > DILI_OMEGA) {
            return -1;
        }
        for(j=prev_limit; j<limit; ++j){
            idx = r[j];
            if(idx >= DILI_N) {
                return -1;
            }
            h[i][idx] = 1;
        }
        prev_limit = limit;
    }
    return 0;  // 成功
}
