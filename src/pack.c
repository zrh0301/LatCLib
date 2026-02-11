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

void pack_z(uint8_t* r, const poly z, int gamma1){
    int i;
    uint32_t t;
    for(i=0; i<DILI_N; ++i){
        t = (uint32_t)(gamma1 - z[i]);
        r[3*i+0] = t & 0xff;
        r[3*i+1] = (t>>8) & 0xff;
        r[3*i+2] = (t>>16) & 0xff;
    }
}

void unpack_z(poly z, const uint8_t* r, int gamma1){
    int i;
    uint32_t t;
    for(i=0; i<DILI_N; ++i){
        t = r[3*i+0];
        t |= (uint32_t)(r[3*i+1]<<8);
        t |= (uint32_t)(r[3*i+2]<<16);
        z[i] = gamma1 - (lssi_t)t;
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
            h[i][idx] = 1;
        }
        prev_limit = limit;
    }
    return 0;
}
