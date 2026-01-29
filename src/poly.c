#include "../inc/types.h"
#include "../inc/rand.h"
#include "../inc/poly.h"
#include <stdio.h>

typedef lssi_t poly[MAX_POLY_COEFS];

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