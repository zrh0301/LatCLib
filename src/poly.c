#include "../inc/types.h"
#include "../inc/rand.h"
#include "../inc/poly.h"
#include <stdio.h>

typedef lssi_t poly[MAX_POLY_COEFS];

void poly_zero(lssi_t* a, int coefs){
    for(int i=0;i<coefs;++i) a[i] = 0;
}

void poly_rand(lssi_t* a, int coefs){
    int tmp_size = MAX_POLY_COEFS * LQBYTS;
    uint8_t tmp [tmp_size];
    rand_bytes(tmp, coefs * LQBYTS);

    int k = 0;
    for(int i=0;i<coefs;++i){
        a[i] = 0;
        for(int j=0;j<LQBYTS;++j){
            a[i] = (a[i]<<8) | tmp[k++];
        }
    }
}

void poly_copy(lssi_t* dst, const lssi_t* src, int coefs){
    for(int i=0;i<coefs;++i){
        dst[i] = src[i];
    }
}

//cmp 0->noe equal; 1->equal
int poly_eq(const lssi_t* a, const lssi_t* b, int coefs){
    for(int i=0;i<coefs;++i){
        if(a[i] != b[i]) return 0;
    }
    return 1;
}

int poly_is_zero(const lssi_t* a, int coefs){
    for(int i=0;i<coefs;++i){
        if(a[i] != 0) return 0;
    }
    return 1;
}

void poly_add(lssi_t* c, const lssi_t* a, const lssi_t* b, int coefs){
    for(int i=0;i<coefs;++i) c[i] = a[i] + b[i];
}

void poly_sub(lssi_t* c, const lssi_t* a, const lssi_t* b, int coefs){
    for(int i=0;i<coefs;++i) c[i] = a[i] - b[i];
}

void poly_lsh_low(lssi_t* c, const lssi_t* a, int b, int coefs){
    for(int i=0;i<coefs;++i) c[i] = a[i]<<b;
}

void poly_rsh_low(lssi_t* c, const lssi_t* a, int b, int coefs){
    for(int i=0;i<coefs;++i) c[i] = a[i]>>b; 
}

void poly_mod_barrett_lssi(lssi_t* c, const lssi_t* a, const lssi_t q, int qbits, int coefs){
    int k = qbits>>1;
    lsdi_t mu = (1ll<<k)/q;
    lssi_t q_half = (q-1)/2;

    for(int i=0;i<coefs;++i){
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
    for(int i=0;i<coefs;++i){
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
    return poly_eq(c, d, coefs);
}

int poly_mod_is_zero(const lssi_t *a, const lssi_t q, int qbits, int coefs){
    poly c;
    poly_mod_barrett_lssi(c, a, q, qbits, coefs);
    return poly_is_zero(c, coefs);
}

void poly_neg(lssi_t *c, const lssi_t *a, const lssi_t q, int qbits, int coefs){
    poly tmp;
    for(int i=0;i<coefs;++i){
        tmp[i] = -a[i];
    }
    poly_mod_barrett_lssi(c, tmp, q, qbits, coefs);
}

int poly_mont_init(lssi_t* qc, const lssi_t q){
    if(!(q&1)) return -1;
    
    dig_t x = q;
    dig_t q1 = (dig_t)q;
    x = (((q1+2)&4)<<1)+q1;
    x *= (2-q1*x);
    x *= (2-q1*x);
    x *= (2-q1*x);

    *qc = (lssi_t)(~x+1);
    return 0;
}

lssi_t lssi_mont_mul(const lssi_t x, const lssi_t y, const lssi_t q, lssi_t qc){
    lssi_t r = 0;
    lsdi_t a = 0;
    a = (lsdi_t)x*y;
    r = (lssi_t)a*qc;
    r = (lssi_t)((a+(lsdi_t)r*q)>>LQBITS);
    if(r >= q) return r-q;
    return r;
}

void poly_mont_mul(lssi_t *r, const lssi_t *x, const lssi_t *y, const lssi_t q, lssi_t qc, int coefs){
	lsdi_t a = 0;
	for(int i=0;i<coefs;++i){
		a = (lsdi_t)x[i]*y[i];
		r[i] = (lssi_t)a*qc;
		r[i] = (lssi_t)((a+(lsdi_t)r[i]*q)>>LQBITS);
	}
}

void poly_mont_rdc(lssi_t* r, const lssi_t* x, const lssi_t q, lssi_t qc, int coefs){
    for(int i=0;i<coefs;++i){
        r[i] = x[i]*qc;
        r[i] = (lssi_t)(((lsdi_t)x[i] + (lsdi_t)r[i]*q )>>LQBITS);
    }
}