#ifndef POLY_H
#define POLY_H

#include "../inc/types.h"
#include "../inc/rand.h"
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

#define MAX_POLY_COEFS 256

#define MONTGOMERY_REDUCE(a, q, qc) ({ \
    int32_t _t = (int64_t)(int32_t)(a) * (qc); \
    _t = ((a) - (int64_t)_t * (q)) >> 32; \
    _t; \
})

#define FROM_MONTGOMERY(a, q, qc) MONTGOMERY_REDUCE((int64_t)(a), q, qc)

#define MOD_INV(a, q) mod_exp(a, (q)-2, q)
#define TO_MONT(a, q) (((lsdi_t)(a)<<32)%(q))

#define POLY_ZERO(a, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (a)[_i] = 0; \
} while(0)

#define POLY_COPY(dst, src, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (dst)[_i] = (src)[_i]; \
} while(0)

#define POLY_ADD(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i] + (b)[_i]; \
} while(0)

#define POLY_SUB(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i] - (b)[_i]; \
} while(0)

#define POLY_LSH_LOW(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i]<<(b); \
} while(0)

#define POLY_RSH_LOW(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i]>>(b); \
} while(0)

#define POLY_EQ(a, b, coefs) ({ \
    int _i; \
    int _result = 1; \
    for(_i=0; _i<(coefs); ++_i) { \
        if((a)[_i] != (b)[_i]) { \
            _result = 0; \
            break; \
        } \
    } \
    _result; \
})

#define POLY_IS_ZERO(a, coefs) ({ \
    int _i; \
    int _result = 1; \
    for(_i=0; _i<(coefs); ++_i) { \
        if((a)[_i] != 0) { \
            _result = 0; \
            break; \
        } \
    } \
    _result; \
})

#define POLY_MONT_INIT(qc, q) ({ \
    int _ret = -1; \
    if((q)&1) { \
        dig_t _x = (q); \
        dig_t _q1 = (dig_t)(q); \
        _x = (((_q1+2)&4)<<1)+_q1; \
        _x *= (2-_q1*_x); \
        _x *= (2-_q1*_x); \
        _x *= (2-_q1*_x); \
        *(qc) = (lssi_t)(~_x+1); \
        _ret = 0; \
    } \
    _ret; \
})

#define LSSI_MONT_MUL(x, y, q, qc) ({ \
    lssi_t _r = 0; \
    lsdi_t _a = 0; \
    _a = (lsdi_t)(x)*(y); \
    _r = (lssi_t)_a*(qc); \
    _r = (lssi_t)((_a+(lsdi_t)_r*(q))>>LQBITS); \
    if(_r >= (q)) _r -= (q); \
    _r; \
})

#define POLY_MONT_MUL(r, x, y, q, qc, coefs) do { \
    lsdi_t _a = 0; \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) { \
        _a = (lsdi_t)(x)[_i]*(y)[_i]; \
        (r)[_i] = (lssi_t)_a*(qc); \
        (r)[_i] = (lssi_t)((_a+(lsdi_t)(r)[_i]*(q))>>LQBITS); \
    } \
} while(0)

#define POLY_MONT_RDC(r, x, q, qc, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) { \
        (r)[_i] = (x)[_i]*(qc); \
        (r)[_i] = (lssi_t)(((lsdi_t)(x)[_i] + (lsdi_t)(r)[_i]*(q))>>LQBITS); \
    } \
} while(0)

typedef enum{
    NTT_FULL_MODE, // for Dilithium
    NTT_INCOMPLETE_MODE // for Kyber
} ntt_mode_t;

typedef struct {
    lssi_t q; // mod number
    lssi_t qc; // montgomery const qc = -q^{-1} mod q
    int n; // poly coef num
    int logn;
    ntt_mode_t mode;
    
    lssi_t n_inv;
    lssi_t zetas[MAX_POLY_COEFS]; // origin root table
    int rev_table[MAX_POLY_COEFS];
} lp_ctx;

void rand_bytes(uint8_t* output, int len);

void poly_rand(lssi_t* a, int coefs);
int poly_mod_equal(const lssi_t *a, const lssi_t *b, const lssi_t q, int qbits, int coefs);
int poly_mod_is_zero(const lssi_t *a, const lssi_t q, int qbits, int coefs);
void poly_neg(lssi_t *c, const lssi_t *a, const lssi_t q, int qbits, int coefs);

void poly_mod_barrett_lssi(lssi_t* c, const lssi_t* a, const lssi_t q, int qbits, int coefs);
void poly_mod_barrett_lsdi(lssi_t* c, const lsdi_t* a, const lssi_t q, int qbits, int coefs);

void init_ntt_dili(lp_ctx* Q);

void poly_ntt(lssi_t* a, lp_ctx* Q);
void poly_invntt(lssi_t* a, lp_ctx* Q);

void poly_base_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q);
void poly_scalar_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q);
#endif // POLY_H