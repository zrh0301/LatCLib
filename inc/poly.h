#ifndef POLY_H
#define POLY_H

#include "../inc/types.h"
#include <sys/random.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

#define MAX_POLY_COEFS 256

typedef struct {
    lssi_t q; // mod number
    lssi_t qc; // montgomery const qc = -q^{-1} mod q
    int n; // poly coef num
    int logn;
    lssi_t n_inv;
    lssi_t w[MAX_POLY_COEFS]; // origin root table
    lssi_t w_inv[MAX_POLY_COEFS];
    int rev_table[MAX_POLY_COEFS];
} lp_ctx;

#define MONTGOMERY_REDUCE(a, q, qc) ({ \
    int32_t _t = (int64_t)(int32_t)(a) * (qc); \
    _t = ((a) - (int64_t)_t * (q)) >> 32; \
    _t; \
})

void rand_bytes(uint8_t* output, int len);
void poly_zero(lssi_t* a, int coefs);
void poly_rand(lssi_t* a, int coefs);

void poly_copy(lssi_t* dst, const lssi_t* src, int coefs);
int poly_eq(const lssi_t* a, const lssi_t* b, int coefs);
int poly_is_zero(const lssi_t* a, int coefs);

void poly_add(lssi_t* c, const lssi_t* a, const lssi_t* b, int coefs);
void poly_sub(lssi_t* c, const lssi_t* a, const lssi_t* b, int coefs);
void poly_lsh_low(lssi_t* c, const lssi_t* a, int b, int coefs);
void poly_rsh_low(lssi_t* c, const lssi_t* a, int b, int coefs);

void poly_mod_barrett_lssi(lssi_t* c, const lssi_t* a, const lssi_t q, int qbits, int coefs);
void poly_mod_barrett_lsdi(lssi_t* c, const lsdi_t* a, const lssi_t q, int qbits, int coefs);
int poly_mod_equal(const lssi_t *a, const lssi_t *b, const lssi_t q, int qbits, int coefs);
int poly_mod_is_zero(const lssi_t *a, const lssi_t q, int qbits, int coefs);
void poly_neg(lssi_t *c, const lssi_t *a, const lssi_t q, int qbits, int coefs);

int poly_mont_init(lssi_t* qc, const lssi_t q);
lssi_t lssi_mont_mul(const lssi_t x, const lssi_t y, const lssi_t q, lssi_t qc);
void poly_mont_mul(lssi_t *r, const lssi_t *x, const lssi_t *y, const lssi_t q, lssi_t qc, int coefs);
void poly_mont_rdc(lssi_t* r, const lssi_t* x, const lssi_t q, lssi_t qc, int coefs);

void poly_ntt(lssi_t* a,lp_ctx* Q);
void poly_intt(lssi_t* a, lp_ctx* Q);

#endif