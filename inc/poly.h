#ifndef POLY_H
#define POLY_H

#include "../inc/types.h"
#include "../inc/rand.h"
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

#define MAX_POLY_COEFS 256

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
    lssi_t w[MAX_POLY_COEFS]; // origin root table
    lssi_t w_inv[MAX_POLY_COEFS];
    int rev_table[MAX_POLY_COEFS];
} lp_ctx;

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

static lssi_t mod_exp(lssi_t base, lssi_t exp, lssi_t q);
static lssi_t mod_inv(lssi_t a, lssi_t q);
static lssi_t to_mont(lssi_t a, lssi_t q);
static void bit_reverse(lssi_t *a, int *rev, int n);
void init_ntt(lp_ctx* Q, int coefs, lssi_t root, ntt_mode_t mode);
void poly_ntt(lssi_t* a, lp_ctx* Q);
void poly_invntt(lssi_t* a, lp_ctx* Q);
void poly_base_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q);

#endif // POLY_H