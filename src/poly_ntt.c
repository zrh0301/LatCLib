#include "../inc/poly.h"
#include "../inc/dili_ntt_tables.h"
#include "../inc/kyber_ntt_tables.h"

typedef lssi_t poly[MAX_POLY_COEFS];

void init_ntt_dili(lp_ctx* Q){
    Q->n = 256;
    Q->q = 8380417;
    Q->mode = NTT_FULL_MODE;
    Q->logn = 8;
    Q->n_inv = DILI_NTT_N_INV;
    Q->q = 8380417;
    Q->qc = 58728449;
    memcpy(Q->zetas, zetas_dili_256, sizeof(zetas_dili_256));
    memcpy(Q->rev_table, DILI_NTT_REV_TABLE, sizeof(DILI_NTT_REV_TABLE));
    Q->n_inv = DILI_NTT_N_INV;
}

void poly_ntt(int32_t a[256], lp_ctx* Q) {
    unsigned int len, start, j, k;
    const int32_t* zetas;
    int32_t zeta, t;
    k = 0;
    zetas = Q->zetas;
    len = Q->n >> 1;

    for(; len > 0; len >>= 1) {
        for(start = 0; start < Q->n; start = j + len) {
            zeta = zetas[++k];
            for(j = start; j < start + len; ++j) {
                t = MONTGOMERY_REDUCE((int64_t)zeta * a[j + len], Q->q, Q->qc);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}

void poly_invntt(int32_t a[256], lp_ctx* Q) {
    unsigned int start, len, j, k;
    int32_t t, zeta;
    int32_t f = 41978;
    const int32_t* zetas = Q->zetas;
    k = 256;
    f = 41978;

    for(len = 1; len < Q->n; len <<= 1) {
        for(start = 0; start < Q->n; start = j + len) {
            zeta = -zetas[--k];
            for(j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = MONTGOMERY_REDUCE((int64_t)zeta * a[j + len], Q->q, Q->qc);
            }
        }
    }
    
    for(j = 0; j < Q->n; ++j) {
        a[j] = MONTGOMERY_REDUCE((int64_t)f * a[j], Q->q, Q->qc);
    }
}


static lssi_t mod_exp(lssi_t base, lssi_t exp, lssi_t q){
    lssi_t res = 1;
    base %= q;
    uint32_t e = (uint32_t)exp;
    while(e > 0){
        if(e & 1) res = (lsdi_t)res*base%q;
        base = (lsdi_t)base*base%q;
        e >>= 1;
    }
    return res;
}

static void bit_reverse_array(lssi_t *a, int *rev, int n){
    int i;
    for (i=0;i<n;++i) {
        if (i<rev[i]){
            lssi_t temp = a[i];
            a[i] = a[rev[i]];
            a[rev[i]] = temp;
        }
    }
}


void poly_base_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
    int i;
    for(i = 0; i < Q->n; ++i) {
        r[i] = LSSI_MONT_MUL(a[i], b[i], Q->q, Q->qc);
    }
}

void poly_scalar_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
    int i;
    for(i=0;i<Q->n;++i){
        r[i] = LSSI_MONT_MUL(a[i], b[i], Q->q, Q->qc);
    }
}