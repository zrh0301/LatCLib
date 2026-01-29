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
    Q->qc = 58728449;
    Q->zetas = zetas_dili_256;
    Q->rev_table = DILI_NTT_REV_TABLE;
}

void poly_ntt(int32_t a[256], lp_ctx* Q) {
    unsigned int n, len, start, j, k, jlen;
    int32_t q, qc;
    const int32_t* zetas;
    int32_t zeta, t;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    k = 0;
    len = n >> 1;

    for(; len > 0; len >>= 1) {
        for(start = 0; start < n; start = j + len) {
            zeta = zetas[++k];
            for(j = start; j < start + len; ++j) {
                jlen = j + len;
                t = MONTGOMERY_REDUCE((int64_t)zeta * a[jlen], q, qc);
                a[jlen] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}

void poly_invntt(int32_t a[256], lp_ctx* Q) {
    unsigned int n, start, len, j, k, jlen;
    int32_t q, qc, f;
    const int32_t* zetas;
    int32_t t, zeta;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    f = 41978;
    k = 256;

    for(len = 1; len < n; len <<= 1) {
        for(start = 0; start < n; start = j + len) {
            zeta = -zetas[--k];
            for(j = start; j < start + len; ++j) {
                jlen = j + len;
                t = a[j];
                a[j] = t + a[jlen];
                a[jlen] = t - a[jlen];
                a[jlen] = MONTGOMERY_REDUCE((int64_t)zeta * a[jlen], q, qc);
            }
        }
    }
    
    for(j = 0; j < n; ++j) {
        a[j] = MONTGOMERY_REDUCE((int64_t)f * a[j], q, qc);
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
    int n, i;
    int32_t q, qc;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    
    for(i = 0; i < n; ++i) {
        r[i] = LSSI_MONT_MUL(a[i], b[i], q, qc);
    }
}

void poly_pointwise(lssi_t* c, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
  int32_t q, qc;
  int i;
  
  q = Q->q;
  qc = Q->qc;
  
  for(i=0;i<256;++i){
    c[i] = MONTGOMERY_REDUCE((lsdi_t)a[i]*b[i], q, qc);
  }
}