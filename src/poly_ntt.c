#include "../inc/poly.h"
#include "../inc/dili_ntt_tables.h"
#include "../inc/kyber_ntt_tables.h"

typedef lssi_t poly[MAX_POLY_COEFS];

void init_ntt_dili(lp_ctx* Q){
    Q->n = 256;
    Q->q = 8380417;
    Q->mode = NTT_FULL_MODE;
    Q->layer = 8;
    Q->n_inv = DILI_NTT_N_INV;
    Q->qc = 58728449;
    Q->zetas = zetas_dili_256;
    Q->rev_table = DILI_NTT_REV_TABLE;
}

void init_ntt_kyber(lp_ctx* Q){
    Q->n = 256;
    Q->q = 3329;
    Q->mode = NTT_INCOMPLETE_MODE;
    Q->layer = 7;
    Q->n_inv = KYBER_NTT_N_INV;
    Q->qc = 2285;
    Q->zetas = zetas_kyber_256;
    Q->rev_table = NULL;
}

void poly_ntt(int32_t a[256], lp_ctx* Q) {
    int n, layer, min_len;
    int32_t q, qc;
    const int32_t* zetas;
    int len, start, j, k;
    int32_t zeta, t;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    layer = Q->layer;
    min_len = (n >> layer);
    k = 0;
    
    for(len = (n >> 1); len >= min_len; len >>= 1) {
        for(start = 0; start < n; start += (len << 1)) {
            zeta = zetas[++k];
            for(j = start; j < start + len; ++j) {
                t = MONTGOMERY_REDUCE((int64_t)a[j + len] * zeta, q, qc);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}

void poly_invntt(int32_t a[256], lp_ctx* Q) {
    int n, layer, step;
    int32_t q, qc, f;
    const int32_t* zetas;
    int start, len, j, k;
    int32_t t, zeta;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    layer = Q->layer;
    f = Q->n_inv;
    k = (n >> (8 - layer));
    
    for(len = (n >> layer); len < n; len <<= 1) {
        step = len << 1;
        for(start = 0; start < n; start += step) {
            zeta = -zetas[--k];
            for(j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = MONTGOMERY_REDUCE((int64_t)zeta * a[j + len], q, qc);
            }
        }
    }
    
    for(j = 0; j < n; ++j) {
        a[j] = MONTGOMERY_REDUCE((int64_t)f * a[j], q, qc);
    }
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
    int i, n;
    lssi_t q, qc;
    const lssi_t* zetas;
    
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    n = Q->n;
    
    for(i = 0; i < n; i+=2){
        BASEMUL(&r[i], &a[i], &b[i], zetas[128 + (i>>1)], q, qc);
    }
}

void poly_pointwise(lssi_t* c, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
    int32_t q, qc;
    int i;
    int n;
    q = Q->q;
    qc = Q->qc;
    n = Q->n;

    for(i=0;i<n;++i){
        c[i] = MONTGOMERY_REDUCE((lsdi_t)a[i]*b[i], q, qc);
    }
}