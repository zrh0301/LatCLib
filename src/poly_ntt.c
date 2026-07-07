#include "../inc/poly.h"

typedef lssi_t poly[MAX_POLY_COEFS];

void ntt_vector(poly* vec, int count, lp_ctx* Q) {
    for(int i = 0; i < count; i++) {
        poly_ntt(vec[i], Q);
    }
}

void ntt_matrix(poly* mat, int rows, int cols, lp_ctx* Q) {
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            poly_ntt(mat[i * cols + j], Q);
        }
    }
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
    f = Q->f * (1 << (8 - layer));
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

void poly_base_pointwise(lssi_t * r, const lssi_t * a, const lssi_t * b, lp_ctx * Q){
    unsigned int i, n, x, y;
    const int32_t* zetas;
    int32_t q, qc, zeta;
    int32_t t0, t1, t2;
    
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    n = Q->n / 4;

    for(i=0; i<n; i++) {
        x = 4*i;
        y = x + 2;
        zeta = zetas[64+i];
        
        t0 = MONTGOMERY_REDUCE((lsdi_t)a[x] * b[x], q, qc);
        t1 = MONTGOMERY_REDUCE((lsdi_t)a[x+1] * b[x+1], q, qc);
        t2 = MONTGOMERY_REDUCE((lsdi_t)(a[x] + a[x+1]) * (b[x] + b[x+1]), q, qc);
        r[x] = MONTGOMERY_REDUCE((lsdi_t)t1 * zeta, q, qc) + t0;
        r[x+1] = t2 - t0 - t1;
        
        t0 = MONTGOMERY_REDUCE((lsdi_t)a[y] * b[y], q, qc);
        t1 = MONTGOMERY_REDUCE((lsdi_t)a[y+1] * b[y+1], q, qc);
        t2 = MONTGOMERY_REDUCE((lsdi_t)(a[y] + a[y+1]) * (b[y] + b[y+1]), q, qc);
        r[y] = MONTGOMERY_REDUCE((lsdi_t)t1 * (-zeta), q, qc) + t0;
        r[y+1] = t2 - t0 - t1;
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