#include "../inc/poly.h"
#include "../inc/dili_ntt_tables.h"
#include "../inc/kyber_ntt_tables.h"

typedef lssi_t poly[MAX_POLY_COEFS];

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

static lssi_t mod_inv(lssi_t a, lssi_t q){
    return mod_exp(a, q-2, q);
}

static lssi_t to_mont(lssi_t a, lssi_t q){
    return ((lsdi_t)a<<32)%q;
}

static void bit_reverse_array(lssi_t *a, int *rev, int n){
    for (int i=0;i<n;++i) {
        if (i<rev[i]){
            lssi_t temp = a[i];
            a[i] = a[rev[i]];
            a[rev[i]] = temp;
        }
    }
}

void init_ntt_kyber(lp_ctx* Q){
    Q->n = 256;
    Q->q = 3329;
    Q->mode = NTT_INCOMPLETE_MODE;
    Q->logn = 8;
    
    poly_mont_init(&Q->qc, Q->q);
    
    memcpy(Q->w, KYBER_NTT_W, sizeof(KYBER_NTT_W));
    memcpy(Q->w_inv, KYBER_NTT_W_INV, sizeof(KYBER_NTT_W_INV));
    memcpy(Q->rev_table, KYBER_NTT_REV_TABLE, sizeof(KYBER_NTT_REV_TABLE));
    Q->n_inv = KYBER_NTT_N_INV;
}

void init_ntt_dili(lp_ctx* Q){
    Q->n = 256;
    Q->q = 8380417;
    Q->mode = NTT_FULL_MODE;
    Q->logn = 8;
    
    poly_mont_init(&Q->qc, Q->q);
    
    memcpy(Q->w, DILI_NTT_W, sizeof(DILI_NTT_W));
    memcpy(Q->w_inv, DILI_NTT_W_INV, sizeof(DILI_NTT_W_INV));
    memcpy(Q->rev_table, DILI_NTT_REV_TABLE, sizeof(DILI_NTT_REV_TABLE));
    Q->n_inv = DILI_NTT_N_INV;
}

void init_ntt(lp_ctx* Q, int coefs, lssi_t root, ntt_mode_t mode){
    if(coefs > MAX_POLY_COEFS) return; 
    Q->n = coefs;
    Q->q = (mode == NTT_FULL_MODE) ? 8380417 : 3329;
    Q->mode = mode;
    
    
    Q->logn = 0;
    while(coefs > (1<<Q->logn)) ++Q->logn;

    poly_mont_init(&(Q->qc), Q->q);

    
    int k;
    for(int i=0; i<coefs; ++i){
        k = 0;
        for(int j=0; j<Q->logn; ++j){
            if(i & (1<<j)) k |= (1<<(Q->logn-1-j));
        }
        Q->rev_table[i] = k;
    }

    
    lssi_t inv_root = mod_inv(root, Q->q);

    if (mode == NTT_FULL_MODE){
        lssi_t w_cur = 1;
        lssi_t w_inv_cur = 1;
        for(int i=0;i<coefs;++i){
            Q->w[i] = to_mont(w_cur, Q->q);
            Q->w_inv[i] = to_mont(w_inv_cur, Q->q);
            w_cur = (lsdi_t)w_cur * root % Q->q;
            w_inv_cur = (lsdi_t)w_inv_cur * inv_root % Q->q;
        }
    } 
    else{
        for(int i=0;i<coefs;++i){
            lssi_t w_pow = mod_exp(root, Q->rev_table[i], Q->q);
            lssi_t w_inv_pow = mod_exp(inv_root, Q->rev_table[i], Q->q);
            
            Q->w[i] = to_mont(w_pow, Q->q);
            Q->w_inv[i] = to_mont(w_inv_pow, Q->q);
        }
    }

    int scale_div = (mode == NTT_FULL_MODE) ? coefs : (coefs >> 1);
    lssi_t n_inv = mod_inv(scale_div, Q->q);
    Q->n_inv = to_mont(n_inv, Q->q);
}

void poly_ntt(lssi_t* a, lp_ctx* Q){
    lssi_t q = Q->q;
    lssi_t qc = Q->qc;
    int coefs = Q->n;
    if (Q->mode == NTT_FULL_MODE){
        bit_reverse_array(a, Q->rev_table, coefs);
        for(int len=1;len<coefs;len<<=1){
            int step = (coefs/len) >> 1;
            for(int i=0;i<coefs;i+=(len<<1)){
                for(int j=0;j<len;++j){
                    lssi_t zeta = Q->w[j*step];
                    lssi_t u = a[i+j];
                    lssi_t v = lssi_mont_mul(a[i+j+len], zeta, q, qc);
                    
                    lssi_t res_add = u + v;
                    lssi_t res_sub = u - v;
                    lssi_t t_add = res_add - q;
                    res_add = t_add + (q & (t_add >> 31));
                    res_sub += q & (res_sub >> 31);

                    a[i+j] = res_add;
                    a[i+j+len] = res_sub;
                }
            }
        }
    } 
    else{
        int len, start, j, k;
        k = 1; 
        for(int len = 128;len>=2;len>>=1){
            for(int start=0;start<256;start+=2*len){
                lssi_t zeta = Q->w[k++];
                
                for(j=start;j<start+len;++j) {
                    lssi_t t = lssi_mont_mul(zeta, a[j + len], q, qc);
                    lssi_t u = a[j];
                    
                    lssi_t res_add = u + t;
                    lssi_t res_sub = u - t;
                    lssi_t t_add = res_add - q;
                    res_add = t_add + (q & (t_add >> 31));
                    res_sub += q & (res_sub >> 31);
                    
                    a[j] = res_add;
                    a[j + len] = res_sub;
                }
            }
        }
    }
}

void poly_invntt(lssi_t* a, lp_ctx* Q){
    lssi_t q = Q->q;
    lssi_t qc = Q->qc;
    int coefs = Q->n;

    if(Q->mode == NTT_FULL_MODE){
        bit_reverse_array(a, Q->rev_table, coefs);
        for(int len=1;len<coefs;len<<=1){
            int step = (coefs/len)>>1;
            for(int i=0;i<coefs;i+=2*len){
                int k = 0;
                for(int j=0;j<len;++j){
                    lssi_t zeta = Q->w_inv[k*step];
                    lssi_t u = a[i+j];
                    lssi_t v = lssi_mont_mul(a[i+j+len], zeta, q, qc);

                    lssi_t res_add = u+v;
                    lssi_t res_sub = u-v;
                    if(res_add >= q) res_add -= q;
                    if(res_sub < 0) res_sub += q;
                    
                    a[i+j] = res_add;
                    a[i+j+len] = res_sub;
                    ++k;
                }
            }
        }
        for(int i=0;i<coefs;++i){
            a[i] = lssi_mont_mul(a[i], Q->n_inv, q, qc);
        }
    }
    else{
        int k = 127;
        for(int len=2;len<=128;len<<=1){
            
            int start_k = k - (256 / (2 * len)) + 1;
            int current_k = start_k;

            for(int start=0;start<256;start+=2*len){
                lssi_t zeta = Q->w_inv[current_k++];
                
                for(int j=start;j<start+len;++j) {
                    lssi_t u = a[j];
                    lssi_t v = a[j + len];
                    
                    lssi_t res_add = u + v;
                    lssi_t res_sub = u - v;
                    
                    if(res_add >= q) res_add -= q;
                    if(res_sub < 0)  res_sub += q;
                    
                    a[j] = res_add;
                    a[j + len] = lssi_mont_mul(res_sub, zeta, q, qc);
                }
            }
            k = start_k - 1;
        }
        for(int i=0; i<coefs; ++i){
            a[i] = lssi_mont_mul(a[i], Q->n_inv, q, qc);
        }
    }
}

void poly_base_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
    lssi_t q = Q->q;
    lssi_t qc = Q->qc;
    
    if (Q->mode == NTT_FULL_MODE){
        for(int i=0; i<Q->n; ++i) {
            r[i] = lssi_mont_mul(a[i], b[i], q, qc);
        }
    } 
    else{
        
        for(int i=0;i<(Q->n>>1);++i){
            lssi_t zeta = Q->w[64 + i]; 
            
            lssi_t a0 = a[2*i];
            lssi_t a1 = a[2*i+1];
            lssi_t b0 = b[2*i];
            lssi_t b1 = b[2*i+1];
            
            lssi_t val0 = lssi_mont_mul(a0, b0, q, qc);
            lssi_t val1 = lssi_mont_mul(a1, b1, q, qc);
            val1 = lssi_mont_mul(val1, zeta, q, qc);
            lssi_t res0 = val0 + val1;

            lssi_t val2 = lssi_mont_mul(a0, b1, q, qc);
            lssi_t val3 = lssi_mont_mul(a1, b0, q, qc);
            lssi_t res1 = val2 + val3;

            if(res0 >= q) res0 -= q;
            if(res1 >= q) res1 -= q;

            r[2*i] = res0;
            r[2*i+1] = res1;
        }
    }
}

void poly_scalar_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
    for(int i=0;i<Q->n;++i){
        r[i] = lssi_mont_mul(a[i], b[i], Q->q, Q->qc);
    }
}