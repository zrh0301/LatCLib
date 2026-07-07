#include "poly.h"

static lssi_t center_mod(lsdi_t a, lssi_t q) {
    a %= q;
    if(a < 0) a += q;
    if(a > q / 2) a -= q;
    return (lssi_t)a;
}

static int bit_reverse_width(int x, int width) {
    int r = 0;
    for(int i = 0; i < width; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

void ntt_generate_zetas(lssi_t zetas[MAX_POLY_COEFS], lssi_t q, lssi_t root,
                        lssi_t mont_r, int bitrev_bits, int used) {
    lssi_t powers[MAX_POLY_COEFS];
    lsdi_t root_power = 1;

    for(int i = 0; i < MAX_POLY_COEFS; i++) {
        powers[i] = center_mod(root_power * mont_r, q);
        root_power = (root_power * root) % q;
    }

    for(int i = 0; i < MAX_POLY_COEFS; i++) {
        zetas[i] = 0;
    }

    for(int i = 1; i < used; i++) {
        zetas[i] = powers[bit_reverse_width(i, bitrev_bits)];
    }
}

void ntt_init_dili(lp_ctx *Q) {
    Q->n = 256;
    Q->q = 8380417;
    Q->mode = NTT_FULL_MODE;
    Q->layer = 8;
    Q->f = 41978;
    Q->qc = 58728449;
    Q->zetas = Q->zetas_storage;
    Q->rev_table = NULL;
    ntt_generate_zetas(Q->zetas_storage, Q->q, 1753, 4193792, 8, 256);
}

void ntt_init_kyber(lp_ctx *Q) {
    Q->n = 256;
    Q->q = 3329;
    Q->mode = NTT_INCOMPLETE_MODE;
    Q->layer = 7;
    Q->f = 1104;
    Q->qc = 1806234369;
    Q->zetas = Q->zetas_storage;
    Q->rev_table = NULL;
    ntt_generate_zetas(Q->zetas_storage, Q->q, 17, 1353, 7, 128);
}
