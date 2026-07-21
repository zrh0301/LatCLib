#include "kyber.h"
#include "kyber_pack.h"
#include "kyber_sample.h"
#include "polyvec.h"
#include <string.h>

static void kyber_matrix_sample_ntt(poly out, const uint8_t seed[34],
                                    lp_ctx *ctx, void *arg) {
    (void)arg;
    kyber_sample_ntt(out, seed, ctx);
}

static void kyber_matrix_expand(poly *mat, int k,
                                const uint8_t rho[KYBER_SYMBYTES],
                                lp_ctx *ctx) {
    polyvec_matrix_expand_with_sampler(mat, k, k, rho, ctx,
                                       kyber_matrix_sample_ntt, NULL);
}

void kyber_init_param(kyber_param *param, int mode) {
    memset(param, 0, sizeof(*param));
    param->mode = mode;

    if(mode == KYBER_MODE_2) {
        param->k = 2;
        param->eta1 = 3;
        param->eta2 = 2;
        param->du = 10;
        param->dv = 4;
    } else if(mode == KYBER_MODE_3) {
        param->k = 3;
        param->eta1 = 2;
        param->eta2 = 2;
        param->du = 10;
        param->dv = 4;
    } else if(mode == KYBER_MODE_5) {
        param->k = 4;
        param->eta1 = 2;
        param->eta2 = 2;
        param->du = 11;
        param->dv = 5;
    } else {
        return;
    }

    param->ek_bytes = param->k * KYBER_POLYBYTES + KYBER_SYMBYTES;
    param->dk_pke_bytes = param->k * KYBER_POLYBYTES;
    param->ct_bytes = param->k * KYBER_N * param->du / 8 + KYBER_N * param->dv / 8;
}

const char *kyber_mode_name(int mode) {
    if(mode == KYBER_MODE_2) {
        return "Kyber2";
    }
    if(mode == KYBER_MODE_3) {
        return "Kyber3";
    }
    if(mode == KYBER_MODE_5) {
        return "Kyber5";
    }
    return "Kyber-unknown";
}

int kyber_pke_keygen(uint8_t *ek_pke, uint8_t *dk_pke,
                     const uint8_t d[KYBER_SYMBYTES], kyber_param param) {
    lp_ctx ctx;
    uint8_t seedbuf[KYBER_SYMBYTES + 1];
    uint8_t gbuf[64];
    uint8_t rho[KYBER_SYMBYTES];
    uint8_t sigma[KYBER_SYMBYTES];
    poly A[KYBER_MAX_K * KYBER_MAX_K];
    poly s[KYBER_MAX_K];
    poly e[KYBER_MAX_K];
    poly s_hat[KYBER_MAX_K];
    poly e_hat[KYBER_MAX_K];
    poly t_hat[KYBER_MAX_K];
    uint8_t nonce = 0;

    if(!kyber_param_is_valid(param)) {
        return -1;
    }

    INIT_NTT_KYBER(&ctx);
    memcpy(seedbuf, d, KYBER_SYMBYTES);
    seedbuf[KYBER_SYMBYTES] = (uint8_t)param.k;
    kyber_g(gbuf, seedbuf, sizeof(seedbuf));
    memcpy(rho, gbuf, KYBER_SYMBYTES);
    memcpy(sigma, gbuf + KYBER_SYMBYTES, KYBER_SYMBYTES);

    kyber_matrix_expand(A, param.k, rho, &ctx);

    for(int i = 0; i < param.k; i++) {
        kyber_sample_noise(s[i], sigma, nonce++, param.eta1);
    }
    for(int i = 0; i < param.k; i++) {
        kyber_sample_noise(e[i], sigma, nonce++, param.eta1);
    }

    polyvec_copy(s_hat, s, param.k);
    polyvec_copy(e_hat, e, param.k);
    polyvec_ntt(s_hat, param.k, &ctx);
    polyvec_ntt(e_hat, param.k, &ctx);

    polyvec_matrix_pointwise(t_hat, A, param.k, param.k, s_hat, &ctx, POLYVEC_MUL_BASE);
    polyvec_tomont(t_hat, param.k, &ctx);
    polyvec_add(t_hat, t_hat, e_hat, param.k);
    polyvec_reduce_mod(t_hat, param.k, KYBER_Q);

    kyber_pack_ek(ek_pke, t_hat, rho, param);
    kyber_pack_dk_pke(dk_pke, s_hat, param);
    return 0;
}

int kyber_pke_encrypt(uint8_t *ct, const uint8_t *ek_pke,
                      const uint8_t m[KYBER_SYMBYTES],
                      const uint8_t r[KYBER_SYMBYTES], kyber_param param) {
    lp_ctx ctx;
    uint8_t rho[KYBER_SYMBYTES];
    poly A[KYBER_MAX_K * KYBER_MAX_K];
    poly t_hat[KYBER_MAX_K];
    poly y[KYBER_MAX_K];
    poly y_hat[KYBER_MAX_K];
    poly e1[KYBER_MAX_K];
    poly e2;
    poly u[KYBER_MAX_K];
    poly v;
    poly mu_bits;
    poly mu;
    uint8_t nonce = 0;

    if(!kyber_param_is_valid(param) || !kyber_check_ek_modulus(ek_pke, param)) {
        return -1;
    }

    INIT_NTT_KYBER(&ctx);
    kyber_unpack_ek(t_hat, rho, ek_pke, param);
    kyber_matrix_expand(A, param.k, rho, &ctx);

    for(int i = 0; i < param.k; i++) {
        kyber_sample_noise(y[i], r, nonce++, param.eta1);
    }
    for(int i = 0; i < param.k; i++) {
        kyber_sample_noise(e1[i], r, nonce++, param.eta2);
    }
    kyber_sample_noise(e2, r, nonce++, param.eta2);

    polyvec_copy(y_hat, y, param.k);
    polyvec_ntt(y_hat, param.k, &ctx);

    polyvec_matrix_transpose_pointwise(u, A, param.k, param.k, y_hat, &ctx, POLYVEC_MUL_BASE);
    polyvec_invntt(u, param.k, &ctx);
    polyvec_add(u, u, e1, param.k);
    polyvec_reduce_mod(u, param.k, KYBER_Q);

    polyvec_pointwise_acc(v, t_hat, y_hat, param.k, &ctx, POLYVEC_MUL_BASE);
    poly_invntt(v, &ctx);
    POLY_ADD(v, v, e2, KYBER_N);

    kyber_byte_decode(mu_bits, m, 1);
    kyber_poly_decompress(mu, mu_bits, 1);
    POLY_ADD(v, v, mu, KYBER_N);
    poly_reduce_mod(v, KYBER_N, KYBER_Q);

    kyber_pack_ct(ct, u, v, param);
    return 0;
}

int kyber_pke_decrypt(uint8_t m[KYBER_SYMBYTES], const uint8_t *dk_pke,
                      const uint8_t *ct, kyber_param param) {
    lp_ctx ctx;
    poly s_hat[KYBER_MAX_K];
    poly u[KYBER_MAX_K];
    poly u_hat[KYBER_MAX_K];
    poly v;
    poly w;
    poly m_poly;

    if(!kyber_param_is_valid(param)) {
        return -1;
    }

    INIT_NTT_KYBER(&ctx);
    kyber_unpack_dk_pke(s_hat, dk_pke, param);
    kyber_unpack_ct(u, v, ct, param);

    polyvec_copy(u_hat, u, param.k);
    polyvec_ntt(u_hat, param.k, &ctx);
    polyvec_pointwise_acc(w, s_hat, u_hat, param.k, &ctx, POLYVEC_MUL_BASE);
    poly_invntt(w, &ctx);

    POLY_SUB(m_poly, v, w, KYBER_N);
    poly_reduce_mod(m_poly, KYBER_N, KYBER_Q);
    kyber_poly_compress(m_poly, m_poly, 1);
    kyber_byte_encode(m, m_poly, 1);
    return 0;
}
