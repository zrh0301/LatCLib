/**
 * @file parameters.h
 * @brief Scloud+ level-128 KEM default parameters.
 *
 * This file contains values fixed by the security level.  Primitive family
 * and backend choices are supplied by compiler definitions from the build entry.
 */
#ifndef SCLOUDPLUS_LEVEL_PARAMETERS_H
#define SCLOUDPLUS_LEVEL_PARAMETERS_H

#include <stdint.h>

#define scloudplus_l 128
#define SCLOUDPLUS_LEVEL_NAME "128"
#define SCLOUDPLUS_DIV_CEIL(num, den) (((num) + (den) - 1) / (den))

#define scloudplus_m 608
#define scloudplus_n 608
#define scloudplus_mbar 8
#define scloudplus_nbar 8
#define scloudplus_tau 3
#define scloudplus_secret_bd 4
#define scloudplus_error_bd 2
#define scloudplus_msg_scale_bits 3
#define scloudplus_code_rep 2

#define scloudplus_logq 10
#define scloudplus_q (1U << scloudplus_logq)
#define scloudplus_q_mask (scloudplus_q - 1U)

#define scloudplus_bw_k 5
#define scloudplus_mu 64
#define scloudplus_bw_n (1U << scloudplus_bw_k)
#define scloudplus_bw_complex (scloudplus_bw_n >> 1)
#define mat_subm SCLOUDPLUS_DIV_CEIL(scloudplus_l, scloudplus_mu)

#define scloudplus_a_bytes_per_row ((scloudplus_n * scloudplus_logq) / 8)
#define scloudplus_a_blocks_per_row SCLOUDPLUS_DIV_CEIL(scloudplus_a_bytes_per_row, 16)
#define scloudplus_a_padded_bytes_per_row (scloudplus_a_blocks_per_row * 16U)

#define scloudplus_seedA_bytes 16
#define scloudplus_hash_bytes 64
#define scloudplus_G_bytes 128
#define scloudplus_kem_z_bytes 64
#define scloudplus_pke_keygen_coins_bytes 64
#define scloudplus_pke_keygen_r1_bytes 64
#define scloudplus_pke_keygen_r2_bytes 64
#define scloudplus_pke_keygen_expand_bytes \
    (scloudplus_seedA_bytes + scloudplus_pke_keygen_r1_bytes + \
     scloudplus_pke_keygen_r2_bytes)
#define scloudplus_pke_enc_coins_bytes 64
#define scloudplus_pke_enc_r1_bytes 64
#define scloudplus_pke_enc_r2_bytes 64
#define scloudplus_pke_enc_expand_bytes \
    (scloudplus_pke_enc_r1_bytes + scloudplus_pke_enc_r2_bytes)
#define scloudplus_ss (scloudplus_l >> 3)

#define scloudplus_pk \
    (SCLOUDPLUS_DIV_CEIL(scloudplus_m * scloudplus_nbar * scloudplus_logq, 8) + \
     scloudplus_seedA_bytes)
#define scloudplus_c1 \
    SCLOUDPLUS_DIV_CEIL(scloudplus_mbar * scloudplus_n * scloudplus_logq, 8)
#define scloudplus_c2 \
    SCLOUDPLUS_DIV_CEIL(scloudplus_mbar * scloudplus_nbar * scloudplus_logq, 8)
#define scloudplus_ctx (scloudplus_c1 + scloudplus_c2)
#define scloudplus_sk_coeff_bits 2
#define scloudplus_sk_coeff_bias 1
#define scloudplus_pke_sk \
    SCLOUDPLUS_DIV_CEIL(scloudplus_n * scloudplus_nbar * scloudplus_sk_coeff_bits, 8)
#define scloudplus_kem_sk \
    (scloudplus_pke_sk + scloudplus_pk + scloudplus_hash_bytes + \
     scloudplus_kem_z_bytes)

#endif /* SCLOUDPLUS_LEVEL_PARAMETERS_H */
