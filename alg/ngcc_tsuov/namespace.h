/* Symbol namespace for the TSUOV reference sources.
 *
 * The ICCS reference package exports several generic names (MGF_*, PRG_*,
 * TSUOV_*) that would collide with other modules when link-time checking or
 * whole-program builds are used. The library target and its KAT driver
 * force-include this header so every global gets a module prefix. */

#ifndef PQC_NGCC_TSUOV_NAMESPACE_H
#define PQC_NGCC_TSUOV_NAMESPACE_H

#define drng_algorithm pqc_ngcc_tsuov_drng_algorithm

#define sig_keygen pqc_ngcc_tsuov_sig_keygen
#define sig_sign pqc_ngcc_tsuov_sig_sign
#define sig_verify pqc_ngcc_tsuov_sig_verify
#define sig_get_pk_len_bytes pqc_ngcc_tsuov_sig_pk_len
#define sig_get_sk_len_bytes pqc_ngcc_tsuov_sig_sk_len
#define sig_get_sn_len_bytes pqc_ngcc_tsuov_sig_sn_len

#define MGF_init pqc_ngcc_tsuov_mgf_init
#define MGF_update pqc_ngcc_tsuov_mgf_update
#define MGF_yield pqc_ngcc_tsuov_mgf_yield
#define MGF_final pqc_ngcc_tsuov_mgf_final
#define MGF_CTX_copy pqc_ngcc_tsuov_mgf_ctx_copy

#define PRG_init_sm3 pqc_ngcc_tsuov_prg_init
#define PRG_yield_sm3 pqc_ngcc_tsuov_prg_yield
#define PRG_final_sm3 pqc_ngcc_tsuov_prg_final
#define PRG_copy_sm3 pqc_ngcc_tsuov_prg_copy
#define RejSampPRG_sm3 pqc_ngcc_tsuov_rejsamp_prg

#define TSUOV_KeyGen pqc_ngcc_tsuov_core_keygen
#define TSUOV_Sign pqc_ngcc_tsuov_core_sign
#define TSUOV_Verify pqc_ngcc_tsuov_core_verify
#define store_TSUOV_P3 pqc_ngcc_tsuov_store_p3
#define restore_TSUOV_P3 pqc_ngcc_tsuov_restore_p3

#endif /* PQC_NGCC_TSUOV_NAMESPACE_H */
