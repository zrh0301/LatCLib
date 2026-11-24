/**
 * @file api_parameters.h
 * @brief Shared API_PKC-facing parameter facade for Scloud+ KEM selections.
 *
 * This file is included by scloudplus_param_common.h after the leaf level
 * parameters and family/backend selection macros have been defined.  It
 * exposes the names required by API_PKC helpers and tests without duplicating
 * aliases in every leaf directory.
 */
#ifndef SCLOUDPLUS_API_PARAMETERS_H
#define SCLOUDPLUS_API_PARAMETERS_H

#define SCLOUDPLUS_SECURITY_BITS scloudplus_l
#define SCLOUDPLUS_SECURITY_BYTES (scloudplus_l >> 3)
#define SCLOUDPLUS_SHARED_SECRET_BYTES scloudplus_ss

#define SCLOUDPLUS_PARAMETER_M scloudplus_m
#define SCLOUDPLUS_PARAMETER_N scloudplus_n
#define SCLOUDPLUS_PARAMETER_MBAR scloudplus_mbar
#define SCLOUDPLUS_PARAMETER_NBAR scloudplus_nbar
#define SCLOUDPLUS_PARAMETER_LOGQ scloudplus_logq
#define SCLOUDPLUS_PARAMETER_Q scloudplus_q
#define SCLOUDPLUS_PARAMETER_TAU scloudplus_tau
#define SCLOUDPLUS_PARAMETER_MU scloudplus_mu
#define SCLOUDPLUS_PARAMETER_CODE_REP scloudplus_code_rep
#define SCLOUDPLUS_SEEDA_BYTES scloudplus_seedA_bytes
#define SCLOUDPLUS_HASH_BYTES scloudplus_hash_bytes

#define SCLOUDPLUS_PUBLIC_KEY_BYTES scloudplus_pk
#define SCLOUDPLUS_SECRET_KEY_BYTES scloudplus_kem_sk
#define SCLOUDPLUS_CIPHERTEXT_BYTES scloudplus_ctx

#endif /* SCLOUDPLUS_API_PARAMETERS_H */
