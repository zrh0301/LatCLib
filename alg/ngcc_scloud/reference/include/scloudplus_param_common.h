/**
 * @file scloudplus_param_common.h
 * @brief Common parameter glue for all submitted Scloud+ KEM entry points.
 *
 * Each leaf KEM entry provides only the level constants in parameters.h.
 * Family and backend are selected by compile definitions from the entry build
 * file.  Keeping this glue in one shared header avoids copying the same
 * family/backend/sample-selection logic into every level directory.
 */
#ifndef SCLOUDPLUS_PARAM_COMMON_H
#define SCLOUDPLUS_PARAM_COMMON_H

#include "parameters.h"

#define SCLOUDPLUS_PARAM_FAMILY_AES 1
#define SCLOUDPLUS_PARAM_FAMILY_SHAKE 2
#define SCLOUDPLUS_PARAM_FAMILY_SM3 3

#if defined(SCLOUDPLUS_FAMILY_AES)
#define scloudplus_param_family SCLOUDPLUS_PARAM_FAMILY_AES
#define SCLOUDPLUS_FAMILY_NAME "AES"
#elif defined(SCLOUDPLUS_FAMILY_SHAKE)
#define scloudplus_param_family SCLOUDPLUS_PARAM_FAMILY_SHAKE
#define SCLOUDPLUS_FAMILY_NAME "SHAKE"
#elif defined(SCLOUDPLUS_FAMILY_SM3)
#define scloudplus_param_family SCLOUDPLUS_PARAM_FAMILY_SM3
#define SCLOUDPLUS_FAMILY_NAME "SM3"
#else
#error "Define one of SCLOUDPLUS_FAMILY_AES, SCLOUDPLUS_FAMILY_SHAKE, or SCLOUDPLUS_FAMILY_SM3"
#endif

#if defined(SCLOUDPLUS_TIER_REFERENCE)
#define SCLOUDPLUS_BACKEND_SUFFIX ""
#elif defined(SCLOUDPLUS_BACKEND_AVX2)
#define SCLOUDPLUS_BACKEND_SUFFIX "-AVX2"
#elif defined(SCLOUDPLUS_BACKEND_NEON)
#define SCLOUDPLUS_BACKEND_SUFFIX "-NEON"
#else
#error "Define SCLOUDPLUS_TIER_REFERENCE or one optimized backend macro"
#endif

#if scloudplus_param_family == SCLOUDPLUS_PARAM_FAMILY_SM3
#define SCLOUDPLUS_SAMPLE_XOF_SM3 1
#endif

#define SCLOUDPLUS_SAMPLE_SECRET_BD scloudplus_secret_bd
#define SCLOUDPLUS_SAMPLE_ERROR_BD scloudplus_error_bd

#if (SCLOUDPLUS_SAMPLE_SECRET_BD != 2) && \
    (SCLOUDPLUS_SAMPLE_SECRET_BD != 4) && \
    (SCLOUDPLUS_SAMPLE_SECRET_BD != 6) && \
    (SCLOUDPLUS_SAMPLE_SECRET_BD != 12)
#error "Unsupported secret sampler parameter"
#endif
#if (SCLOUDPLUS_SAMPLE_ERROR_BD != 2) && \
    (SCLOUDPLUS_SAMPLE_ERROR_BD != 6) && \
    (SCLOUDPLUS_SAMPLE_ERROR_BD != 12)
#error "Unsupported error sampler parameter"
#endif

#define SYSTEM_NAME \
    "Scloudplus-" SCLOUDPLUS_LEVEL_NAME "-" SCLOUDPLUS_FAMILY_NAME \
    "-packed10" SCLOUDPLUS_BACKEND_SUFFIX

#include "api_parameters.h"

#endif /* SCLOUDPLUS_PARAM_COMMON_H */
