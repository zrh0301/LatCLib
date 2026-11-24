// tsuov_params.h

#pragma once


#ifndef TSUOV_VARIANT
#define TSUOV_VARIANT 128
#endif

#if (TSUOV_VARIANT == 128)
#define ALGORITHM_INSTANCE    "TSUOV_128"
#define TSUOV_q 31
#define TSUOV_V 60
#define TSUOV_O 4
#define TSUOV_L 2
#define TSUOV_k 11
#define TSUOV_m1 61
#define TSUOV_m2 88
#define TSUOV_SEED_LEN 16
#define TSUOV_Pi1_RS_LEN 2469 
#define TSUOV_Pi2_RS_LEN 359  
#define TSUOV_kv_RS_LEN 923   
#define TSUOV_t_RS_LEN 86     
#define F_COEF {4, 3, 1, 0}         // f(z) = z^88 + 0*z^3 + 30*z^2 + 28*z + 27

#elif (TSUOV_VARIANT == 256)
#define ALGORITHM_INSTANCE    "TSUOV_256"
#define TSUOV_q 31
#define TSUOV_V 114
#define TSUOV_O 5
#define TSUOV_L 2
#define TSUOV_k 16
#define TSUOV_m1 114
#define TSUOV_m2 160
#define TSUOV_SEED_LEN 32
#define TSUOV_Pi1_RS_LEN 8738 
#define TSUOV_Pi2_RS_LEN 841  
#define TSUOV_kv_RS_LEN 2517  
#define TSUOV_t_RS_LEN 162    
#define F_COEF {6, 4, 1, 0}        // f(z)=z^160+ 0*z^3 +30*z^2+ 27*z+ 25

#elif (TSUOV_VARIANT == 512)
#define ALGORITHM_INSTANCE    "TSUOV_512"
#define TSUOV_q 31
#define TSUOV_V 218
#define TSUOV_O 12
#define TSUOV_L 2
#define TSUOV_k 10
#define TSUOV_m1 218
#define TSUOV_m2 240
#define TSUOV_SEED_LEN 64
#define TSUOV_Pi1_RS_LEN 31536
#define TSUOV_Pi2_RS_LEN 3665 
#define TSUOV_kv_RS_LEN 3082  
#define TSUOV_t_RS_LEN 264    
#define F_COEF {15, 12, 1, 0}       // f(z) = z^240 + 0*z^3 + 30*z^2 + 19*z + 16

#else
#  error "Unsupported: TSUOV_VARIANT == " # TSUOV_VARIANT
#endif

/* Keep legacy TSUOV_security_strength_category in sync with TSUOV_VARIANT
 * (multiply_E.c, tsuov.h EVP_* selection still use the category macro). */
#undef TSUOV_security_strength_category
#if (TSUOV_VARIANT == 128)
#define TSUOV_security_strength_category 1
#elif (TSUOV_VARIANT == 256)
#define TSUOV_security_strength_category 3
#elif (TSUOV_VARIANT == 512)
#define TSUOV_security_strength_category 5
#endif

// --------------------------------------------------------------------------------

// define other useful params

#define TSUOV_ceil_log_2_q 5
/// length of salt for a signature, in # bytes
#define TSUOV_SALT_LEN TSUOV_SEED_LEN

#define TSUOV_MU_LEN 64

#define TSUOV_N (TSUOV_V + TSUOV_O)
#define TSUOV_n (TSUOV_N * TSUOV_L)
#define TSUOV_o (TSUOV_O * TSUOV_L)

#define TSUOV_Pi1_LEN (TSUOV_L * TSUOV_V * (TSUOV_V + 1)/2)  
#define TSUOV_Pi2_LEN (TSUOV_L * TSUOV_V * TSUOV_O)
#define TSUOV_Pi3_LEN (TSUOV_L * TSUOV_O * (TSUOV_O + 1)/2)
#define TSUOV_kv_LEN (TSUOV_L * TSUOV_V * TSUOV_k)
#define TSUOV_t_LEN (TSUOV_m2)

#define ALIGN_8BIT   __attribute((__aligned__( 1))) ; //   8 bit
#define ALIGN_16BIT  __attribute((__aligned__( 2))) ; //  16 bit
#define ALIGN_32BIT  __attribute((__aligned__( 4))) ; //  32 bit
#define ALIGN_64BIT  __attribute((__aligned__( 8))) ; //  64 bit
#define ALIGN_128BIT __attribute((__aligned__(16))) ; // 128 bit
#define ALIGN_256BIT __attribute((__aligned__(32))) ; // 256 bit


#define TYPEDEF_STRUCT(TYPE_NAME, BODY)             \
typedef struct TYPE_NAME ## _t {                    \
  BODY                                              \
} TYPE_NAME ## _s, * TYPE_NAME ## _p, TYPE_NAME [1]



