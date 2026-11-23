#ifndef PQC_AIGIS_PARAMS_H
#define PQC_AIGIS_PARAMS_H

#define NTT_DIM 128
#define PARAM_Q 3329
#define QINV 62209
#define QBITS 12
#define BITS_PK 10
#define BITS_C1 10
#define MESSAGE_COPIES 4

#if PARAMS == 1
#define PARAM_N 512
#define ETA_S 5
#define ETA_E 6
#define BITS_C2 4
#define SEED_BYTES 16
#elif PARAMS == 2
#define PARAM_N 1024
#define ETA_S 2
#define ETA_E 6
#define BITS_C2 3
#define SEED_BYTES 32
#elif PARAMS == 3
#define PARAM_N 2048
#define ETA_S 2
#define ETA_E 3
#define BITS_C2 4
#define SEED_BYTES 64
#else
#error "PARAMS must be in {1,2,3}"
#endif

#define ETA_MAX (ETA_S > ETA_E ? ETA_S : ETA_E)
#define MSG_BYTES SEED_BYTES

#define PKE_SK_BYTES (PARAM_N * QBITS / 8)
#define PK_COEFF_BYTES (PARAM_N * BITS_PK / 8)
#define CT_U_BYTES (PARAM_N * BITS_C1 / 8)
#define PK_BYTES (PK_COEFF_BYTES + SEED_BYTES)
#define SK_BYTES (PKE_SK_BYTES + PK_BYTES + 2 * SEED_BYTES)
#define CT_BYTES (CT_U_BYTES + PARAM_N * BITS_C2 / 8)

#endif
