#ifndef VDOO_CONFIG_H
#define VDOO_CONFIG_H

#if defined(VDOO_128)
#define VDOO_Q 16
#define VDOO_V1 67
#define VDOO_D 16
#define VDOO_O1 14
#define VDOO_O2 40

#elif defined(VDOO_256)
#define VDOO_Q 256
#define VDOO_V1 205
#define VDOO_D 10
#define VDOO_O1 17
#define VDOO_O2 68

#elif defined(VDOO_512)
#define VDOO_Q 256
#define VDOO_V1 260
#define VDOO_D 10
#define VDOO_O1 10
#define VDOO_O2 175
#endif

#define SUM_K(n) (n*(n+1)/2)
#define A_B_SUM_K(a, b) ((b+a)*(b-a+1)/2)
#define SUM_OF_SQUARES(n) ((n) * ((n) + 1) * (2*(n) + 1) / 6)
#define A_B_SUM_K_SQUARE(a, b) (SUM_OF_SQUARES(b) - SUM_OF_SQUARES((a) - 1))

#define VDOO_V2 (VDOO_V1 + VDOO_D)
#define VDOO_V3 (VDOO_V1 + VDOO_D + VDOO_O1)
#define VDOO_N (VDOO_V1 + VDOO_D + VDOO_O1 + VDOO_O2)
#define VDOO_M (VDOO_D + VDOO_O1 + VDOO_O2)

#define VDOO_DIAG_VD (A_B_SUM_K(VDOO_V1, VDOO_V1 + VDOO_D - 1))
#define VDOO_DIAG_VV (((A_B_SUM_K(VDOO_V1, VDOO_V1 + VDOO_D - 1) + A_B_SUM_K_SQUARE(VDOO_V1, VDOO_V1 + VDOO_D - 1)) / 2))
#define VDOO_O1_VV (SUM_K(VDOO_V2))
#define VDOO_O1_VO (VDOO_V2 * VDOO_O1)
#define VDOO_O2_VV (SUM_K(VDOO_V3))
#define VDOO_O2_VO (VDOO_V3 * VDOO_O2)

#if VDOO_Q == 16
#define VDOO_ELEMS_PER_BYTE 2
#define VDOO_V1_BYTE ((VDOO_V1 + 1) / 2)
#define VDOO_D_BYTE ((VDOO_D + 1) / 2)
#define VDOO_O1_BYTE ((VDOO_O1 + 1) / 2)
#define VDOO_O2_BYTE ((VDOO_O2 + 1) / 2)
#define VDOO_V2_BYTE ((VDOO_V2 + 1) / 2)
#define VDOO_V3_BYTE ((VDOO_V3 + 1) / 2)
#define VDOO_N_BYTE ((VDOO_N + 1) / 2)
#define VDOO_M_BYTE ((VDOO_M + 1) / 2)
#define VDOO_DIAG_VD_BYTE ((VDOO_DIAG_VD + 1) / 2)
#define VDOO_DIAG_VV_BYTE ((VDOO_DIAG_VV + 1) / 2)
#define VDOO_O1_VV_BYTE ((VDOO_O1_VV + 1) / 2)
#define VDOO_O1_VO_BYTE ((VDOO_O1_VO + 1) / 2)
#define VDOO_O2_VV_BYTE ((VDOO_O2_VV + 1) / 2)
#define VDOO_O2_VO_BYTE ((VDOO_O2_VO + 1) / 2)
#else
#define VDOO_ELEMS_PER_BYTE 1
#define VDOO_V1_BYTE VDOO_V1
#define VDOO_D_BYTE VDOO_D
#define VDOO_O1_BYTE VDOO_O1
#define VDOO_O2_BYTE VDOO_O2
#define VDOO_V2_BYTE VDOO_V2
#define VDOO_V3_BYTE VDOO_V3
#define VDOO_N_BYTE VDOO_N
#define VDOO_M_BYTE VDOO_M
#define VDOO_DIAG_VD_BYTE VDOO_DIAG_VD
#define VDOO_DIAG_VV_BYTE VDOO_DIAG_VV
#define VDOO_O1_VV_BYTE VDOO_O1_VV
#define VDOO_O1_VO_BYTE VDOO_O1_VO
#define VDOO_O2_VV_BYTE VDOO_O2_VV
#define VDOO_O2_VO_BYTE VDOO_O2_VO
#endif

#define LEN_SKSEED 32
#define LEN_PKSEED 32

#define SALT_BYTES 16
#define HASH_LEN 32
#define CRYPTO_BYTES (VDOO_N_BYTE + SALT_BYTES)

#define PK_ELEMENTS (VDOO_M * SUM_K(VDOO_N))
#define PK_BYTES (PK_ELEMENTS / VDOO_ELEMS_PER_BYTE)

#endif /* !VDOO_CONFIG_H */
