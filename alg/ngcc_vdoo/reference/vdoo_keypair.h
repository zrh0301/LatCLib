#ifndef VDOO_KEYPAIR_H
#define VDOO_KEYPAIR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gf_config.h"
#include "ngcc/drng.h"
#include "vdoo_config.h"

/* The reference entry points and the public adapter use the same stream. */
extern DRNG_ctx pqc_ngcc_vdoo_drng;

typedef struct vdoo_secretkey
{
    unsigned char sk_seed[LEN_SKSEED];

    //unsigned char S[VDOO_M_BYTE * VDOO_M];
    unsigned char invS[VDOO_M_BYTE * VDOO_M];        // S map
    //unsigned char T[VDOO_N_BYTE * VDOO_N];
    unsigned char invT[VDOO_N_BYTE * VDOO_N];        // T map
    
    unsigned char F[VDOO_M_BYTE * SUM_K(VDOO_N)];
    //unsigned char  F_d[VDOO_DIAG_VD_BYTE + VDOO_DIAG_VV_BYTE]; // diagonal layer of F
    //unsigned char F_o1[VDOO_O1_VV_BYTE * VDOO_O1_VO_BYTE];     // layer oil1 of F
    //unsigned char F_o2[VDOO_O2_VV_BYTE * VDOO_O2_VO_BYTE];     // layer oil2 of F
} sk_t;

typedef struct vdoo_publickey
{
    unsigned char pk[VDOO_M_BYTE * SUM_K(VDOO_N)];
} pk_t;

int generate_keypair(pk_t *pk, sk_t *sk, const unsigned char *sk_seed);

#endif /* !VDOO_KEYPAIR_H */
