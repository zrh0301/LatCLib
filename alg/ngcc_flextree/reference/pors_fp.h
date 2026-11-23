#ifndef SPX_PORS_FP_H
#define SPX_PORS_FP_H

#include <stdint.h>

#include "params.h"

/**
 * Signs a message digest using PORS+FP.
 */
void pors_fp_sign(unsigned char *sig, unsigned char *pk,
                  const unsigned char *m,
                  const spx_ctx* ctx,
                  const uint32_t pors_fp_addr[8],
                  const unsigned char *R);

/**
 * Recomputes the PORS+FP root from a signature.
 */
void pors_fp_pk_from_sig(unsigned char *pk,
                         const unsigned char *sig, const unsigned char *m,
                         const spx_ctx* ctx,
                         const uint32_t pors_fp_addr[8],
                         const unsigned char *R);

#endif
