#ifndef VDOO_VERIF_H
#define VDOO_VERIF_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "vdoo_config.h"
#include "vdoo_keypair.h"
#include "gf_config.h"

void vdoo_optimized_evaluation(unsigned char *y, const unsigned char *pk, const unsigned char *x);
int vdoo_verify(const uint8_t *digest, const uint8_t *signature, const pk_t *pk);

#endif /* !VDOO_VERIF_H */