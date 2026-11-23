#ifndef VDOO_SIGN_H
#define VDOO_SIGN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gf_config.h"
#include "vdoo_keypair.h"
#include "vdoo_config.h"

int vdoo_sign(uint8_t * signature, const sk_t * sk, const uint8_t * _digest);

#endif /* !VDOO_SIGN_H */