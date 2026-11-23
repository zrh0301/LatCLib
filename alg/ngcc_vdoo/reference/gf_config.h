#ifndef PQC_VDOO_GF_CONFIG_H
#define PQC_VDOO_GF_CONFIG_H

#include "vdoo_config.h"
#include "gf/composite_matrix.h"
#include "gf/gf16_vector.h"

#if VDOO_Q == 16
#define VDOO_FIELD GF_COMPOSITE_16
#define gfv_get_ele gf16_packed_get
#define gfv_set_ele gf16_packed_set
#define gfv_mul gf16_composite_mul
#define gfv_inv gf16_composite_inv
#else
#define VDOO_FIELD GF_COMPOSITE_256
#define gfv_get_ele(v, i) ((v)[i])
#define gfv_set_ele(v, i, value) ((v)[i] = (value))
#define gfv_mul gf256_composite_mul
#define gfv_inv gf256_composite_inv
#endif

#endif
