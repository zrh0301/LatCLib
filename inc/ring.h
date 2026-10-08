#ifndef PQC_RING_H
#define PQC_RING_H

#include "fp32.h"

typedef enum {
    RING_XN_PLUS_ONE = 1
} ring_modulus;

/* Mathematical ring only: F_q[X]/(X^n+1). No transform tables or RNG.
 * Value-owned context: copying does not retain pointers to the original.
 * Treat all fields as immutable after successful initialization. */
typedef struct {
    fp32_ctx field;
    int n;
    ring_modulus modulus;
} ring_ctx;

/* Validate the supplied field and positive dimension. Failure leaves out
 * unchanged. Individual arithmetic backends may support fewer dimensions. */
int ring_init_negacyclic(ring_ctx *out, const fp32_ctx *field, int n);

#endif
