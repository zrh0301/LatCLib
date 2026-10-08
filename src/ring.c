#include "ring.h"
#include <stddef.h>

int ring_init_negacyclic(ring_ctx *out, const fp32_ctx *field, int n)
{
    fp32_ctx checked;
    if (out == NULL || field == NULL || n <= 0 ||
        fp32_init(&checked, field->q) != 0 ||
        checked.qinv != field->qinv || checked.mont_r != field->mont_r ||
        checked.mont_r2 != field->mont_r2) {
        return -1;
    }
    ring_ctx result = {checked, n, RING_XN_PLUS_ONE};
    *out = result;
    return 0;
}
