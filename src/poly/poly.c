#include "poly.h"

/* Legacy qbits argument is retained for source compatibility. The field
 * reference backend provides exact reduction for the full input width. */
void poly_mod_barrett_lssi(lssi_t *c, const lssi_t *a, lssi_t q,
                           int qbits, int coefs) {
    (void)qbits;
    for (int i = 0; i < coefs; ++i) c[i] = fp32_center_raw(a[i], q);
}

void poly_mod_barrett_lsdi(lssi_t *c, const lsdi_t *a, lssi_t q,
                           int qbits, int coefs) {
    (void)qbits;
    for (int i = 0; i < coefs; ++i) c[i] = fp32_center_raw(a[i], q);
}

int poly_mod_equal(const lssi_t *a, const lssi_t *b, const lssi_t q, int qbits, int coefs){
    poly c, d;
    poly_mod_barrett_lssi(c, a, q, qbits, coefs);
    poly_mod_barrett_lssi(d, b, q, qbits, coefs);
    return POLY_EQ(c, d, coefs);
}

int poly_mod_is_zero(const lssi_t *a, const lssi_t q, int qbits, int coefs){
    poly c;
    poly_mod_barrett_lssi(c, a, q, qbits, coefs);
    return POLY_IS_ZERO(c, coefs);
}

void poly_neg(lssi_t *c, const lssi_t *a, lssi_t q, int qbits, int coefs) {
    (void)qbits;
    for (int i = 0; i < coefs; ++i) c[i] = fp32_center_raw(-(int64_t)a[i], q);
}
