#ifndef PQC_CT_H
#define PQC_CT_H

#include <stddef.h>
#include <stdint.h>

/* Constant-time byte-string helpers. Both touch every byte of their inputs, so
 * their run time depends only on the lengths, not on the data or on where the
 * first difference is. Neither is a compiler barrier: callers that compare
 * secrets should keep the result out of data-dependent branches. */
int pqc_ct_equal(const uint8_t *a, const uint8_t *b, size_t n); /* 1 if equal */
/* out = b ? x : out, with b restricted to 0 or 1. out and x may not overlap
 * unless they are identical. */
void pqc_ct_select(uint8_t *out, const uint8_t *x, size_t n, uint8_t b);

#endif
