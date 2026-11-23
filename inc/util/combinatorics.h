#ifndef PQC_COMBINATORICS_H
#define PQC_COMBINATORICS_H

#include "util/u256.h"

/* Consecutive binomial rows C(n,k), n >= k, stored as little-endian limbs.
 * offsets[n-k] and offsets[n-k+1] delimit one row in data. Each row has at
 * most four limbs. A table array is indexed by k, with entries 0 and 1 unused. */
typedef struct {
    const uint32_t *offsets;
    const uint64_t *data;
} pqc_binomial_table;

/* Compute C(n,k). Return -1 on overflow of 256 bits, leaving out unchanged.
 * k > n produces zero. All inputs and table accesses are public. */
int pqc_u256_binomial(pqc_u256 *out, uint32_t n, uint32_t k);

/* Read C(n,k) from tables[k], or compute it when tables is NULL. The caller
 * supplies tables for every k used and rows through n. */
int pqc_u256_binomial_get(pqc_u256 *out, uint32_t n, uint32_t k,
                          const pqc_binomial_table *tables);

/* Decode the zero-based lexicographic rank of a k-element subset of [0,n).
 * out has k entries, which become strictly increasing. Return -1 for an
 * invalid rank, k > n or overflowing binomial coefficient. rank and out
 * must be disjoint. tables may be NULL. k == 0 permits out == NULL. */
int pqc_combination_unrank(uint32_t *out, uint32_t n, uint32_t k,
                           const pqc_u256 *rank, const pqc_binomial_table *tables);

/* Reduce a public byte string using a positional table of len*256 residues.
 * table[i][b] is the residue contributed by byte b at position i. Every
 * residue must be smaller than the nonzero modulus. out must be disjoint
 * from modulus and table. len == 0 permits NULL input and table. */
void pqc_u256_reduce_bytes_table(pqc_u256 *out, const uint8_t *in, size_t len,
                                 const pqc_u256 table[][256], const pqc_u256 *modulus);

#endif
