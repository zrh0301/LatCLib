#ifndef PQC_RING4171_H
#define PQC_RING4171_H

#include <stdint.h>

/* Ring F_q[X]/(X^512+1) with q = 4171777 = 2^22 - 11*2^11 + 1, the parameter
 * set of the NGCC Aigis-Sig+ signature.
 *
 * Coefficients are int32_t in the reference representation: NTT outputs and
 * products are not fully reduced, and the centred range is used by rounding.
 * The routines reproduce the submitted arithmetic bit for bit so the published
 * KAT still applies. */

#define PQC_Q4171 4171777
#define PQC_Q4171_N 512

/* a / 2^32 modulo q, result in (-q, q). */
int32_t pqc_q4171_montgomery_reduce(int64_t a);
/* Barrett reduction, valid for |a| < 460138740. */
int32_t pqc_q4171_barrett_reduce(int32_t a);
/* Reduction valid for the whole int32 range. */
int32_t pqc_q4171_general_reduce(int32_t a);
/* Canonical [0, q) and centred (-q/2, q/2] representatives, |a| < q. */
int32_t pqc_q4171_amod(int32_t a);
int32_t pqc_q4171_cmod(int32_t a);

/* In-place forward NTT and inverse NTT times the Montgomery factor.
 * Both operate on exactly PQC_Q4171_N coefficients. */
void pqc_q4171_ntt(int32_t p[PQC_Q4171_N]);
void pqc_q4171_invntt_frominvmont(int32_t p[PQC_Q4171_N]);

/* Coefficient-wise helpers over n coefficients. */
void pqc_q4171_add(int32_t *c, const int32_t *a, const int32_t *b, int n);
void pqc_q4171_sub(int32_t *c, const int32_t *a, const int32_t *b, int n);
void pqc_q4171_sub_scaled(int32_t *c, const int32_t *a, const int32_t *w,
                          int32_t scale, int n);
void pqc_q4171_shiftl(int32_t *a, int n, unsigned int k);
void pqc_q4171_amod_vec(int32_t *a, int n);
void pqc_q4171_cmod_vec(int32_t *a, int n);
void pqc_q4171_barrett_vec(int32_t *a, int n);
void pqc_q4171_general_vec(int32_t *a, int n);
/* Montgomery product of two coefficient arrays. */
void pqc_q4171_pointwise_invmontgomery(int32_t *c, const int32_t *a,
                                       const int32_t *b, int n);
/* Returns 1 when some |a[i]| is outside [-bound+1, bound], 0 otherwise. */
int pqc_q4171_chknorm(const int32_t *a, int n, int32_t bound);

#endif
