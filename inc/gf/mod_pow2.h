#ifndef PQC_MOD_POW2_H
#define PQC_MOD_POW2_H

#include <stdint.h>

/* Arithmetic modulo a power of two, q = 2^logq.
 *
 * Reduction is a public bit mask (no secret-dependent branch), which makes
 * these helpers usable on ciphertext coefficients. They are shared by the
 * module-LWR style candidates (currently Scloud+). */

uint16_t pqc_mod_pow2(uint32_t value, unsigned int logq);

/* Round a q-ary coefficient down to `bits` bits, with round-half-up. */
int32_t pqc_round_pow2(uint16_t value, unsigned int logq, unsigned int bits);

#endif
