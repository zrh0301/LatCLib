#ifndef PQC_NGCC_AUXFUNC_H
#define PQC_NGCC_AUXFUNC_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ICCS NGCC placeholder hash primitives.
 *
 * These are the reference auxiliary functions that every NGCC submission is
 * required to use for correctness and preliminary performance evaluation.
 * They are built on SM3 (GB/T 32905-2016) and carry no security claim; later
 * NGCC rounds replace them with functions derived from the new hash track.
 *
 * All lengths are expressed in BITS, matching the ICCS interface. */

/* SM3 over a bit string. msg_len_bits may be any length. */
void pqc_ngcc_sm3_bit(const unsigned char *msg, unsigned long long msg_len_bits,
                      unsigned char *digest);

/* SM3 digest. digest_len_bits must be 256. */
int sm3hash(int digest_len_bits, const unsigned char *msg,
            unsigned long long msg_len_bits, unsigned char *digest);

/* SM3/HMAC-SM3 derived pseudo-hash. digest_len_bits is 512, 768 or 1024. */
int pseudohash(int digest_len_bits, const unsigned char *msg,
               unsigned long long msg_len_bits, unsigned char *digest);

/* Counter-mode SM3 XOF (GB/T 32918.4-2016 5.4.3). output_len_bits must be
 * below 2^40 - 2^8. */
int pseudoXOF(unsigned long long output_len_bits, const unsigned char *msg,
              unsigned long long msg_len_bits, unsigned char *output);

#ifdef __cplusplus
}
#endif

#endif
