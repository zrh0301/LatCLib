/**
 * @file encode.h
 * @brief Message encoding and decoding interfaces for BW32/BW128/RBW128.
 */

#ifndef SCLOUDPLUS_ENCODE_H
#define SCLOUDPLUS_ENCODE_H
#include <stdint.h>

/**
 * @brief Integer complex coordinate used by BW/RBW message encoding.
 *
 * `real` and `imag` are never floating-point values.  In `msg_encode` they are
 * ordinary Gaussian-integer coordinates.  In `msg_decode` they are integer
 * numerators with a public implicit denominator `2^s`, exactly as in the
 * integer BDD specification; the scale is passed through the decoder rather
 * than stored in this struct.
 */
typedef struct
{
	int32_t real;
	int32_t imag;
} Complex;

/**
 * @brief Encode a KEM message into the `mbar x nbar` q-ary message matrix.
 *
 * The active instance selects BW32, BW128, or RBW128 at compile time through
 * `parameters.h`. The output matrix is row-major and contains coefficients
 * modulo `q = 2^10`, ready to be added to the PKE `C2` component.
 *
 * @param[in]  msg     Message byte string of length `scloudplus_ss`.
 * @param[out] matrix  Row-major message matrix with `mbar*nbar` coefficients.
 */
void msg_encode(const uint8_t *msg, uint16_t *matrix);

/**
 * @brief Decode a noisy q-ary message matrix back to the KEM message bytes.
 *
 * The decoder first converts q-ary coefficients to centered integer
 * numerators and then applies the selected BW/RBW integer nearest-point
 * decoder with the public scale `s = log2(q) - tau` (plus one for RBW).
 *
 * @param[in]  matrix  Row-major `mbar*nbar` matrix after PKE decryption.
 * @param[out] msg     Recovered message byte string of length `scloudplus_ss`.
 */
void msg_decode(const uint16_t *matrix, uint8_t *msg);
#endif
