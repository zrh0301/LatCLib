/**
 * @file sample.h
 * @brief Short-noise sampling interface for S, S prime, E, E1, and E2.
 */

#ifndef SCLOUDPLUS_SAMPLE_H
#define SCLOUDPLUS_SAMPLE_H
#include <stddef.h>
#include <stdint.h>

/**
 * @brief Sample the PKE secret matrix S.
 *
 * The active BD distribution is selected by `parameters.h`. The output is
 * row-major and uses the internal q-ary representation for the short
 * coefficients.
 */
void sample_s(const uint8_t *seed, size_t seedlen, uint16_t *s);

/**
 * @brief Sample the ephemeral secret matrix S' used in encryption.
 */
void sample_sp(const uint8_t *seed, size_t seedlen, uint16_t *sp);

/**
 * @brief Sample the PKE key-generation error matrix E.
 */
void sample_e(const uint8_t *seed, size_t seedlen, uint16_t *e);

/**
 * @brief Sample the encryption error matrices E1 and E2 from one XOF stream.
 */
void sample_e12(const uint8_t *seed, size_t seedlen, uint16_t *e1, uint16_t *e2);
#endif
