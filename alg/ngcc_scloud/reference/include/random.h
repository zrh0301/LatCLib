/**
 * @file random.h
 * @brief randombytes declaration used by PKE key generation and tests.
 */

#ifndef __RANDOM_H__
#define __RANDOM_H__

/**
 * @brief Fill a buffer with operating-system randomness.
 *
 * Correctness tests and the benchmark use this routine for non-KAT random
 * sampling. KAT generation uses the deterministic API_PKC DRNG instead.
 *
 * @return 0 on success, negative values on platform-specific failures.
 */
int randombytes(unsigned char *buffer, unsigned int size);
#endif
