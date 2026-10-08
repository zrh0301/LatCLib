#ifndef PQC_RAND_H
#define PQC_RAND_H

#include <stddef.h>
#include <stdint.h>

/* Fill every requested byte and return 0, or return nonzero on failure.
 * A callback may have partially written out before reporting failure. */
typedef int (*pqc_rng_func)(void *rng_ctx, uint8_t *out, size_t out_len);

/* OS entropy callback; rng_ctx is ignored. Returns 0 or -1. */
int pqc_system_random(void *rng_ctx, uint8_t *out, size_t out_len);

/* NULL rng selects pqc_system_random. Returns 0 or -1 and clears the entire
 * output on failure, including failures after a partial callback write.
 * NULL out is accepted only when out_len == 0. Zero-length calls succeed
 * without invoking a callback. There is no global seed or deterministic mode. */
int pqc_randombytes(uint8_t *out, size_t out_len, pqc_rng_func rng, void *rng_ctx);

/* Legacy fail-closed wrappers: invalid input or OS RNG failure calls abort().
 * New code should use pqc_randombytes and propagate its error. rand_bits stores
 * ceil(bits / 8) bytes, with unused high bits of the last byte cleared. */
void rand_bits(uint8_t *buf, int bits);
void rand_bytes(uint8_t *buf, int size);

#endif
