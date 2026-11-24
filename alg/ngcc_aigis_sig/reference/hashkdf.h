#ifndef hashkdf_h
#define hashkdf_h

/* Key derivation and hashing for Aigis-Sig+, restricted to the ICCS build
 * (USE_ICCS) that generated the published KAT vector.
 *
 * Local changes to the reference layer:
 *   - the SHA3/SM3/DRBG alternatives are dropped; the ICCS pseudoXOF,
 *     pseudohash and sm3hash primitives from inc/ngcc/auxfunc.h are used;
 *   - all lengths are 64-bit, so |M| >= 2^29 no longer wraps the bit length
 *     (sign-01-9);
 *   - kdf_init reports DRNG instantiation failure instead of returning void
 *     (sign-01-10). */

#include <stdint.h>
#include "params.h"
#include "ngcc/auxfunc.h"
#include "ngcc/drng.h"

#define KDF_RATE 32
#define KDF128RATE 32
#define KDF256RATE 32
#define KDF512RATE 32
typedef DRNG_ctx kdfstate;

int kdf_init(kdfstate *state, const uint8_t *input, unsigned long long inlen);
int kdf_squeezeblocks(uint8_t *output, unsigned long long nblocks,
                      kdfstate *state);

int kdf128(uint8_t *out, unsigned long long outlen, const uint8_t *in,
           unsigned long long inlen);
int kdf256(uint8_t *out, unsigned long long outlen, const uint8_t *in,
           unsigned long long inlen);
int kdf512(uint8_t *out, unsigned long long outlen, const uint8_t *in,
           unsigned long long inlen);

int hash256(uint8_t *out, const uint8_t *in, unsigned long long inlen);
int hash512(uint8_t *out, const uint8_t *in, unsigned long long inlen);
int hash1024(uint8_t *out, const uint8_t *in, unsigned long long inlen);

#endif
