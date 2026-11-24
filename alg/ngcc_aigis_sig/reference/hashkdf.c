#include "hashkdf.h"

/* ICCS build: the pseudo-random functions are the submission's auxfunc
 * primitives. The reference wrapped them without checking their status; the
 * wrappers now return it so callers can abort instead of hashing from an
 * all-zero DRBG state (sign-01-10) or from a failed allocation (sign-01-9). */

int kdf_init(kdfstate *state, const uint8_t *input, unsigned long long inlen)
{
    return init_random_number(state, input, inlen);
}

int kdf_squeezeblocks(uint8_t *output, unsigned long long nblocks,
                      kdfstate *state)
{
    return get_random_number(state, output, nblocks * 256);
}

static int xof_check(unsigned long long outlen, unsigned long long inlen)
{
    /* pseudoXOF takes bit lengths in unsigned long long; reject values that
     * cannot be represented instead of wrapping them. */
    if (outlen > (~0ULL) / 8 || inlen > (~0ULL) / 8)
        return -1;
    return 0;
}

int kdf128(uint8_t *out, unsigned long long outlen, const uint8_t *in,
           unsigned long long inlen)
{
    if (xof_check(outlen, inlen) != 0)
        return -1;
    return pseudoXOF(outlen * 8, in, inlen * 8, out);
}

int kdf256(uint8_t *out, unsigned long long outlen, const uint8_t *in,
           unsigned long long inlen)
{
    if (xof_check(outlen, inlen) != 0)
        return -1;
    return pseudoXOF(outlen * 8, in, inlen * 8, out);
}

int kdf512(uint8_t *out, unsigned long long outlen, const uint8_t *in,
           unsigned long long inlen)
{
    if (xof_check(outlen, inlen) != 0)
        return -1;
    return pseudoXOF(outlen * 8, in, inlen * 8, out);
}

int hash256(uint8_t *out, const uint8_t *in, unsigned long long inlen)
{
    if (inlen > (~0ULL) / 8)
        return -1;
    return sm3hash(256, in, inlen * 8, out);
}

int hash512(uint8_t *out, const uint8_t *in, unsigned long long inlen)
{
    if (inlen > (~0ULL) / 8)
        return -1;
    return pseudohash(512, in, inlen * 8, out);
}

int hash1024(uint8_t *out, const uint8_t *in, unsigned long long inlen)
{
    if (inlen > (~0ULL) / 8)
        return -1;
    return pseudohash(1024, in, inlen * 8, out);
}
