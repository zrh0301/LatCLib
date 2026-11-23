#include "api.h"
#include "owcpa.h"

#include "ngcc/auxfunc.h"
#include "ngcc/drng.h"
#include "util/ct.h"

#include <string.h>

extern DRNG_ctx pqc_ngcc_aigis_drng;

static int hash_public_key(uint8_t *out, const uint8_t *pk)
{
#if SEED_BYTES == 16
    return pseudoXOF(SEED_BYTES * 8, pk, PK_BYTES * 8, out);
#elif SEED_BYTES == 32
    return sm3hash(256, pk, PK_BYTES * 8, out);
#else
    return pseudohash(512, pk, PK_BYTES * 8, out);
#endif
}

/* The submission selects the hash by the length of the seed pair. */
static int derive_pair(uint8_t *out, const uint8_t *in, size_t len)
{
#if SEED_BYTES == 16
    return sm3hash(256, in, len * 8, out);
#else
    return pseudohash(2 * SEED_BYTES * 8, in, len * 8, out);
#endif
}

int mkem_keygen(uint8_t *pk, uint8_t *sk)
{
    uint8_t seed[SEED_BYTES], seeds[2 * SEED_BYTES];

    if (get_random_number(&pqc_ngcc_aigis_drng, seed, sizeof(seed) * 8) != 0 ||
        derive_pair(seeds, seed, sizeof(seed)) != 0 ||
        pqc_aigis_owcpa_keypair(pk, sk, seeds) != 0)
        return -1;

    memcpy(sk + PKE_SK_BYTES, pk, PK_BYTES);
    if (hash_public_key(sk + SK_BYTES - 2 * SEED_BYTES, pk) != 0)
        return -1;
    return get_random_number(&pqc_ngcc_aigis_drng, sk + SK_BYTES - SEED_BYTES,
                             SEED_BYTES * 8);
}

int mkem_enc(uint8_t *pk, uint8_t *ss, uint8_t *ct)
{
    uint8_t rnd[SEED_BYTES];
    if (get_random_number(&pqc_ngcc_aigis_drng, rnd, sizeof(rnd) * 8) != 0)
        return -1;
    return mdkem_enc(pk, rnd, ss, ct);
}

int mdkem_enc(uint8_t *pk, uint8_t *rnd, uint8_t *ss, uint8_t *ct)
{
    uint8_t key_and_coins[2 * SEED_BYTES];
    uint8_t input[2 * SEED_BYTES];

    memcpy(input, rnd, SEED_BYTES);
    if (hash_public_key(input + SEED_BYTES, pk) != 0 ||
        derive_pair(key_and_coins, input, sizeof(input)) != 0 ||
        pqc_aigis_owcpa_enc(ct, input, pk, key_and_coins + SEED_BYTES) != 0)
        return -1;

    memcpy(ss, key_and_coins, SEED_BYTES);
    return 0;
}

int mkem_dec(uint8_t *sk, uint8_t *ct, uint8_t *ss)
{
    uint8_t candidate[CT_BYTES];
    uint8_t input[2 * SEED_BYTES];
    uint8_t reject_input[SEED_BYTES + CT_BYTES];
    uint8_t reject_key[SEED_BYTES];
    uint8_t key_and_coins[2 * SEED_BYTES];
    const uint8_t *pk = sk + PKE_SK_BYTES;

    pqc_aigis_owcpa_dec(input, ct, sk);
    memcpy(input + SEED_BYTES, sk + SK_BYTES - 2 * SEED_BYTES, SEED_BYTES);
    if (derive_pair(key_and_coins, input, sizeof(input)) != 0 ||
        pqc_aigis_owcpa_enc(candidate, input, pk, key_and_coins + SEED_BYTES) != 0)
        return -1;

    const uint8_t rejected = (uint8_t)!pqc_ct_equal(ct, candidate, CT_BYTES);
    memcpy(reject_input, sk + SK_BYTES - SEED_BYTES, SEED_BYTES);
    memcpy(reject_input + SEED_BYTES, ct, CT_BYTES);
    if (pseudoXOF(SEED_BYTES * 8, reject_input, sizeof(reject_input) * 8,
                  reject_key) != 0)
        return -1;

    memcpy(ss, key_and_coins, SEED_BYTES);
    pqc_ct_select(ss, reject_key, SEED_BYTES, rejected);
    return 0;
}
