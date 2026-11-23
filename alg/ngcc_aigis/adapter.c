/* Adapter from the ICCS NGCC KEM ABI to the library KEM API.
 *
 * Reference sources: ICCS NGCC round-1 submission package, "Aigis-Enc+"
 * (encryption track), Reference_Implementation/Aigis-Enc+-I. Only the first
 * parameter set is shipped; the other two are separate directories whose
 * sources differ only in params.h.
 *
 * Local changes to the reference layer:
 *   - kem.c reads the implicit-rejection secret from sk + SK_BYTES -
 *     SEED_BYTES instead of sk - SEED_BYTES (kem-01-2);
 *   - kem.c overwrites the returned shared secret rather than the pre-key
 *     when re-encryption fails, restoring implicit rejection (kem-01-1);
 *   - sampling and hashing call the shared polynomial and NGCC primitives;
 *   - all external symbols carry the pqc_aigis_ prefix (namespace.h).
 *
 * The KAT driver seeds pqc_ngcc_aigis_drng directly. The adapter seeds it from the
 * system RNG before key generation and encapsulation (decapsulation draws no
 * randomness after the fixes above), so library callers never share the KAT
 * stream. */

#include "api.h"
#include "core/pqc_kem.h"
#include "util/rand.h"
#include "ngcc/drng.h"

/* Referenced by the reference KEM and seeded directly by the KAT. */
DRNG_ctx pqc_ngcc_aigis_drng;

#define PQC_AIGIS_NAME "Aigis-Enc+-I"

static int seed_drng(void)
{
    uint8_t seed[32];
    if (pqc_randombytes(seed, sizeof(seed), NULL, NULL) != 0)
        return -1;
    if (init_random_number(&pqc_ngcc_aigis_drng, seed, sizeof(seed)) != 0)
        return -1;
    return 0;
}

static int keypair(uint8_t *pk, uint8_t *sk)
{
    if (seed_drng() != 0)
        return -1;
    return mkem_keygen(pk, sk);
}

static int encaps(uint8_t *ct, uint8_t *ss, const uint8_t *pk)
{
    if (seed_drng() != 0)
        return -1;
    return mkem_enc((uint8_t *)pk, ss, ct);
}

static int decaps(uint8_t *ss, const uint8_t *ct, const uint8_t *sk)
{
    return mkem_dec((uint8_t *)sk, (uint8_t *)ct, ss);
}

static const pqc_kem algorithm = {
    PQC_AIGIS_NAME,
    KEM_PUBLICKEYBYTES,
    KEM_SECRETKEYBYTES,
    KEM_CIPHERTEXTBYTES,
    KEM_BYTES,
    keypair,
    encaps,
    decaps,
};

const pqc_kem *pqc_ngcc_aigis_algorithm(void)
{
    return &algorithm;
}
