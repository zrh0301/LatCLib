#ifndef PQC_NGCC_DRNG_H
#define PQC_NGCC_DRNG_H

#ifdef __cplusplus
extern "C" {
#endif

/* ICCS NGCC deterministic random number generator.
 *
 * Every NGCC submission draws its randomness from this DRNG so that the KAT
 * schedules are reproducible. It is a test/benchmark generator, not a
 * production CSPRNG. */
#define PQC_NGCC_DRNG_SEEDLEN 55

typedef struct {
    unsigned char V[PQC_NGCC_DRNG_SEEDLEN];
    unsigned char C[PQC_NGCC_DRNG_SEEDLEN];
    unsigned char reseed_counter[PQC_NGCC_DRNG_SEEDLEN];
} DRNG_ctx;

/* seed_len_bytes must not exceed the DRNG block size. */
int init_random_number(DRNG_ctx *drng, const unsigned char *seed,
                       unsigned long long seed_len_bytes);

/* random_number_len_bits counts output bits. */
int get_random_number(DRNG_ctx *drng, unsigned char *random_number,
                      unsigned long long random_number_len_bits);

#ifdef __cplusplus
}
#endif

#endif
