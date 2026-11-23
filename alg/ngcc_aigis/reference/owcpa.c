#include "owcpa.h"

#include "ngcc/auxfunc.h"
#include "ngcc/drng.h"
#include "poly/ntt16.h"
#include "poly/poly16.h"

#include <string.h>

#define NTT_BLOCKS (PARAM_N / NTT_DIM)
#define MESSAGE_CHUNK_BYTES (NTT_DIM / (8 * MESSAGE_COPIES))

static const fp16_ctx field = {
    PARAM_Q, QINV, 32768u / PARAM_Q,
    (uint16_t)(((65536u % PARAM_Q) * (65536u % PARAM_Q)) % PARAM_Q)
};

_Static_assert(PARAM_N == MSG_BYTES * 8 * MESSAGE_COPIES, "message layout");
_Static_assert(MSG_BYTES % MESSAGE_CHUNK_BYTES == 0, "whole message chunks");
_Static_assert(NTT_BLOCKS == 4 || NTT_BLOCKS == 8 || NTT_BLOCKS == 16,
               "supported split NTT sizes");

/* Each quarter has a separate DRNG seeded with seed || quarter. The
 * submission requests 96 or 32 bytes per draw and discards partial pairs. */
static int expand_a(int16_t a[PARAM_N], const uint8_t seed[SEED_BYTES])
{
    const size_t quarter = PARAM_N / 4;
    uint8_t tagged_seed[SEED_BYTES + 1];
    uint8_t bytes[96];

    memcpy(tagged_seed, seed, SEED_BYTES);
    for (uint8_t part = 0; part < 4; ++part) {
        DRNG_ctx rng;
        tagged_seed[SEED_BYTES] = part;
        if (init_random_number(&rng, tagged_seed, sizeof(tagged_seed)) != 0)
            return -1;
        for (size_t count = 0; count < quarter;) {
            const size_t length = quarter - count > 21 ? sizeof(bytes) : 32;
            if (get_random_number(&rng, bytes, length * 8) != 0)
                return -1;
            count += fp16_poly_rej_uniform12(a + part * quarter + count,
                                             quarter - count, bytes, length, &field);
        }
    }
    return 0;
}

static int sample_noise(int16_t out[PARAM_N], const uint8_t seed[SEED_BYTES],
                        uint8_t nonce, unsigned eta)
{
    uint8_t tagged_seed[SEED_BYTES + 1];
    uint8_t bytes[PARAM_N * ETA_MAX / 4];
    const size_t length = PARAM_N * eta / 4;

    memcpy(tagged_seed, seed, SEED_BYTES);
    tagged_seed[SEED_BYTES] = nonce;
    if (pseudoXOF(length * 8, tagged_seed, sizeof(tagged_seed) * 8, bytes) != 0)
        return -1;
    fp16_poly_cbd(out, PARAM_N, bytes, eta);
    return 0;
}

int pqc_aigis_owcpa_keypair(uint8_t *pk, uint8_t *sk, const uint8_t *seeds)
{
    const uint8_t *noise_seed = seeds;
    const uint8_t *public_seed = seeds + SEED_BYTES;
    int16_t a[PARAM_N], s[PARAM_N], e[PARAM_N], t[PARAM_N];

    if (expand_a(a, public_seed) != 0 ||
        sample_noise(s, noise_seed, 0, ETA_S) != 0 ||
        sample_noise(e, noise_seed, 1, ETA_E) != 0)
        return -1;

    pqc_ntt128_blocks(s, NTT_BLOCKS);
    pqc_ntt128_blocks_mont_mul(t, s, a, NTT_BLOCKS);
    fp16_poly_tomont(t, PARAM_N, &field);
    pqc_invntt128_blocks(t, NTT_BLOCKS);
    fp16_poly_add(t, t, e, PARAM_N);
    fp16_poly_caddq(t, PARAM_N, &field);
    fp16_poly_caddq(s, PARAM_N, &field);

    /* PKE secret key is the encoded NTT secret. Public key is t || seed. */
    fp16_poly_encode(sk, s, PARAM_N, QBITS);
    fp16_poly_compress(pk, t, PARAM_N, BITS_PK, &field);
    memcpy(pk + PK_COEFF_BYTES, public_seed, SEED_BYTES);
    return 0;
}

int pqc_aigis_owcpa_enc(uint8_t *ct, const uint8_t *m, const uint8_t *pk,
                       const uint8_t *coins)
{
    int16_t a[PARAM_N], t[PARAM_N], r[PARAM_N], e1[PARAM_N], e2[PARAM_N];
    int16_t u[PARAM_N], v[PARAM_N], encoded[PARAM_N];

    if (expand_a(a, pk + PK_COEFF_BYTES) != 0 ||
        sample_noise(r, coins, 0, ETA_S) != 0 ||
        sample_noise(e1, coins, 1, ETA_E) != 0 ||
        sample_noise(e2, coins, 2, ETA_E) != 0)
        return -1;

    fp16_poly_decompress(t, pk, PARAM_N, BITS_PK, &field);
    fp16_poly_encode_repeated_bits(encoded, m, MSG_BYTES, MESSAGE_CHUNK_BYTES,
                                   MESSAGE_COPIES, &field);
    pqc_ntt128_blocks(t, NTT_BLOCKS);
    pqc_ntt128_blocks(r, NTT_BLOCKS);

    pqc_ntt128_blocks_mont_mul(u, r, a, NTT_BLOCKS);
    fp16_poly_tomont(u, PARAM_N, &field);
    pqc_invntt128_blocks(u, NTT_BLOCKS);
    fp16_poly_add(u, u, e1, PARAM_N);
    fp16_poly_caddq(u, PARAM_N, &field);

    pqc_ntt128_blocks_mont_mul(v, t, r, NTT_BLOCKS);
    fp16_poly_tomont(v, PARAM_N, &field);
    pqc_invntt128_blocks(v, NTT_BLOCKS);
    fp16_poly_add(v, v, e2, PARAM_N);
    fp16_poly_sub(v, v, encoded, PARAM_N);
    fp16_poly_caddq2(v, PARAM_N, &field);

    /* Ciphertext concatenates the two rounded coefficient strings. */
    fp16_poly_compress(ct, u, PARAM_N, BITS_C1, &field);
    fp16_poly_compress(ct + CT_U_BYTES, v, PARAM_N, BITS_C2, &field);
    return 0;
}

void pqc_aigis_owcpa_dec(uint8_t *m, const uint8_t *ct, const uint8_t *sk)
{
    int16_t u[PARAM_N], v[PARAM_N], s[PARAM_N], decoded[PARAM_N];

    fp16_poly_decompress(u, ct, PARAM_N, BITS_C1, &field);
    fp16_poly_decompress(v, ct + CT_U_BYTES, PARAM_N, BITS_C2, &field);
    fp16_poly_decode(s, sk, PARAM_N, QBITS);
    pqc_ntt128_blocks(u, NTT_BLOCKS);
    pqc_ntt128_blocks_mont_mul(decoded, s, u, NTT_BLOCKS);
    fp16_poly_tomont(decoded, PARAM_N, &field);
    pqc_invntt128_blocks(decoded, NTT_BLOCKS);
    fp16_poly_sub(decoded, decoded, v, PARAM_N);
    fp16_poly_caddq2(decoded, PARAM_N, &field);
    fp16_poly_decode_repeated_bits(m, decoded, MSG_BYTES, MESSAGE_CHUNK_BYTES,
                                   MESSAGE_COPIES, &field);
}
