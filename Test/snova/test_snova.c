#include "api.h"
#include "keccak.h"
#include "symmetric.h"

#include <stdio.h>
#include <string.h>

typedef struct {
    uint8_t next;
    int fail;
} test_rng_state;

static int test_rng(void *opaque, uint8_t *out, size_t len) {
    test_rng_state *state = (test_rng_state *)opaque;
    if (state->fail) {
        return -1;
    }
    for (size_t i = 0; i < len; ++i) {
        out[i] = state->next++;
    }
    return 0;
}

static int check_hash_adapter(void) {
    static const uint8_t first[] = {0x00, 0x01, 0x02, 0x03};
    static const uint8_t second[] = {0x04, 0x05, 0x06, 0x07, 0x08};
    uint8_t joined[sizeof(first) + sizeof(second)];
    uint8_t expected[64];
    uint8_t actual[64];
    shake_t ctx;

    memcpy(joined, first, sizeof(first));
    memcpy(joined + sizeof(first), second, sizeof(second));
    shake256(expected, sizeof(expected), joined, sizeof(joined));

    snova_shake256_init(&ctx);
    snova_shake_absorb(&ctx, first, sizeof(first));
    snova_shake_absorb(&ctx, second, sizeof(second));
    snova_shake_finalize(&ctx);
    snova_shake_squeeze(&ctx, actual, 17);
    snova_shake_squeeze(&ctx, actual + 17, sizeof(actual) - 17);
    return memcmp(expected, actual, sizeof(actual)) == 0;
}

int main(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t signed_message[CRYPTO_BYTES + 32];
    uint8_t recovered[32];
    const uint8_t message[] = "SNOVA RNG adapter";
    unsigned long long signed_len = 0;
    unsigned long long recovered_len = 0;
    test_rng_state state = {0, 0};

    if (!check_hash_adapter()) {
        fprintf(stderr, "SNOVA SHAKE adapter mismatch\n");
        return 1;
    }
    memset(pk, 0xA5, sizeof(pk));
    memset(sk, 0xA5, sizeof(sk));
    state.fail = 1;
    if (crypto_sign_keypair_with_rng(pk, sk, test_rng, &state) == 0) {
        fprintf(stderr, "SNOVA key generation ignored RNG failure\n");
        return 1;
    }
    for (size_t i = 0; i < sizeof(pk); ++i) {
        if (pk[i] != 0xA5) {
            fprintf(stderr, "SNOVA modified public key after RNG failure\n");
            return 1;
        }
    }
    for (size_t i = 0; i < sizeof(sk); ++i) {
        if (sk[i] != 0xA5) {
            fprintf(stderr, "SNOVA modified secret key after RNG failure\n");
            return 1;
        }
    }
    state.fail = 0;
    if (crypto_sign_keypair_with_rng(pk, sk, test_rng, &state) != 0) {
        fprintf(stderr, "SNOVA key generation failed\n");
        return 1;
    }
    state.next = 91;
    if (crypto_sign_with_rng(signed_message, &signed_len, message, sizeof(message), sk,
                             test_rng, &state) != 0 ||
        signed_len != CRYPTO_BYTES + sizeof(message)) {
        fprintf(stderr, "SNOVA signing failed\n");
        return 1;
    }
    if (crypto_sign_open(recovered, &recovered_len, signed_message, signed_len, pk) != 0 ||
        recovered_len != sizeof(message) || memcmp(recovered, message, sizeof(message)) != 0) {
        fprintf(stderr, "SNOVA signature verification failed\n");
        return 1;
    }

    memset(signed_message, 0xA5, sizeof(signed_message));
    signed_len = 123;
    state.fail = 1;
    if (crypto_sign_with_rng(signed_message, &signed_len, message, sizeof(message), sk,
                             test_rng, &state) == 0 ||
        signed_len != 123) {
        fprintf(stderr, "SNOVA RNG failure was not propagated\n");
        return 1;
    }
    for (size_t i = 0; i < sizeof(signed_message); ++i) {
        if (signed_message[i] != 0xA5) {
            fprintf(stderr, "SNOVA modified output after RNG failure\n");
            return 1;
        }
    }

    puts("SNOVA SHAKE and RNG checks passed");
    return 0;
}
