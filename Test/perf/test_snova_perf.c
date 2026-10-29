#include "test.h"
#include "api.h"
#include "symmetric.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define KEYGEN_ITERATIONS 5
#define SIGN_ITERATIONS 20
#define VERIFY_ITERATIONS 20
#define SHAKE_ITERATIONS 10000
#define CORE_ITERATIONS 20

typedef struct {
    uint32_t state;
} perf_rng_state;

static int perf_rng(void *opaque, uint8_t *out, size_t out_len) {
    perf_rng_state *state = (perf_rng_state *)opaque;
    for (size_t i = 0; i < out_len; ++i) {
        uint32_t x = state->state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state->state = x;
        out[i] = (uint8_t)x;
    }
    return 0;
}

static int bench_keygen(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    perf_rng_state rng = {UINT32_C(0x6d5a56a9)};
    test_perf_mark start = test_perf_begin();

    for (unsigned int i = 0; i < KEYGEN_ITERATIONS; ++i) {
        if (crypto_sign_keypair_with_rng(pk, sk, perf_rng, &rng) != 0) return 0;
    }
    test_perf_print("SNOVA key generation", test_perf_end(start), KEYGEN_ITERATIONS);
    return 1;
}

static int bench_sign_or_verify(const char *operation) {
    static const uint8_t message[] = "SNOVA performance benchmark";
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t signed_message[CRYPTO_BYTES + sizeof(message)];
    uint8_t recovered[sizeof(message)];
    unsigned long long signed_len = 0;
    unsigned long long recovered_len = 0;
    perf_rng_state rng = {UINT32_C(0x6d5a56a9)};
    test_perf_mark start;

    if (crypto_sign_keypair_with_rng(pk, sk, perf_rng, &rng) != 0) return 0;
    if (strcmp(operation, "sign") == 0) {
        start = test_perf_begin();
        for (unsigned int i = 0; i < SIGN_ITERATIONS; ++i) {
            if (crypto_sign_with_rng(signed_message, &signed_len, message,
                                     sizeof(message), sk, perf_rng, &rng) != 0) return 0;
        }
        test_perf_print("SNOVA sign", test_perf_end(start), SIGN_ITERATIONS);
        return signed_len == CRYPTO_BYTES + sizeof(message);
    }

    if (crypto_sign_with_rng(signed_message, &signed_len, message, sizeof(message),
                             sk, perf_rng, &rng) != 0) return 0;
    start = test_perf_begin();
    for (unsigned int i = 0; i < VERIFY_ITERATIONS; ++i) {
        if (crypto_sign_open(recovered, &recovered_len, signed_message,
                             signed_len, pk) != 0 ||
            recovered_len != sizeof(message) ||
            memcmp(recovered, message, sizeof(message)) != 0) return 0;
    }
    test_perf_print("SNOVA verify", test_perf_end(start), VERIFY_ITERATIONS);
    return 1;
}

static void bench_shake(void) {
    uint8_t input[168];
    uint8_t output[64];
    shake_t shake;
    memset(input, 0x5a, sizeof(input));
    test_perf_mark start = test_perf_begin();
    for (unsigned int i = 0; i < SHAKE_ITERATIONS; ++i) {
        input[0] = (uint8_t)i;
        snova_shake256_init(&shake);
        snova_shake_absorb(&shake, input, sizeof(input));
        snova_shake_finalize(&shake);
        snova_shake_squeeze(&shake, output, sizeof(output));
    }
    volatile uint8_t checksum = output[0];
    (void)checksum;
    test_perf_print("SNOVA SHAKE256 (168B)", test_perf_end(start), SHAKE_ITERATIONS);
}

static int bench_core(const char *operation) {
    static const uint8_t digest[BYTES_DIGEST] = {0x53, 0x4e, 0x4f, 0x56, 0x41};
    uint8_t seed[SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE] = {0};
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t salt[BYTES_SALT] = {0};
    expanded_SK skx;
    expanded_PK pkx;
    test_perf_mark start;

    for (size_t i = 0; i < sizeof(seed); ++i) seed[i] = (uint8_t)(i * 29u + 7u);
    if (SNOVA_NAMESPACE(genkeys)(pk, sk, seed) != 0 ||
        SNOVA_NAMESPACE(sk_expand)(&skx, sk) != 0 ||
        SNOVA_NAMESPACE(pk_expand)(&pkx, pk) != 0 ||
        SNOVA_NAMESPACE(sign)(&skx, sig, digest, sizeof(digest), salt) != 0 ||
        SNOVA_NAMESPACE(verify)(&pkx, sig, digest, sizeof(digest)) != 0) {
        fprintf(stderr, "SNOVA core benchmark setup failed\n");
        return 0;
    }

    if (strcmp(operation, "genkeys") == 0) {
        start = test_perf_begin();
        for (unsigned int i = 0; i < CORE_ITERATIONS; ++i)
            if (SNOVA_NAMESPACE(genkeys)(pk, sk, seed) != 0) return 0;
        test_perf_print("SNOVA core genkeys", test_perf_end(start), CORE_ITERATIONS);
    } else if (strcmp(operation, "sk_expand") == 0) {
        start = test_perf_begin();
        for (unsigned int i = 0; i < CORE_ITERATIONS; ++i)
            if (SNOVA_NAMESPACE(sk_expand)(&skx, sk) != 0) return 0;
        test_perf_print("SNOVA core sk_expand", test_perf_end(start), CORE_ITERATIONS);
    } else if (strcmp(operation, "pk_expand") == 0) {
        start = test_perf_begin();
        for (unsigned int i = 0; i < CORE_ITERATIONS; ++i)
            if (SNOVA_NAMESPACE(pk_expand)(&pkx, pk) != 0) return 0;
        test_perf_print("SNOVA core pk_expand", test_perf_end(start), CORE_ITERATIONS);
    } else if (strcmp(operation, "sign_core") == 0) {
        start = test_perf_begin();
        for (unsigned int i = 0; i < CORE_ITERATIONS; ++i)
            if (SNOVA_NAMESPACE(sign)(&skx, sig, digest, sizeof(digest), salt) != 0) return 0;
        test_perf_print("SNOVA core sign", test_perf_end(start), CORE_ITERATIONS);
    } else if (strcmp(operation, "verify_core") == 0) {
        start = test_perf_begin();
        for (unsigned int i = 0; i < CORE_ITERATIONS; ++i)
            if (SNOVA_NAMESPACE(verify)(&pkx, sig, digest, sizeof(digest)) != 0) return 0;
        test_perf_print("SNOVA core verify", test_perf_end(start), CORE_ITERATIONS);
    } else {
        return -1;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: algorithm_snova_performance keygen|sign|verify|shake|genkeys|sk_expand|pk_expand|sign_core|verify_core\n");
        return 2;
    }
    if (strcmp(argv[1], "keygen") == 0) return bench_keygen() ? 0 : 1;
    if (strcmp(argv[1], "sign") == 0 || strcmp(argv[1], "verify") == 0)
        return bench_sign_or_verify(argv[1]) ? 0 : 1;
    if (strcmp(argv[1], "shake") == 0) {
        bench_shake();
        return 0;
    }
    if (strcmp(argv[1], "genkeys") == 0 || strcmp(argv[1], "sk_expand") == 0 ||
        strcmp(argv[1], "pk_expand") == 0 || strcmp(argv[1], "sign_core") == 0 ||
        strcmp(argv[1], "verify_core") == 0) {
        int result = bench_core(argv[1]);
        return result == 1 ? 0 : 1;
    }
    fprintf(stderr, "Unknown SNOVA benchmark: %s\n", argv[1]);
    return 2;
}
