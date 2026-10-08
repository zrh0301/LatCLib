#include "rand.h"
#include "poly_sample.h"
#include "dili.h"
#include "kyber.h"
#include <stdio.h>
#include <string.h>

#ifdef PQC_TEST_WRAP_GETRANDOM
#include <errno.h>
#include <sys/random.h>

enum { OS_REAL, OS_SHORT_READS, OS_PARTIAL_FAILURE, OS_ZERO_READ };
static int os_test_mode;
static size_t os_test_calls;
static size_t os_test_written;

ssize_t __real_getrandom(void *buf, size_t buflen, unsigned int flags);

ssize_t __wrap_getrandom(void *buf, size_t buflen, unsigned int flags) {
    if(os_test_mode == OS_REAL) {
        return __real_getrandom(buf, buflen, flags);
    }
    ++os_test_calls;
    if(os_test_mode == OS_ZERO_READ) {
        return 0;
    }
    if(os_test_mode == OS_SHORT_READS && os_test_calls == 1) {
        errno = EINTR;
        return -1;
    }
    if(os_test_mode == OS_PARTIAL_FAILURE && os_test_calls > 1) {
        errno = EIO;
        return -1;
    }
    size_t len = buflen < 7 ? buflen : 7;
    for(size_t i = 0; i < len; ++i) {
        ((uint8_t *)buf)[i] = (uint8_t)(os_test_written++);
    }
    return (ssize_t)len;
}

static void set_os_test_mode(int mode) {
    os_test_mode = mode;
    os_test_calls = 0;
    os_test_written = 0;
}
#endif

#define CHECK(condition) do { \
    if(!(condition)) { \
        fprintf(stderr, "RNG check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while(0)

/* Determinism belongs to the test, never the production entropy source. */
typedef struct {
    uint8_t seed[32];
    size_t calls;
    size_t requested;
    int fail;
    int partial;
} fixed_rng;

static int fixed_random(void *arg, uint8_t *out, size_t len) {
    fixed_rng *rng = arg;
    rng->calls++;
    rng->requested = len;
    if(rng->fail) {
        if(rng->partial && len != 0) {
            memset(out, 0x3c, len / 2 + 1);
        }
        return 7; /* A nonzero callback result is consistently mapped to -1. */
    }
    if(len > sizeof(rng->seed)) {
        return -1;
    }
    memcpy(out, rng->seed, len);
    return 0;
}

static int all_bytes(const void *data, size_t len, uint8_t value) {
    const uint8_t *p = data;
    for(size_t i = 0; i < len; ++i) {
        if(p[i] != value) {
            return 0;
        }
    }
    return 1;
}

static fixed_rng make_rng(void) {
    fixed_rng rng = {0};
    for(size_t i = 0; i < sizeof(rng.seed); ++i) {
        rng.seed[i] = (uint8_t)(3 * i + 1);
    }
    return rng;
}

#ifdef PQC_TEST_WRAP_GETRANDOM
static int check_system_failures(void) {
    const dili_param param = {
        0, 4, 4, 2, 39, 78, 80, 1 << 17, (DILI_Q - 1) / 88, 192, 32
    };
    uint8_t guarded[1026];
    uint8_t signature[2420];
    uint8_t ek[KYBER_MAX_EK_BYTES], dk[KYBER_MAX_DK_PKE_BYTES];
    dili_keypair_t keypair;
    kyber_param kyber;

    set_os_test_mode(OS_SHORT_READS);
    memset(guarded, 0xa5, sizeof(guarded));
    CHECK(pqc_system_random(NULL, guarded + 1, 1024) == 0);
    CHECK(os_test_calls > 2 && os_test_written == 1024);
    CHECK(guarded[0] == 0xa5 && guarded[1025] == 0xa5);
    for(size_t i = 0; i < 1024; ++i) {
        CHECK(guarded[i + 1] == (uint8_t)i);
    }

    for(int mode = OS_PARTIAL_FAILURE; mode <= OS_ZERO_READ; ++mode) {
        memset(guarded, 0xa5, sizeof(guarded));
        set_os_test_mode(mode);
        CHECK(pqc_system_random(NULL, guarded + 1, 1024) == -1);
        CHECK(all_bytes(guarded + 1, 1024, 0));
        CHECK(guarded[0] == 0xa5 && guarded[1025] == 0xa5);
    }

    /* NULL callback must propagate system failures at each algorithm boundary. */
    memset(&keypair, 0xa5, sizeof(keypair));
    set_os_test_mode(OS_PARTIAL_FAILURE);
    CHECK(dili_keygen_with_rng(&keypair, param, NULL, NULL) == -1);
    CHECK(all_bytes(&keypair, sizeof(keypair), 0xa5));
    memset(signature, 0xa5, sizeof(signature));
    set_os_test_mode(OS_PARTIAL_FAILURE);
    CHECK(dili_sign(signature, NULL, 0, &keypair.sk, param) == -1);
    CHECK(all_bytes(signature, sizeof(signature), 0xa5));
    kyber_init_param(&kyber, KYBER_MODE_2);
    memset(ek, 0xa5, sizeof(ek));
    memset(dk, 0xa5, sizeof(dk));
    set_os_test_mode(OS_PARTIAL_FAILURE);
    CHECK(kyber_pke_keygen_with_rng(ek, dk, kyber, NULL, NULL) == -1);
    CHECK(all_bytes(ek, sizeof(ek), 0xa5) && all_bytes(dk, sizeof(dk), 0xa5));
    set_os_test_mode(OS_REAL);
    return 0;
}
#endif

static int check_random_source(void) {
    uint8_t output[32];
    uint8_t guarded[1026];
    fixed_rng rng = make_rng();

    CHECK(pqc_randombytes(output, sizeof(output), fixed_random, &rng) == 0);
    CHECK(rng.calls == 1 && rng.requested == sizeof(output));
    CHECK(memcmp(output, rng.seed, sizeof(output)) == 0);
    for(int partial = 0; partial <= 1; ++partial) {
        rng.fail = 1;
        rng.partial = partial;
        memset(output, 0xa5, sizeof(output));
        CHECK(pqc_randombytes(output, sizeof(output), fixed_random, &rng) == -1);
        CHECK(all_bytes(output, sizeof(output), 0));
    }
    rng.calls = 0;
    CHECK(pqc_randombytes(NULL, 0, fixed_random, &rng) == 0);
    CHECK(pqc_randombytes(NULL, 1, fixed_random, &rng) == -1);
    CHECK(rng.calls == 0);
    CHECK(pqc_system_random(NULL, NULL, 0) == 0);
    CHECK(pqc_system_random(NULL, NULL, 1) == -1);

    /* Exercise more than one OS read and verify exact output boundaries. */
    memset(guarded, 0xa5, sizeof(guarded));
    CHECK(pqc_randombytes(guarded + 1, 1024, NULL, NULL) == 0);
    CHECK(guarded[0] == 0xa5 && guarded[1025] == 0xa5);
    CHECK(pqc_system_random(NULL, output, sizeof(output)) == 0);
    rand_bytes(NULL, 0);
    rand_bits(NULL, 0);
    for(int bits = 1; bits <= 17; ++bits) {
        size_t bytes = (size_t)(bits + 7) / 8;
        memset(output, 0xa5, sizeof(output));
        rand_bits(output, bits);
        if(bits % 8 != 0) {
            CHECK((output[bytes - 1] >> (bits % 8)) == 0);
        }
        CHECK(output[bytes] == 0xa5);
    }
    return 0;
}

static int check_dilithium(void) {
    const dili_param params[] = {
        {0, 4, 4, 2, 39, 78, 80, 1 << 17, (DILI_Q - 1) / 88, 192, 32},
        {1, 6, 5, 4, 49, 196, 55, 1 << 19, (DILI_Q - 1) / 32, 128, 48},
        {2, 8, 7, 2, 60, 120, 75, 1 << 19, (DILI_Q - 1) / 32, 128, 64}
    };
    const uint8_t message[] = "injected RNG regression";
    dili_keypair_t seeded, injected;
    uint8_t signature[4627], expected[4627];
    fixed_rng rng = make_rng();

    for(size_t mode = 0; mode < sizeof(params) / sizeof(params[0]); ++mode) {
        dili_param param = params[mode];
        size_t sig_len = (size_t)param.c_tilde_bytes +
                         (size_t)param.l * DILI_POLYZ_PACKEDBYTES(param.gamma1) +
                         (size_t)param.omega + (size_t)param.k;
        memset(&seeded, 0, sizeof(seeded));
        memset(&injected, 0, sizeof(injected));
        rng.fail = 0;
        rng.calls = 0;
        CHECK(dili_keygen_with_seed(&seeded, param, rng.seed) == 0);
        CHECK(dili_keygen_with_rng(&injected, param, fixed_random, &rng) == 0);
        CHECK(rng.calls == 1 && rng.requested == 32);
        CHECK(memcmp(&seeded, &injected, sizeof(seeded)) == 0);
        CHECK(dili_sign_with_rnd(expected, message, sizeof(message) - 1,
                                 &seeded.sk, param, rng.seed) == 0);
        CHECK(dili_sign_with_rng(signature, message, sizeof(message) - 1,
                                 &seeded.sk, param, fixed_random, &rng) == 0);
        CHECK(memcmp(signature, expected, sig_len) == 0);
        CHECK(dili_verify(signature, message, sizeof(message) - 1, &seeded.pk, param) == 0);

        for(int partial = 0; partial <= 1; ++partial) {
            rng.fail = 1;
            rng.partial = partial;
            memset(&injected, 0xa5, sizeof(injected));
            CHECK(dili_keygen_with_rng(&injected, param, fixed_random, &rng) == -1);
            CHECK(all_bytes(&injected, sizeof(injected), 0xa5));
            memset(signature, 0xa5, sizeof(signature));
            CHECK(dili_sign_with_rng(signature, message, sizeof(message) - 1,
                                     &seeded.sk, param, fixed_random, &rng) == -1);
            CHECK(all_bytes(signature, sizeof(signature), 0xa5));
        }
    }

    rng.calls = 0;
    CHECK(dili_keygen_with_seed(NULL, params[0], rng.seed) == -1);
    CHECK(dili_keygen_with_seed(&injected, params[0], NULL) == -1);
    CHECK(dili_keygen_with_rng(NULL, params[0], fixed_random, &rng) == -1);
    CHECK(dili_sign_with_rng(NULL, message, sizeof(message), &seeded.sk,
                             params[0], fixed_random, &rng) == -1);
    CHECK(dili_sign_with_rng(signature, NULL, 1, &seeded.sk,
                             params[0], fixed_random, &rng) == -1);
    CHECK(dili_sign_with_rng(signature, message, -1, &seeded.sk,
                             params[0], fixed_random, &rng) == -1);
    CHECK(dili_sign_with_rng(signature, message, sizeof(message), NULL,
                             params[0], fixed_random, &rng) == -1);
    CHECK(dili_sign_with_rnd(signature, message, sizeof(message), &seeded.sk,
                             params[0], NULL) == -1);
    CHECK(rng.calls == 0);
    /* Exercise the default OS RNG at the algorithm boundary as well. */
    CHECK(dili_keygen_with_rng(&injected, params[0], NULL, NULL) == 0);
    CHECK(dili_sign(signature, message, sizeof(message), &injected.sk, params[0]) == 0);
    CHECK(dili_verify(signature, message, sizeof(message), &injected.pk, params[0]) == 0);
    return 0;
}

static int check_kyber(void) {
    const int modes[] = {KYBER_MODE_2, KYBER_MODE_3, KYBER_MODE_5};
    uint8_t ek[KYBER_MAX_EK_BYTES], dk[KYBER_MAX_DK_PKE_BYTES];
    uint8_t expected_ek[KYBER_MAX_EK_BYTES], expected_dk[KYBER_MAX_DK_PKE_BYTES];
    fixed_rng rng = make_rng();
    kyber_param param;

    for(size_t mode = 0; mode < sizeof(modes) / sizeof(modes[0]); ++mode) {
        kyber_init_param(&param, modes[mode]);
        rng.fail = 0;
        rng.calls = 0;
        CHECK(kyber_pke_keygen(expected_ek, expected_dk, rng.seed, param) == 0);
        CHECK(kyber_pke_keygen_with_rng(ek, dk, param, fixed_random, &rng) == 0);
        CHECK(rng.calls == 1 && rng.requested == 32);
        CHECK(memcmp(ek, expected_ek, (size_t)param.ek_bytes) == 0);
        CHECK(memcmp(dk, expected_dk, (size_t)param.dk_pke_bytes) == 0);
        for(int partial = 0; partial <= 1; ++partial) {
            rng.fail = 1;
            rng.partial = partial;
            memset(ek, 0xa5, sizeof(ek));
            memset(dk, 0xa5, sizeof(dk));
            CHECK(kyber_pke_keygen_with_rng(ek, dk, param, fixed_random, &rng) == -1);
            CHECK(all_bytes(ek, sizeof(ek), 0xa5) && all_bytes(dk, sizeof(dk), 0xa5));
        }
    }
    rng.calls = 0;
    CHECK(kyber_pke_keygen(NULL, dk, rng.seed, param) == -1);
    CHECK(kyber_pke_keygen(ek, NULL, rng.seed, param) == -1);
    CHECK(kyber_pke_keygen(ek, dk, NULL, param) == -1);
    CHECK(kyber_pke_keygen_with_rng(NULL, dk, param, fixed_random, &rng) == -1);
    CHECK(kyber_pke_keygen_with_rng(ek, NULL, param, fixed_random, &rng) == -1);
    CHECK(rng.calls == 0);
    CHECK(kyber_pke_keygen_with_rng(ek, dk, param, NULL, NULL) == 0);
    return 0;
}

static int patterned_random(void *arg, uint8_t *out, size_t len) {
    (void)arg;
    for (size_t i = 0; i < len; ++i) out[i] = (uint8_t)(i * 17 + 0x83);
    return 0;
}

static int check_poly_sampling_rng(void) {
    poly a, b;
    lp_ctx plan;
    fixed_rng rng = make_rng();
    INIT_NTT_KYBER(&plan);
    CHECK(poly_rand_with_rng(a, 256, patterned_random, NULL) == 0);
    CHECK(poly_rand_with_rng(b, 256, patterned_random, NULL) == 0);
    CHECK(memcmp(a, b, sizeof(a)) == 0);
    CHECK(poly_uniform_CBD_with_rng(&a, &plan, patterned_random, NULL) == 0);
    for (int i = 0; i < 256; ++i) CHECK(a[i] >= -2 && a[i] <= 2);
    memset(a, 0xa5, sizeof(a));
    rng.fail = 1;
    rng.partial = 1;
    CHECK(poly_rand_with_rng(a, 256, fixed_random, &rng) == -1);
    CHECK(all_bytes(a, sizeof(a), 0xa5));
    CHECK(poly_uniform_CBD_with_rng(&a, &plan, fixed_random, &rng) == -1);
    CHECK(all_bytes(a, sizeof(a), 0xa5));
    rng.calls = 0;
    CHECK(poly_rand_with_rng(a, 257, fixed_random, &rng) == -1);
    CHECK(poly_rand_with_rng(a, -1, fixed_random, &rng) == -1);
    CHECK(poly_uniform_CBD_with_rng(&a, NULL, fixed_random, &rng) == -1);
    CHECK(rng.calls == 0);
    return 0;
}

int main(void) {
#ifdef PQC_TEST_WRAP_GETRANDOM
    if(check_system_failures()) {
        return 1;
    }
#endif
    if(check_random_source() || check_dilithium() || check_kyber() ||
       check_poly_sampling_rng()) {
        return 1;
    }
    puts("RNG injection, error propagation, and seeded reproducibility: PASS");
    return 0;
}
