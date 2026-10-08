#include "test.h"
#include "dili_sample.h"
#include "kyber_sample.h"
#include "kyber_pack.h"

/* Range/count checks for the existing sampling primitives. These do not
 * assert standards conformance for the legacy gamma helper functions. */
static int sampled_range(const char *name, const poly output, unsigned int count,
                          int32_t low, int32_t high) {
    if (count != TEST_COEFS || !check_poly_range((int32_t *)output, TEST_COEFS, low, high)) {
        fprintf(stderr, "[FAIL] %s count/range: count=%u, expected=%d\n", name, count, TEST_COEFS);
        return 0;
    }
    printf("[PASS] %s count/range\n", name);
    return 1;
}

static int test_sampling_primitives(void) {
    lp_ctx ctx;
    uint8_t bytes[2048];
    poly output;
    unsigned int count;
    int ok = 1;
    INIT_NTT_DILI(&ctx);
    test_random_bytes(bytes, sizeof(bytes));
    count = sample_uniform(output, TEST_COEFS, bytes, sizeof(bytes), &ctx);
    ok &= sampled_range("uniform", output, count, 0, DILI_Q - 1);
    for (int eta = 2; eta <= 4; eta += 2) {
        count = sample_eta(output, TEST_COEFS, bytes, sizeof(bytes), &ctx, eta);
        ok &= sampled_range("eta", output, count, -eta, eta);
    }
    count = sample_CBD(output, TEST_COEFS, bytes, sizeof(bytes), &ctx);
    ok &= sampled_range("CBD eta=2", output, count, -2, 2);
    count = sample_gamma1(output, TEST_COEFS, bytes, sizeof(bytes), &ctx);
    ok &= sampled_range("legacy gamma1 helper", output, count, -(1 << 17), (1 << 17) - 1);
    count = sample_gamma2(output, TEST_COEFS, bytes, sizeof(bytes), &ctx);
    ok &= sampled_range("legacy gamma2 helper", output, count, 0, 15);

    memset(bytes, 0, sizeof(bytes));
    count = sample_uniform(output, TEST_COEFS, bytes, sizeof(bytes), &ctx);
    ok &= sampled_range("uniform all-zero input", output, count, 0, 0);
    count = sample_CBD(output, TEST_COEFS, bytes, sizeof(bytes), &ctx);
    ok &= sampled_range("CBD all-zero input", output, count, 0, 0);
    memset(bytes, 0xff, sizeof(bytes));
    if (sample_uniform(output, TEST_COEFS, bytes, sizeof(bytes), &ctx) != 0 ||
        sample_eta(output, TEST_COEFS, bytes, sizeof(bytes), &ctx, 2) != 0 ||
        sample_eta(output, TEST_COEFS, bytes, sizeof(bytes), &ctx, 4) != 0) {
        fprintf(stderr, "[FAIL] rejection samplers accepted out-of-range input\n");
        ok = 0;
    }
    output[0] = 12345;
    if (sample_uniform(output, TEST_COEFS, bytes, 0, &ctx) != 0 ||
        sample_eta(output, TEST_COEFS, bytes, 0, &ctx, 2) != 0 ||
        sample_CBD(output, TEST_COEFS, bytes, 0, &ctx) != 0 ||
        sample_gamma1(output, TEST_COEFS, bytes, 0, &ctx) != 0 ||
        sample_gamma2(output, TEST_COEFS, bytes, 0, &ctx) != 0 || output[0] != 12345) {
        fprintf(stderr, "[FAIL] empty sampling input changed output\n");
        ok = 0;
    }
    return ok;
}

static int test_sampling_short_buffers(void) {
    lp_ctx ctx;
    const uint8_t short_uniform[2] = {0, 0};
    const uint8_t short_gamma[8] = {0};
    const uint8_t full_gamma[9] = {0};
    int32_t output[5] = {12345, 12345, 12345, 12345, 12345};
    INIT_NTT_DILI(&ctx);
    if (sample_uniform(output, 4, short_uniform, sizeof(short_uniform), &ctx) != 0 ||
        output[0] != 12345 ||
        sample_gamma1(output, 4, short_gamma, sizeof(short_gamma), &ctx) != 0 ||
        output[0] != 12345) {
        fprintf(stderr, "[FAIL] short sampler input was consumed\n");
        return 0;
    }
    if (sample_gamma1(output, 4, full_gamma, sizeof(full_gamma), &ctx) != 4 ||
        output[4] != 12345) {
        fprintf(stderr, "[FAIL] complete gamma group count/output bound\n");
        return 0;
    }
    for (int i = 0; i < 4; ++i) {
        if (output[i] != -(1 << 17)) {
            fprintf(stderr, "[FAIL] complete zero gamma group value\n");
            return 0;
        }
    }
    puts("[PASS] sampler short-buffer boundaries");
    return 1;
}

static void fill_kyber_bytes(uint8_t *buf, size_t len, uint8_t seed) {
    for(size_t i = 0; i < len; i++) {
        buf[i] = (uint8_t)(seed + 29u * (uint8_t)i + (uint8_t)(i >> 2));
    }
}

static int kyber_poly_equal(const poly a, const poly b) {
    for(int i = 0; i < KYBER_N; i++) {
        if(a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static int kyber_centered_value(lssi_t x) {
    x = kyber_reduce(x);
    if(x > KYBER_Q / 2) {
        x -= KYBER_Q;
    }
    return x;
}

int test_kyber_sampling_functions(void) {
    printf("\n=== Test Kyber SampleNTT / CBD / PRF noise ===\n");

    lp_ctx ctx;
    uint8_t seed[34];
    uint8_t seed2[34];
    uint8_t sigma[KYBER_SYMBYTES];
    uint8_t zero_buf[64 * 3];
    uint8_t pos_eta2_buf[64 * 2];
    uint8_t neg_eta2_buf[64 * 2];
    poly a1, a2, a3;
    poly e1, e2, e3;
    int ok = 1;

    INIT_NTT_KYBER(&ctx);
    fill_kyber_bytes(seed, sizeof(seed), 0x19);
    memcpy(seed2, seed, sizeof(seed));
    seed2[33] ^= 0x01;
    fill_kyber_bytes(sigma, sizeof(sigma), 0xA7);

    kyber_sample_ntt(a1, seed, &ctx);
    kyber_sample_ntt(a2, seed, &ctx);
    kyber_sample_ntt(a3, seed2, &ctx);

    if(!kyber_poly_equal(a1, a2)) {
        printf("  [FAIL] Kyber SampleNTT is not deterministic for the same seed\n");
        ok = 0;
    }
    if(kyber_poly_equal(a1, a3)) {
        printf("  [FAIL] Kyber SampleNTT did not change after seed domain changed\n");
        ok = 0;
    }
    for(int i = 0; i < KYBER_N; i++) {
        if(a1[i] < 0 || a1[i] >= KYBER_Q) {
            printf("  [FAIL] Kyber SampleNTT out of range at %d: %d\n", i, a1[i]);
            ok = 0;
            break;
        }
    }

    memset(zero_buf, 0, sizeof(zero_buf));
    kyber_sample_cbd(e1, zero_buf, 2);
    kyber_sample_cbd(e2, zero_buf, 3);
    for(int i = 0; i < KYBER_N; i++) {
        if(e1[i] != 0 || e2[i] != 0) {
            printf("  [FAIL] Kyber CBD zero-buffer test failed at %d\n", i);
            ok = 0;
            break;
        }
    }

    memset(pos_eta2_buf, 0x33, sizeof(pos_eta2_buf));
    memset(neg_eta2_buf, 0xCC, sizeof(neg_eta2_buf));
    kyber_sample_cbd(e1, pos_eta2_buf, 2);
    kyber_sample_cbd(e2, neg_eta2_buf, 2);
    for(int i = 0; i < KYBER_N; i++) {
        if(e1[i] != 2 || e2[i] != KYBER_Q - 2) {
            printf("  [FAIL] Kyber CBD eta=2 exact-pattern test failed at %d: %d %d\n",
                   i, e1[i], e2[i]);
            ok = 0;
            break;
        }
    }

    kyber_sample_noise(e1, sigma, 7, 3);
    kyber_sample_noise(e2, sigma, 7, 3);
    kyber_sample_noise(e3, sigma, 8, 3);
    if(!kyber_poly_equal(e1, e2)) {
        printf("  [FAIL] Kyber noise sampling is not deterministic for same sigma/nonce\n");
        ok = 0;
    }
    if(kyber_poly_equal(e1, e3)) {
        printf("  [FAIL] Kyber noise sampling did not change after nonce changed\n");
        ok = 0;
    }
    for(int i = 0; i < KYBER_N; i++) {
        int centered = kyber_centered_value(e1[i]);
        if(centered < -3 || centered > 3) {
            printf("  [FAIL] Kyber eta=3 noise out of centered range at %d: %d\n",
                   i, centered);
            ok = 0;
            break;
        }
    }

    printf("  %s Kyber sampling unit tests\n", ok ? "[PASS]" : "[FAIL]");
    return ok;
}

int main(void) {
    int ok;
    test_init();
    ok = test_sampling_primitives();
    ok &= test_sampling_short_buffers();
    ok &= test_kyber_sampling_functions();
    return ok ? 0 : 1;
}
