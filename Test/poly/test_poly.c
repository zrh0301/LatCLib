#include "test.h"
#include "polyvec.h"
#include "fp32.h"

static lssi_t reference_mod(lssi_t a, lssi_t q) {
    a %= q;
    if(a < 0) {
        a += q;
    }
    return a;
}

static int test_dili_fast_reduce(void) {
    lssi_t input[TEST_COEFS];
    lssi_t output[TEST_COEFS];
    uint32_t state = 0x6D5A56A9u;

    printf("Testing if Dilithium fast reduction matches %% q...");
    for(int round = 0; round < 100; round++) {
        input[0] = INT32_MIN;
        input[1] = INT32_MAX;
        input[2] = -TEST_Q_DILI;
        input[3] = TEST_Q_DILI;
        input[4] = (lssi_t)((lsdi_t)TEST_Q_DILI * 255);
        input[5] = (lssi_t)(-(lsdi_t)TEST_Q_DILI * 256);
        for(int i = 6; i < TEST_COEFS; i++) {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            input[i] = (lssi_t)state;
        }
        memcpy(output, input, sizeof(output));
        poly_reduce_mod(output, TEST_COEFS, TEST_Q_DILI);
        for(int i = 0; i < TEST_COEFS; i++) {
            lssi_t expected = reference_mod(input[i], TEST_Q_DILI);
            if(output[i] != expected || poly_reduce_mod_one(input[i], TEST_Q_DILI) != expected) {
                printf(" [FAIL] round=%d index=%d input=%d got=%d expected=%d\n", round, i, input[i], output[i], expected);
                return 0;
            }
        }
    }
    printf(" [PASS]\n");
    return 1;
}

static void fill_polynomial(poly a, int32_t q) {
    for (int i = 0; i < TEST_COEFS; ++i)
        a[i] = (int32_t)(test_random_u32() % (uint32_t)q);
}

static int test_correctness(void) {
    poly a, b, c, d;
    fp32_ctx field;
    int ok = 1;
    if (fp32_init(&field, TEST_Q_KYBER) != 0) return 0;
    fill_polynomial(a, field.q);
    fill_polynomial(b, field.q);

    /* Bounded inputs respect the signed addition/subtraction contract. */
    POLY_ADD(c, a, b, TEST_COEFS);
    POLY_SUB(d, c, b, TEST_COEFS);
    if (memcmp(a, d, sizeof(a)) != 0) {
        fprintf(stderr, "[FAIL] polynomial add/sub round trip\n");
        ok = 0;
    }
    {
        int64_t products[TEST_COEFS];
        for (int i = 0; i < TEST_COEFS; ++i)
            products[i] = (int64_t)a[i] * b[i];
        poly_mod_barrett_lsdi(c, products, field.q, TEST_QBITS, TEST_COEFS);
        for (int i = 0; i < TEST_COEFS; ++i) {
            int64_t expected = products[i] % field.q;
            if (expected > (field.q - 1) / 2) expected -= field.q;
            if (c[i] != expected) {
                fprintf(stderr, "[FAIL] centered reduction at coefficient %d\n", i);
                ok = 0;
                break;
            }
        }
    }
    printf("%s polynomial arithmetic\n", ok ? "[PASS]" : "[FAIL]");
    return ok;
}

static void montgomery_product(poly c, const poly a, const poly b,
                                const fp32_ctx *field) {
    for (int i = 0; i < TEST_COEFS; ++i)
        c[i] = fp32_mont_reduce((int64_t)a[i] * b[i], field);
}

static int test_performance(void) {
    poly a, b, c;
    fp32_ctx field;
    if (fp32_init(&field, TEST_Q_KYBER) != 0) return 0;
    fill_polynomial(a, field.q);
    fill_polynomial(b, field.q);
    BENCH_START("poly_mont_mul", 10) {
        BENCH_ITEM(montgomery_product(c, a, b, &field), 10000);
    } BENCH_FINAL(10, 10000);
    BENCH_START("poly_mod_barrett_lssi", 10) {
        BENCH_ITEM(poly_mod_barrett_lssi(c, a, field.q, TEST_QBITS, TEST_COEFS), 10000);
    } BENCH_FINAL(10, 10000);
    printf("Benchmark checksum: %d\n", c[0]);
    return 1;
}

int main(int argc, char **argv) {
    int ok;
    test_init();
    ok = test_correctness();
    ok &= test_dili_fast_reduce();
    if (argc == 2 && strcmp(argv[1], "--bench") == 0)
        ok &= test_performance();
    return ok ? 0 : 1;
}
