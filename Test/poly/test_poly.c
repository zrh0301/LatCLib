#include "test.h"
#include "polyvec.h"

typedef struct {
    const char *name;
    int (*run)(void);
} test_case;

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

static int test_correctness() {
    int ok = 1;

    printf("=== Running Correctness Tests ===\n");
    
    lssi_t a[TEST_COEFS], b[TEST_COEFS], c[TEST_COEFS], d[TEST_COEFS];
    lssi_t qc;
    lssi_t q = TEST_Q_KYBER;
    int ret;

    TEST_START("poly_rand generates random polynomials", 1) {
        poly_rand(a, TEST_COEFS);
        poly_rand(b, TEST_COEFS);
        TEST_ITEM(1, end_rand);
    } TEST_FINAL
    end_rand:

    TEST_START("poly_add/sub: (a+b)-b == a", 1) {
        POLY_ADD(c, a, b, TEST_COEFS);
        POLY_SUB(d, c, b, TEST_COEFS);
        TEST_ITEM(memcmp(a, d, TEST_COEFS * sizeof(lssi_t)) == 0, end_addsub);
    } TEST_FINAL
    end_addsub:

    TEST_START("poly_mont_init initializes correctly", 1) {
        ret = POLY_MONT_INIT(&qc, q);
        TEST_ITEM(ret == 0, end_mont);
        printf("(qc = %d) ", qc);
    } TEST_FINAL
    end_mont:

    TEST_START("poly_mod_barrett_lsdi reduces to valid range", 1) {
        lsdi_t big_a[TEST_COEFS];
        int i;
        for(i=0; i<TEST_COEFS; i++) big_a[i] = (lsdi_t)a[i] * b[i];
        
        poly_mod_barrett_lsdi(c, big_a, q, TEST_QBITS, TEST_COEFS);
        
        TEST_ITEM(check_poly_range(c, TEST_COEFS, -(q-1)/2, (q-1)/2), end_barrett);
    } TEST_FINAL
    end_barrett:

    return ok;
}

static int test_performance() {
    printf("\n=== Running Performance Tests ===\n");
    
    lssi_t a[TEST_COEFS], b[TEST_COEFS], c[TEST_COEFS];
    lssi_t q = TEST_Q_KYBER;
    lssi_t qc;
    POLY_MONT_INIT(&qc, q);
    poly_rand(a, TEST_COEFS);
    poly_rand(b, TEST_COEFS);

    BENCH_START("poly_mont_mul", 10) {
        BENCH_ITEM(POLY_MONT_MUL(c, a, b, q, qc, TEST_COEFS), 100000);
    } BENCH_FINAL(10, 100000);

    BENCH_START("poly_mod_barrett_lssi", 10) {
        BENCH_ITEM(poly_mod_barrett_lssi(c, a, q, TEST_QBITS, TEST_COEFS), 100000);
    } BENCH_FINAL(10, 100000);

    return 1;
}

int main() {
    const test_case poly_cases[] = {
        {"poly correctness", test_correctness},
        {"Dilithium fast reduction", test_dili_fast_reduce},
        {"poly performance", test_performance},
    };
    int ok = 1;

    test_init();

    for(size_t i = 0; i < sizeof(poly_cases) / sizeof(poly_cases[0]); i++) {
        ok &= poly_cases[i].run();
    }

    return ok ? 0 : 1;
}
