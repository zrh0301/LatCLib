#include "../inc/test_common.h"

void test_correctness() {
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
        poly_add(c, a, b, TEST_COEFS);
        poly_sub(d, c, b, TEST_COEFS);
        TEST_ITEM(memcmp(a, d, TEST_COEFS * sizeof(lssi_t)) == 0, end_addsub);
    } TEST_FINAL
    end_addsub:

    TEST_START("poly_mont_init initializes correctly", 1) {
        ret = poly_mont_init(&qc, q);
        TEST_ITEM(ret == 0, end_mont);
        printf("(qc = %d) ", qc);
    } TEST_FINAL
    end_mont:

    TEST_START("poly_mod_barrett_lsdi reduces to valid range", 1) {
        lsdi_t big_a[TEST_COEFS];
        for(int i=0; i<TEST_COEFS; i++) big_a[i] = (lsdi_t)a[i] * b[i];
        
        poly_mod_barrett_lsdi(c, big_a, q, TEST_QBITS, TEST_COEFS);
        
        TEST_ITEM(check_poly_range(c, TEST_COEFS, -(q-1)/2, (q-1)/2), end_barrett);
    } TEST_FINAL
    end_barrett:
}

void test_performance() {
    printf("\n=== Running Performance Tests ===\n");
    
    lssi_t a[TEST_COEFS], b[TEST_COEFS], c[TEST_COEFS];
    lssi_t q = TEST_Q_KYBER;
    lssi_t qc;
    poly_mont_init(&qc, q);
    poly_rand(a, TEST_COEFS);
    poly_rand(b, TEST_COEFS);

    BENCH_START("poly_mont_mul", 10) {
        BENCH_ITEM(poly_mont_mul(c, a, b, q, qc, TEST_COEFS), 100000);
    } BENCH_FINAL(10, 100000);

    BENCH_START("poly_mod_barrett_lssi", 10) {
        BENCH_ITEM(poly_mod_barrett_lssi(c, a, q, TEST_QBITS, TEST_COEFS), 100000);
    } BENCH_FINAL(10, 100000);
}

int main() {
    test_init();

    test_correctness();
    test_performance();

    return 0;
}