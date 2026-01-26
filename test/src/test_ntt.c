#include "../inc/test_common.h"

static int check_eq(lssi_t *a, lssi_t *b, int n, lssi_t q) {
    for(int i=0; i<n; i++) {
        lssi_t va = (a[i] % q + q) % q;
        lssi_t vb = (b[i] % q + q) % q;
        if (va != vb) {
            return 0;
        }
    }
    return 1;
}

void test_ntt_dilithium_table(){
    printf("\n============================================\n");
    printf("  Test Suite: using Dilithium table\n");
    printf("============================================\n");
    lp_ctx ctx;
    init_ntt_dili(&ctx);
    lssi_t q = ctx.q;
    int n = ctx.n;

    lssi_t a[256], b[256], res[256], original[256];
    poly_rand(a, n); 
    poly_rand(b, n);
    
    for(int i=0; i<n; i++) {
        a[i] = (a[i] % q + q) % q;
        b[i] = (b[i] % q + q) % q;
        original[i] = a[i];
    }

    TEST_START("NTT/INTT reversibility (Dilithium)", 1) {
        poly_ntt(a, &ctx);
        poly_invntt(a, &ctx);
        TEST_ITEM(check_eq(a, original, n, q), end_dili_correct);
    } TEST_FINAL
    end_dili_correct:

    printf("\n[Benchmark] Dilithium NTT Performance:\n");

    BENCH_START("INIT", 10) {
        BENCH_ITEM(init_ntt_dili(&ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("NTT (8-layer)", 10) {
        BENCH_ITEM(poly_ntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("INTT (8-layer)", 10) {
        BENCH_ITEM(poly_invntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("Scalar Mul", 10) {
        BENCH_ITEM(poly_base_pointwise(res, a, b, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

void test_dilithium() {
    printf("\n============================================\n");
    printf("  Test Suite: Dilithium (Full NTT Mode)\n");
    printf("============================================\n");

    lssi_t q = TEST_Q_DILI;
    lssi_t root = 1753; 
    int n = 256;
    lp_ctx ctx;
    init_ntt(&ctx, n, root, NTT_FULL_MODE);

    lssi_t a[256], b[256], res[256], original[256];
    poly_rand(a, n); 
    poly_rand(b, n);
    
    for(int i=0; i<n; i++) {
        a[i] = (a[i] % q + q) % q;
        b[i] = (b[i] % q + q) % q;
        original[i] = a[i];
    }

    TEST_START("NTT/INTT reversibility (Dilithium)", 1) {
        poly_ntt(a, &ctx);
        poly_invntt(a, &ctx);
        TEST_ITEM(check_eq(a, original, n, q), end_dili_correct);
    } TEST_FINAL
    end_dili_correct:

    printf("\n[Benchmark] Dilithium NTT Performance:\n");

    BENCH_START("INIT", 10) {
        BENCH_ITEM(init_ntt(&ctx, n, root, NTT_FULL_MODE), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("NTT (8-layer)", 10) {
        BENCH_ITEM(poly_ntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("INTT (8-layer)", 10) {
        BENCH_ITEM(poly_invntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("Scalar Mul", 10) {
        BENCH_ITEM(poly_base_pointwise(res, a, b, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

void test_ntt_kyber_table() {
    printf("\n============================================\n");
    printf("  Test Suite: using Kyber table\n");
    printf("============================================\n");

    
    lp_ctx ctx;
    init_ntt_kyber(&ctx);
    lssi_t q = ctx.q;
    int n = ctx.n;

    lssi_t a[256], b[256], res[256], original[256];
    poly_rand(a, n);
    poly_rand(b, n);
    
    for(int i=0; i<n; i++) {
        a[i] = (a[i] % q + q) % q;
        b[i] = (b[i] % q + q) % q;
        original[i] = a[i];
    }

    TEST_START("NTT/INTT reversibility (Kyber)", 1) {
        poly_ntt(a, &ctx);
        poly_invntt(a, &ctx);
        TEST_ITEM(check_eq(a, original, n, q), end_kyber_correct);
    } TEST_FINAL
    end_kyber_correct:

    printf("\n[Benchmark] Kyber NTT Performance:\n");

    BENCH_START("INIT", 10) {
        BENCH_ITEM(init_ntt_kyber(&ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    BENCH_START("NTT (7-layer)", 10) {
        BENCH_ITEM(poly_ntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("INTT (7-layer)", 10) {
        BENCH_ITEM(poly_invntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("Basecase Mul", 10) {
        BENCH_ITEM(poly_base_pointwise(res, a, b, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

void test_kyber() {
    printf("\n============================================\n");
    printf("  Test Suite: Kyber (Incomplete NTT Mode)\n");
    printf("============================================\n");

    lssi_t q = TEST_Q_KYBER;
    lssi_t root = 17; 
    int n = 256;
    lp_ctx ctx;
    init_ntt(&ctx, n, root, NTT_INCOMPLETE_MODE);

    lssi_t a[256], b[256], res[256], original[256];
    poly_rand(a, n);
    poly_rand(b, n);
    
    for(int i=0; i<n; i++) {
        a[i] = (a[i] % q + q) % q;
        b[i] = (b[i] % q + q) % q;
        original[i] = a[i];
    }

    TEST_START("NTT/INTT reversibility (Kyber)", 1) {
        poly_ntt(a, &ctx);
        poly_invntt(a, &ctx);
        TEST_ITEM(check_eq(a, original, n, q), end_kyber_correct);
    } TEST_FINAL
    end_kyber_correct:

    printf("\n[Benchmark] Kyber NTT Performance:\n");

    BENCH_START("INIT", 10) {
        BENCH_ITEM(init_ntt(&ctx, n, root, NTT_INCOMPLETE_MODE), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("NTT (7-layer)", 10) {
        BENCH_ITEM(poly_ntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("INTT (7-layer)", 10) {
        BENCH_ITEM(poly_invntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    BENCH_START("Basecase Mul", 10) {
        BENCH_ITEM(poly_base_pointwise(res, a, b, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

void test_no_table() {
    printf("\n========================================================================================\n");
    printf("  Test Suite: No Precomputed Table\n");
    printf("========================================================================================\n");
    
    test_dilithium();
    test_kyber();
}

void test_using_tables() {
    printf("\n========================================================================================\n");
    printf("  Test Suite: Using Precomputed Tables\n");
    printf("========================================================================================\n");

    test_ntt_dilithium_table();
    test_ntt_kyber_table();
}

int main() {
    test_init();
    test_no_table();
    test_using_tables();
    return 0;
}