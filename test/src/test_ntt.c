#include "../inc/test_common.h"

static int check_eq(lssi_t *a, lssi_t *b, int n, lssi_t q) {
    int i;
    for(i=0; i<n; i++) {
        lssi_t va = (a[i] % q + q) % q;
        lssi_t vb = (b[i] % q + q) % q;
        if (va != vb) {
            return 0;
        }
    }
    return 1;
}

void test_ntt_dilithium(){
    printf("\n============================================\n");
    printf("  Test Suite: NTT-8\n");
    printf("============================================\n");
    lp_ctx ctx;
    int i;
    init_ntt_dili(&ctx);
    lssi_t q = ctx.q;
    int n = ctx.n;

    lssi_t a[256], b[256], res[256], original[256];
    poly_rand(a, n); 
    poly_rand(b, n);
    
    for(i=0; i<n; i++) {
        a[i] = (a[i] % q + q) % q;
        b[i] = (b[i] % q + q) % q;
        original[i] = a[i];
    }

    TEST_START("NTT/INTT reversibility (Dilithium)", 1) {
        poly_ntt(a, &ctx);
        poly_invntt(a, &ctx);
        for(i=0; i<n; i++) {
            a[i] = FROM_MONTGOMERY(a[i], q, ctx.qc);
        }
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

void test_using_tables() {
    printf("\n========================================================================================\n");
    printf("  Test Suite: Using Precomputed Tables\n");
    printf("========================================================================================\n");

    test_ntt_dilithium();
    //test_ntt_kyber_table();
}

static void poly_mult_ntt(lssi_t* c, lssi_t* a_tmp, lssi_t* b_tmp, lp_ctx* Q) {
    poly_ntt(a_tmp, Q);
    poly_ntt(b_tmp, Q);
    poly_pointwise(c, a_tmp, b_tmp, Q);
    poly_invntt(c, Q);
}

void test_poly_mul_ntt() {
    printf("\n============================================\n");
    printf("  Test Suite: Polynomial Multiplication via NTT-8\n");
    printf("============================================\n");
    
    lp_ctx ctx;
    int i;
    init_ntt_dili(&ctx);
    lssi_t q = ctx.q;
    int n = ctx.n;
    
    lssi_t a[256], b[256], c[256];
    lssi_t a_tmp[256], b_tmp[256];
    TEST_START("Zero polynomial multiplication", 1) {
        POLY_ZERO(a, n);
        poly_rand(b, n);
        for(i=0; i<n; i++) b[i] = (b[i] % q + q) % q;
        
        POLY_COPY(a_tmp, a, n);
        POLY_COPY(b_tmp, b, n);
        poly_mult_ntt(c, a_tmp, b_tmp, &ctx);
        
        int is_zero = 1;
        for(i=0; i<n; i++) {
            lssi_t val = (c[i] % q + q) % q;
            if(val != 0) {
                is_zero = 0;
                break;
            }
        }
        TEST_ITEM(is_zero, end_zero);
    } TEST_FINAL
    end_zero:
    
    TEST_START("Identity polynomial multiplication", 1) {
        POLY_ZERO(a, n);
        a[0] = 5; a[1] = 10; a[2] = 15;
        
        POLY_ZERO(b, n);
        b[0] = 1;
        
        POLY_COPY(a_tmp, a, n);
        POLY_COPY(b_tmp, b, n);
        poly_mult_ntt(c, a_tmp, b_tmp, &ctx);
        
        TEST_ITEM(check_eq(a, c, n, q), end_identity);
    } TEST_FINAL
    end_identity:
    
    TEST_START("Commutativity (a*b = b*a)", 1) {
        lssi_t c1[256], c2[256];
        poly_rand(a, n);
        poly_rand(b, n);
        for(i=0; i<n; i++) {
            a[i] = (a[i] % q + q) % q;
            b[i] = (b[i] % q + q) % q;
        }
        
        POLY_COPY(a_tmp, a, n);
        POLY_COPY(b_tmp, b, n);
        poly_mult_ntt(c1, a_tmp, b_tmp, &ctx);
        
        POLY_COPY(a_tmp, b, n);
        POLY_COPY(b_tmp, a, n);
        poly_mult_ntt(c2, a_tmp, b_tmp, &ctx);
        
        TEST_ITEM(check_eq(c1, c2, n, q), end_comm);
    } TEST_FINAL
    end_comm:
    
    TEST_START("Small polynomial multiplication", 1) {
        POLY_ZERO(a, n);
        a[0] = 1; a[1] = 2;
        POLY_ZERO(b, n);
        b[0] = 3; b[1] = 4;
        
        POLY_COPY(a_tmp, a, n);
        POLY_COPY(b_tmp, b, n);
        poly_mult_ntt(c, a_tmp, b_tmp, &ctx);
        
        int correct = ((c[0] % q + q) % q == 3) && 
                      ((c[1] % q + q) % q == 10) && 
                      ((c[2] % q + q) % q == 8);
        for(i=3; i<n; i++) {
            if((c[i] % q + q) % q != 0) {
                correct = 0;
                break;
            }
        }
        TEST_ITEM(correct, end_small);
    } TEST_FINAL
    end_small:
    
    TEST_START("Multiplication mod x^256+1", 1) {
        POLY_ZERO(a, n);
        a[255] = 1;
        POLY_ZERO(b, n);
        b[1] = 1;
        
        POLY_COPY(a_tmp, a, n);
        POLY_COPY(b_tmp, b, n);
        poly_mult_ntt(c, a_tmp, b_tmp, &ctx);
        
        int correct = ((c[0] % q + q) % q == ((q - 1) % q));
        for(i=1; i<n && correct; i++) {
            if((c[i] % q + q) % q != 0) correct = 0;
        }
        TEST_ITEM(correct, end_mulx);
    } TEST_FINAL
    end_mulx:
    
    printf("\n[Benchmark] Polynomial Multiplication Performance:\n");
    
    poly_rand(a, n);
    poly_rand(b, n);
    for(i=0; i<n; i++) {
        a[i] = (a[i] % q + q) % q;
        b[i] = (b[i] % q + q) % q;
    }
    
    BENCH_START("Poly Mul (Complete)", 10) {
        POLY_COPY(a_tmp, a, n);
        POLY_COPY(b_tmp, b, n);
        BENCH_ITEM(poly_mult_ntt(c, a_tmp, b_tmp, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    lssi_t a_ntt[256], b_ntt[256];
    POLY_COPY(a_ntt, a, n);
    POLY_COPY(b_ntt, b, n);
    
    BENCH_START("  - NTT (single)", 10) {
        BENCH_ITEM(poly_ntt(a_ntt, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    BENCH_START("  - NTT x2 (for a&b)", 10) {
        POLY_COPY(a_ntt, a, n);
        POLY_COPY(b_ntt, b, n);
        BENCH_ITEM({poly_ntt(a_ntt, &ctx); poly_ntt(b_ntt, &ctx);}, 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_ntt(a_ntt, &ctx);
    poly_ntt(b_ntt, &ctx);
    
    BENCH_START("  - Pointwise mul", 10) {
        BENCH_ITEM(poly_pointwise(c, a_ntt, b_ntt, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    BENCH_START("  - INTT", 10) {
        BENCH_ITEM(poly_invntt(c, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

int main() {
    test_init();
    //test_no_table();
    test_using_tables();
    test_poly_mul_ntt();
    return 0;
}