#include "../inc/test_common.h"
#include "../../inc/poly.h"

#define NTESTS 5000
typedef lssi_t poly[MAX_POLY_COEFS];

int poly_copy(lssi_t* dst, const lssi_t* src, int coefs){
    for(int i=0;i<coefs;++i){
        dst[i] = src[i];
    }
    return 0;
}

int poly_eq(const lssi_t* a, const lssi_t* b, int coefs) {
    for (int i = 0; i < coefs; ++i) {
        if (a[i] != b[i])
            return 0;
    }
    return 1;
}

static void poly_naivemul(poly *c, const poly *a, const poly *b, lp_ctx Q) {
    int i, j;
    int n = Q.n;
    int64_t q = Q.q;
    int64_t r[2 * MAX_POLY_COEFS] = {0};
    
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            r[i+j] = (r[i+j] + (int64_t)(*a)[i] * (*b)[j]) % q;
    
    for(i = n; i < 2*n; i++)
        r[i-n] = (r[i-n] - r[i]) % q;
    
    for(i = 0; i < n; i++) {
        (*c)[i] = r[i];
        if((*c)[i] < 0) (*c)[i] += q;
    }
}

void poly_ntt_mul8(poly *c, const poly *a, const poly *b, lp_ctx Q) {
    lp_ctx ctx = Q;
    lssi_t * a_tmp = (lssi_t*)malloc(sizeof(lssi_t) * ctx.n);
    lssi_t * b_tmp = (lssi_t*)malloc(sizeof(lssi_t) * ctx.n);
    poly_copy(a_tmp, *a, ctx.n);
    poly_copy(b_tmp, *b, ctx.n);
    poly_ntt(a_tmp, &ctx);
    poly_ntt(b_tmp, &ctx);
    poly_pointwise(*c, a_tmp, b_tmp, &ctx);
    poly_invntt(*c, &ctx);
    free(a_tmp);
    free(b_tmp);
}

void poly_ntt_mul7(poly *c, const poly *a, const poly *b, lp_ctx Q) {
    lp_ctx ctx = Q;
    lssi_t * a_tmp = (lssi_t*)malloc(sizeof(lssi_t) * ctx.n);
    lssi_t * b_tmp = (lssi_t*)malloc(sizeof(lssi_t) * ctx.n);
    poly_copy(a_tmp, *a, ctx.n);
    poly_copy(b_tmp, *b, ctx.n);
    poly_ntt(a_tmp, &ctx);
    poly_ntt(b_tmp, &ctx);
    poly_base_pointwise(*c, a_tmp, b_tmp, &ctx);
    ctx.layer = 7;
    poly_invntt(*c, &ctx);
    free(a_tmp);
    free(b_tmp);
}

void correctness_check_7_dili() {
    printf("\n============================================\n");
    printf("  Test Suite: point_mul_7_dili\n");
    printf("============================================\n");
    int i = 0, pass = 0;
    lp_ctx ctx;
    INIT_NTT_DILI(&ctx);
    ctx.layer = 7;

    lssi_t a[MAX_POLY_COEFS], b[MAX_POLY_COEFS];
    lssi_t c_naive[MAX_POLY_COEFS], c_ntt[MAX_POLY_COEFS];
    memset(a, 0, sizeof(a));
    memset(b, 0, sizeof(b));
    for (int i = 0; i < NTESTS; i++)
    {
        for(int j=0;j<256;++j){
            a[j] = rand() % ctx.q;
            b[j] = rand() % ctx.q;
        }

        poly_naivemul((poly*)&c_naive, (poly*)&a, (poly*)&b, ctx);
        poly_ntt_mul7((poly*)&c_ntt, (poly*)&a, (poly*)&b, ctx);

        for(int j=0;j<ctx.n;++j){
            c_naive[j] = (c_naive[j]%ctx.q + ctx.q)%ctx.q;
            c_ntt[j] = (c_ntt[j]%ctx.q + ctx.q)%ctx.q;
        }

        if(poly_eq((const lssi_t*)&c_naive, (const lssi_t*)&c_ntt, ctx.n)){
            pass++;
            //printf("Test %d passed!\n", i+1);
        } else {
            printf("Test %d failed!\n", i+1);
            for(int k=0; k<ctx.n; ++k) {
                lssi_t va = (c_naive[k] % ctx.q + ctx.q) % ctx.q;
                lssi_t vb = (c_ntt[k] % ctx.q + ctx.q) % ctx.q;
                if(va != vb) {
                printf("First mismatch at index %d: Naive=%lld, NTT=%lld\n", k, (long long)va, (long long)vb);
                break;
                }
            }
            print_poly("Naive Result", c_naive, ctx.n);
            print_poly("NTT Result", c_ntt, ctx.n);
        }
    }
    printf("Passed %d/%d tests.", pass, NTESTS);
    if(pass == NTESTS) {
        test_pass();
    }

    int n = ctx.n;
    printf("\n[Benchmark] Polynomial Multiplication Components:\n");
    lssi_t c[256];
    lssi_t a_ntt[256], b_ntt[256];
  // 7-layer components
    ctx.layer = 7;
    poly_copy(a_ntt, a, n);
    poly_copy(b_ntt, b, n);
    
    BENCH_START("7-layer: NTT*2", 10) {
        BENCH_ITEM({
            poly_copy(a_ntt, a, n);
            poly_copy(b_ntt, b, n);
            poly_ntt(a_ntt, &ctx);
            poly_ntt(b_ntt, &ctx);
        }, 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_ntt(a_ntt, &ctx);
    poly_ntt(b_ntt, &ctx);
    BENCH_START("7-layer: Base Pointwise", 10) {
        BENCH_ITEM(poly_base_pointwise(c, a_ntt, b_ntt, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_base_pointwise(c, a_ntt, b_ntt, &ctx);
    ctx.layer = 7;
    BENCH_START("7-layer: INTT", 10) {
        BENCH_ITEM(poly_invntt(c, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    ctx.layer = 7;
    BENCH_START("7-layer: Total Mul", 10) {
        BENCH_ITEM(poly_ntt_mul7((poly*)&c, (poly*)&a, (poly*)&b, ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

void correctness_check_7_kyber() {
    printf("\n============================================\n");
    printf("  Test Suite: point_mul_7_kyber\n");
    printf("============================================\n");
    int i = 0, pass = 0;
    lp_ctx ctx;
    INIT_NTT_KYBER(&ctx);
    ctx.layer = 7;

    lssi_t a[MAX_POLY_COEFS], b[MAX_POLY_COEFS];
    lssi_t c_naive[MAX_POLY_COEFS], c_ntt[MAX_POLY_COEFS];
    memset(a, 0, sizeof(a));
    memset(b, 0, sizeof(b));
    for (int i = 0; i < NTESTS; i++)
    {
        for(int j=0;j<256;++j){
            a[j] = rand() % ctx.q;
            b[j] = rand() % ctx.q;
        }

        poly_naivemul((poly*)&c_naive, (poly*)&a, (poly*)&b, ctx);
        poly_ntt_mul7((poly*)&c_ntt, (poly*)&a, (poly*)&b, ctx);

        for(int j=0;j<ctx.n;++j){
            c_naive[j] = (c_naive[j]%ctx.q + ctx.q)%ctx.q;
            c_ntt[j] = (c_ntt[j]%ctx.q + ctx.q)%ctx.q;
        }

        if(poly_eq((const lssi_t*)&c_naive, (const lssi_t*)&c_ntt, ctx.n)){
            pass++;
            //printf("Test %d passed!\n", i+1);
        } else {
            printf("Test %d failed!\n", i+1);
            for(int k=0; k<ctx.n; ++k) {
                lssi_t va = (c_naive[k] % ctx.q + ctx.q) % ctx.q;
                lssi_t vb = (c_ntt[k] % ctx.q + ctx.q) % ctx.q;
                if(va != vb) {
                printf("First mismatch at index %d: Naive=%lld, NTT=%lld\n", k, (long long)va, (long long)vb);
                break;
                }
            }
            print_poly("Naive Result", c_naive, ctx.n);
            print_poly("NTT Result", c_ntt, ctx.n);
        }
    }
    printf("Passed %d/%d tests.", pass, NTESTS);
    if(pass == NTESTS) {
        test_pass();
    }

    int n = ctx.n;
    printf("\n[Benchmark] Polynomial Multiplication Components:\n");
    lssi_t c[256];
    lssi_t a_ntt[256], b_ntt[256];
  // 7-layer components
    ctx.layer = 7;
    poly_copy(a_ntt, a, n);
    poly_copy(b_ntt, b, n);
    
    BENCH_START("7-layer: NTT*2", 10) {
        BENCH_ITEM({
            poly_copy(a_ntt, a, n);
            poly_copy(b_ntt, b, n);
            poly_ntt(a_ntt, &ctx);
            poly_ntt(b_ntt, &ctx);
        }, 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_ntt(a_ntt, &ctx);
    poly_ntt(b_ntt, &ctx);
    BENCH_START("7-layer: Base Pointwise", 10) {
        BENCH_ITEM(poly_base_pointwise(c, a_ntt, b_ntt, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_base_pointwise(c, a_ntt, b_ntt, &ctx);
    ctx.layer = 7;
    BENCH_START("7-layer: INTT", 10) {
        BENCH_ITEM(poly_invntt(c, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    ctx.layer = 7;
    BENCH_START("7-layer: Total Mul", 10) {
        BENCH_ITEM(poly_ntt_mul7((poly*)&c, (poly*)&a, (poly*)&b, ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

void correctness_check_8_dili() {
    printf("\n============================================\n");
    printf("  Test Suite: point_mul_8_dili\n");
    printf("============================================\n");
    int i = 0, pass = 0;
    lp_ctx ctx;
    INIT_NTT_DILI(&ctx);
    lssi_t a[MAX_POLY_COEFS], b[MAX_POLY_COEFS];
    lssi_t c_naive[MAX_POLY_COEFS], c_ntt[MAX_POLY_COEFS];
    memset(a, 0, sizeof(a));
    memset(b, 0, sizeof(b));
    for (int i = 0; i < NTESTS; i++)
    {
        for(int j=0;j<256;++j){
        a[j] = rand() % ctx.q;
        b[j] = rand() % ctx.q;
    }
    
    poly_naivemul((poly*)&c_naive, (poly*)&a, (poly*)&b, ctx);
    poly_ntt_mul8((poly*)&c_ntt, (poly*)&a, (poly*)&b, ctx);

    for(int j=0;j<ctx.n;++j){
      c_naive[j] = (c_naive[j]%ctx.q + ctx.q)%ctx.q;
      c_ntt[j] = (c_ntt[j]%ctx.q + ctx.q)%ctx.q;
    }

    if(poly_eq((const lssi_t*)&c_naive, (const lssi_t*)&c_ntt, ctx.n)){
        pass++;
        //printf("Test %d passed!\n", i+1);
    } else {
        printf("Test %d failed!\n", i+1);
        for(int k=0; k<ctx.n; ++k) {
            lssi_t va = (c_naive[k] % ctx.q + ctx.q) % ctx.q;
            lssi_t vb = (c_ntt[k] % ctx.q + ctx.q) % ctx.q;
            if(va != vb) {
                printf("First mismatch at index %d: Naive=%lld, NTT=%lld\n", k, (long long)va, (long long)vb);
                break;
            }
        }
        print_poly("Naive Result", c_naive, ctx.n);
        print_poly("NTT Result", c_ntt, ctx.n);
    }
  }
  printf("Passed %d/%d tests.", pass, NTESTS);
  if(pass == NTESTS) {
      test_pass();
  }

  int n = ctx.n;
  printf("\n[Benchmark] Polynomial Multiplication Components:\n");
    lssi_t c[256];
    lssi_t a_ntt[256], b_ntt[256];
    poly_rand(a, n);
    poly_rand(b, n);
    ctx.layer = 8;
    poly_copy(a_ntt, a, n);
    poly_copy(b_ntt, b, n);
    
    BENCH_START("8-layer: NTT*2", 10) {
        BENCH_ITEM({
            poly_copy(a_ntt, a, n);
            poly_copy(b_ntt, b, n);
            poly_ntt(a_ntt, &ctx);
            poly_ntt(b_ntt, &ctx);
        }, 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_ntt(a_ntt, &ctx);
    poly_ntt(b_ntt, &ctx);
    BENCH_START("8-layer: Pointwise", 10) {
        BENCH_ITEM(poly_pointwise(c, a_ntt, b_ntt, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    poly_pointwise(c, a_ntt, b_ntt, &ctx);
    BENCH_START("8-layer: INTT", 10) {
        BENCH_ITEM(poly_invntt(c, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    BENCH_START("8-layer: Total Mul", 10) {
        BENCH_ITEM(poly_ntt_mul8((poly*)&c, (poly*)&a, (poly*)&b, ctx), 10000);
    } BENCH_FINAL(10, 10000);
}

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

void test_ntt(){
    printf("\n============================================\n");
    printf("  Test Suite: NTT\n");
    printf("============================================\n");
    lp_ctx ctx;
    int i;
    INIT_NTT_DILI(&ctx);
    lssi_t q = ctx.q;
    int n = ctx.n;
    ctx.layer = 7;
    lssi_t a[256], res[256], original[256];

    int passed = 0;
    TEST_START("NTT/INTT reversibility (Dilithium)", 1) {
    for(int j=0;j<100;j++){
        for(i=0; i<n; i++) {
            a[i] = rand() % q;
            original[i] = a[i];
        }
        poly_ntt(a, &ctx);
        poly_invntt(a, &ctx);
        for(i=0; i<n; i++) {
            a[i] = FROM_MONTGOMERY(a[i], q, ctx.qc);
        }
        
        if(!check_eq(a, original, n, q)) {
            printf(" Reversibility failed on iteration %d\n", j);
            print_poly("Original", original, n);
            print_poly("After NTT->INTT", a, n);
        } else {
            passed++;
        }
    } 
    }TEST_FINAL
    printf("Passed %d/100 reversibility tests.\n", passed);

    printf("\n[Benchmark] Dilithium NTT Performance:\n");
    BENCH_START("INIT", 10) {
        BENCH_ITEM(INIT_NTT_DILI(&ctx), 10000);
    } BENCH_FINAL(10, 10000);

    ctx.layer = 8;
    BENCH_START("NTT (8-layer)", 10) {
        BENCH_ITEM(poly_ntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    ctx.layer = 8;
    BENCH_START("INTT (8-layer)", 10) {
        BENCH_ITEM(poly_invntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    ctx.layer = 7;
    BENCH_START("NTT (7-layer)", 10) {
        BENCH_ITEM(poly_ntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);

    ctx.layer = 7;
    BENCH_START("INTT (7-layer)", 10) {
        BENCH_ITEM(poly_invntt(a, &ctx), 10000);
    } BENCH_FINAL(10, 10000);
    
    
}

int main() {
    test_init();
    test_ntt();
    correctness_check_8_dili();
    correctness_check_7_dili();
    correctness_check_7_kyber();
    return 0;
}