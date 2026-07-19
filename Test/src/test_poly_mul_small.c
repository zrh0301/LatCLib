#include "test_common.h"

/*
 * poly_mul_small 正确性测试
 *
 * 思路：用一份朴素的 O(n^2) 负循环卷积（schoolbook negacyclic convolution）
 * 作为参考实现，和 poly_mul_small 在随机输入上逐系数比较，凡是不一致就判为
 * 失败并打印出第一处不一致的位置，方便定位 bug。
 */

/* 参考实现里用到的、非常数时间的朴素取模居中，仅用于测试，不追求侧信道安全 */
static lssi_t ref_center_mod(lsdi_t x, lssi_t q) {
    lsdi_t r = x % q;
    if (r < 0) r += q;
    if (r > (q - 1) / 2) r -= q;
    return (lssi_t)r;
}

/* 朴素负循环卷积：u = c * a mod (x^n + 1, q)，用 int64 累加避免溢出 */
static void ref_poly_mul(lssi_t* u, const lssi_t* c, const lssi_t* a, int n, lssi_t q) {
    lsdi_t acc[TEST_COEFS];
    int i, j, k;

    for (k = 0; k < n; k++) acc[k] = 0;

    for (i = 0; i < n; i++) {
        if (c[i] == 0) continue;
        for (j = 0; j < n; j++) {
            lsdi_t term = (lsdi_t)c[i] * (lsdi_t)a[j];
            k = i + j;
            if (k >= n) {           /* x^n == -1，回绕时变号 */
                k -= n;
                term = -term;
            }
            acc[k] += term;
        }
    }

    for (k = 0; k < n; k++) u[k] = ref_center_mod(acc[k], q);
}

/* 生成恰有 tau 个非零系数（值 +-1，随机符号、随机位置）的三值多项式 */
static void gen_ternary_poly(lssi_t* c, int n, int tau) {
    int i, placed;

    for (i = 0; i < n; i++) c[i] = 0;

    placed = 0;
    while (placed < tau) {
        uint8_t rb[2];
        int pos;
        rand_bytes(rb, 2);
        pos = rb[0] % n;
        if (c[pos] != 0) continue;
        c[pos] = (rb[1] & 1) ? 1 : -1;
        placed++;
    }
}

/* 生成系数落在 [-U, U] 内的随机小系数多项式 */
static void gen_small_poly(lssi_t* a, int n, lssi_t U) {
    int i;
    for (i = 0; i < n; i++) {
        if (U == 0) {
            a[i] = 0;
            continue;
        }
        uint32_t r;
        rand_bytes((uint8_t*)&r, sizeof(r));
        a[i] = (lssi_t)(r % (uint32_t)(2 * U + 1)) - U;
    }
}

/* 跑一组随机 (c, a) 用例，返回 1 表示 poly_mul_small 与参考实现结果一致 */
static int run_case(lp_ctx* Q, int tau, lssi_t U, const char* label) {
    lssi_t c[TEST_COEFS], a[TEST_COEFS], u_ref[TEST_COEFS], u_fast[TEST_COEFS];
    int n = Q->n;
    int i, mismatches;

    gen_ternary_poly(c, n, tau);
    gen_small_poly(a, n, U);

    ref_poly_mul(u_ref, c, a, n, Q->q);
    poly_mul_small(u_fast, c, a, U, Q);

    if (poly_equals(u_ref, u_fast, n)) {
        return 1;
    }

    mismatches = 0;
    printf("\n  [%s] mismatch (n=%d, q=%d, tau=%d, U=%d):\n", label, n, Q->q, tau, U);
    for (i = 0; i < n && mismatches < 8; i++) {
        if (u_ref[i] != u_fast[i]) {
            printf("    coef[%d]: ref=%d  poly_mul_small=%d\n", i, u_ref[i], u_fast[i]);
            mismatches++;
        }
    }
    return 0;
}

void test_poly_mul_small_correctness(void) {
    lp_ctx Qd, Qk;
    int all_ok = 1;

    ntt_init_dili(&Qd);
    ntt_init_kyber(&Qk);

    printf("=== poly_mul_small correctness (vs. schoolbook reference) ===\n");

    TEST_START("poly_mul_small matches reference (Dilithium q, tau=39, U=2 i.e. eta=2)", 20) {
        int ok = run_case(&Qd, 39, 2, "dili-eta2");
        all_ok &= ok;
        TEST_ITEM(ok, end1);
    } TEST_FINAL
    end1:

    TEST_START("poly_mul_small matches reference (Dilithium q, tau=39, U=4096 i.e. t0-like)", 10) {
        int ok = run_case(&Qd, 39, 4096, "dili-t0");
        all_ok &= ok;
        TEST_ITEM(ok, end2);
    } TEST_FINAL
    end2:

    TEST_START("poly_mul_small matches reference (Kyber q, tau=39, U=3 i.e. eta=3)", 20) {
        int ok = run_case(&Qk, 39, 3, "kyber-eta3");
        all_ok &= ok;
        TEST_ITEM(ok, end3);
    } TEST_FINAL
    end3:

    TEST_START("edge case: c all-zero (tau=0) => u should be all-zero", 5) {
        int ok = run_case(&Qd, 0, 5, "tau0");
        all_ok &= ok;
        TEST_ITEM(ok, end4);
    } TEST_FINAL
    end4:

    TEST_START("edge case: c fully dense (tau=n)", 5) {
        int ok = run_case(&Qd, Qd.n, 4, "tau-full");
        all_ok &= ok;
        TEST_ITEM(ok, end5);
    } TEST_FINAL
    end5:

    TEST_START("edge case: a all-zero (U=0)", 5) {
        int ok = run_case(&Qd, 39, 0, "U0");
        all_ok &= ok;
        TEST_ITEM(ok, end6);
    } TEST_FINAL
    end6:

    printf(all_ok ? "\n>>> poly_mul_small: ALL correctness cases PASSED\n"
                   : "\n>>> poly_mul_small: SOME correctness cases FAILED\n");
}

int main(void) {
    test_init();
    test_poly_mul_small_correctness();
    return 0;
}
