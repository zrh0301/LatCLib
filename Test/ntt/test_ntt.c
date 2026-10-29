#include "test.h"
#include "polyvec.h"

/* Exhaustive mathematical contracts live in test_contexts, linked only to
 * algebra. Here check the compatibility presets used by existing callers. */
static int check_preset(int kyber)
{
    lp_ctx plan;
    poly a, b, original;
    if (kyber) INIT_NTT_KYBER(&plan);
    else INIT_NTT_DILI(&plan);
    int32_t q = plan.ring.field.q;
    if (plan.ring.n != 256 || q != (kyber ? 3329 : 8380417) ||
        plan.layers != (kyber ? 7 : 8) ||
        plan.block_degree != (kyber ? 2 : 1)) return 0;
    if (plan.zetas[1] != (kyber ? -950 : 25847)) return 0;
    for (int i = 0; i < 256; ++i) original[i] = (i * 197 + 23) % q;
    memcpy(a, original, sizeof(a));
    poly_ntt(a, &plan);
    memcpy(b, a, sizeof(b));
    poly_invntt(b, &plan);
    for (int i = 0; i < 256; ++i) {
        if (fp32_from_mont(b[i], &plan.ring.field) != original[i]) return 0;
    }
    return 1;
}

static void benchmark_preset(int kyber)
{
    lp_ctx plan;
    poly a, b, c, transformed_a, transformed_b;
    if (kyber) INIT_NTT_KYBER(&plan);
    else INIT_NTT_DILI(&plan);
    for (int i = 0; i < 256; ++i) {
        a[i] = i;
        b[i] = 255 - i;
    }
    memcpy(transformed_a, a, sizeof(a));
    memcpy(transformed_b, b, sizeof(b));
    poly_ntt(transformed_a, &plan);
    poly_ntt(transformed_b, &plan);

    printf("\n%s polynomial operators\n", kyber ? "Kyber" : "Dilithium");
    BENCH_START("poly_ntt", 10) {
        BENCH_ITEM({
            poly input;
            memcpy(input, a, sizeof(input));
            input[_r & 255]++;
            poly_ntt(input, &plan);
            c[_r & 255] = input[_r & 255];
        }, 1000);
    } BENCH_FINAL(10, 1000);
    BENCH_START("poly_mul_ntt", 10) {
        BENCH_ITEM({
            poly_mul_ntt(c, transformed_a, transformed_b, &plan);
            c[_r & 255] ^= _r;
        }, 1000);
    } BENCH_FINAL(10, 1000);
    BENCH_START("poly_invntt", 10) {
        BENCH_ITEM({
            memcpy(c, transformed_a, sizeof(c));
            c[_r & 255] ^= _r;
            poly_invntt(c, &plan);
        }, 1000);
    } BENCH_FINAL(10, 1000);
}

int main(int argc, char **argv)
{
    int ok = check_preset(0) && check_preset(1);
    printf("NTT compatibility presets: %s\n", ok ? "[PASS]" : "[FAIL]");
    if (argc == 2 && strcmp(argv[1], "--bench") == 0) {
        benchmark_preset(0);
        benchmark_preset(1);
    }
    return ok ? 0 : 1;
}
