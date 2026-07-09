#include "test_common.h"
#include "kyber.h"
#include "kyber_pack.h"
#include "kyber_sample.h"
#include <string.h>
#include <stdlib.h>

#define KYBER_PERF_ITERATIONS 10000

static void fill_bytes(uint8_t *buf, size_t len, uint8_t seed) {
    for(size_t i = 0; i < len; i++) {
        buf[i] = (uint8_t)(seed + 17u * (uint8_t)i + (uint8_t)(i >> 3));
    }
}

static int bytes_eq(const uint8_t *a, const uint8_t *b, size_t len) {
    return memcmp(a, b, len) == 0;
}

static void print_case(const char *name, int ok) {
    printf("%-72s%s\n", name, ok ? "[PASS]" : "[FAIL]");
}

static inline double kyber_now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static int check_param_one(int mode, int k, int eta1, int eta2,
                           int du, int dv, int ek, int dk_pke, int ct) {
    kyber_param param;

    kyber_init_param(&param, mode);
    return param.mode == mode && param.k == k && param.eta1 == eta1 &&
           param.eta2 == eta2 && param.du == du && param.dv == dv &&
           param.ek_bytes == ek && param.dk_pke_bytes == dk_pke &&
           param.ct_bytes == ct;
}

static int test_params(void) {
    return check_param_one(KYBER_MODE_2, 2, 3, 2, 10, 4, 800, 768, 768) &&
           check_param_one(KYBER_MODE_3, 3, 2, 2, 10, 4, 1184, 1152, 1088) &&
           check_param_one(KYBER_MODE_5, 4, 2, 2, 11, 5, 1568, 1536, 1568);
}

static int test_byte_encode_decode_d(int d) {
    poly a;
    poly b;
    uint8_t buf[KYBER_POLYBYTES];
    int mask = (1 << d) - 1;

    for(int i = 0; i < KYBER_N; i++) {
        a[i] = (d == KYBER_QBITS) ? (i * 29 + 7) % KYBER_Q
                                  : ((i * 29 + 7) & mask);
    }

    kyber_byte_encode(buf, a, d);
    kyber_byte_decode(b, buf, d);

    for(int i = 0; i < KYBER_N; i++) {
        if(b[i] != a[i]) {
            printf("  ByteEncode%d mismatch at %d: got=%d expect=%d\n",
                   d, i, b[i], a[i]);
            return 0;
        }
    }
    return 1;
}

static int test_byte_encode_decode(void) {
    int ds[] = {1, 4, 5, 10, 11, 12};

    for(size_t i = 0; i < sizeof(ds) / sizeof(ds[0]); i++) {
        if(!test_byte_encode_decode_d(ds[i])) {
            return 0;
        }
    }
    return 1;
}

static int test_compress_decompress(void) {
    int ds[] = {1, 4, 5, 10, 11};

    for(size_t i = 0; i < sizeof(ds) / sizeof(ds[0]); i++) {
        int d = ds[i];
        int bound = 1 << d;
        for(int x = 0; x < bound; x++) {
            lssi_t y = kyber_decompress_one((lssi_t)x, d);
            lssi_t z = kyber_compress_one(y, d);
            if(z != x) {
                printf("  Compress(Decompress(.)) failed d=%d x=%d z=%d\n",
                       d, x, z);
                return 0;
            }
        }
    }
    return 1;
}

static int test_sampling(void) {
    lp_ctx ctx;
    uint8_t seed[34];
    uint8_t sigma[KYBER_SYMBYTES];
    poly a;
    poly e;

    INIT_NTT_KYBER(&ctx);
    fill_bytes(seed, sizeof(seed), 0x31);
    fill_bytes(sigma, sizeof(sigma), 0xA5);

    kyber_sample_ntt(a, seed, &ctx);
    for(int i = 0; i < KYBER_N; i++) {
        if(a[i] < 0 || a[i] >= KYBER_Q) {
            printf("  SampleNTT out of range at %d: %d\n", i, a[i]);
            return 0;
        }
    }

    kyber_sample_noise(e, sigma, 7, 3);
    for(int i = 0; i < KYBER_N; i++) {
        if(e[i] < 0 || e[i] >= KYBER_Q) {
            printf("  CBD out of range at %d: %d\n", i, e[i]);
            return 0;
        }
    }
    return 1;
}

static int test_ek_modulus_check(void) {
    kyber_param param;
    uint8_t ek[KYBER_MAX_EK_BYTES];

    kyber_init_param(&param, KYBER_MODE_2);
    memset(ek, 0, sizeof(ek));
    ek[0] = 0xFF;
    ek[1] = 0x0F;
    return !kyber_check_ek_modulus(ek, param);
}

static int test_pke_mode(int mode) {
    kyber_param param;
    uint8_t d[KYBER_SYMBYTES];
    uint8_t r[KYBER_SYMBYTES];
    uint8_t m[KYBER_SYMBYTES];
    uint8_t m2[KYBER_SYMBYTES];
    uint8_t ek[KYBER_MAX_EK_BYTES];
    uint8_t dk_pke[KYBER_MAX_DK_PKE_BYTES];
    uint8_t ct[KYBER_MAX_CT_BYTES];

    kyber_init_param(&param, mode);
    fill_bytes(d, sizeof(d), (uint8_t)(mode & 0xFF));
    fill_bytes(r, sizeof(r), (uint8_t)((mode >> 1) & 0xFF));
    fill_bytes(m, sizeof(m), (uint8_t)((mode >> 2) & 0xFF));

    if(kyber_pke_keygen(ek, dk_pke, d, param) != 0) {
        printf("  %s keygen failed\n", kyber_mode_name(mode));
        return 0;
    }
    if(kyber_pke_encrypt(ct, ek, m, r, param) != 0) {
        printf("  %s encrypt failed\n", kyber_mode_name(mode));
        return 0;
    }
    if(kyber_pke_decrypt(m2, dk_pke, ct, param) != 0) {
        printf("  %s decrypt failed\n", kyber_mode_name(mode));
        return 0;
    }
    if(!bytes_eq(m, m2, sizeof(m))) {
        printf("  %s PKE message mismatch\n", kyber_mode_name(mode));
        return 0;
    }
    return 1;
}

static int test_pke_all_modes(void) {
    return test_pke_mode(KYBER_MODE_2) &&
           test_pke_mode(KYBER_MODE_3) &&
           test_pke_mode(KYBER_MODE_5);
}

static int prepare_perf_inputs(uint8_t (*d_list)[KYBER_SYMBYTES],
                               uint8_t (*r_list)[KYBER_SYMBYTES],
                               uint8_t (*m_list)[KYBER_SYMBYTES],
                               int mode) {
    for(int i = 0; i < KYBER_PERF_ITERATIONS; i++) {
        fill_bytes(d_list[i], KYBER_SYMBYTES, (uint8_t)(0x41 + mode + i));
        fill_bytes(r_list[i], KYBER_SYMBYTES, (uint8_t)(0x53 + mode + i));
        fill_bytes(m_list[i], KYBER_SYMBYTES, (uint8_t)(0x67 + mode + i));
    }
    return 1;
}

static void test_perf_mode(int mode) {
    kyber_param param;
    uint8_t (*d_list)[KYBER_SYMBYTES];
    uint8_t (*r_list)[KYBER_SYMBYTES];
    uint8_t (*m_list)[KYBER_SYMBYTES];
    uint8_t (*ct_list)[KYBER_MAX_CT_BYTES];
    uint8_t ek[KYBER_MAX_EK_BYTES];
    uint8_t dk_pke[KYBER_MAX_DK_PKE_BYTES];
    uint8_t ct[KYBER_MAX_CT_BYTES];
    uint64_t cyc_start, cyc_end;
    double ms_start, ms_end;

    kyber_init_param(&param, mode);

    printf("\n============================================\n");
    printf("  Test: %s Performance (%d iterations)\n",
           kyber_mode_name(mode), KYBER_PERF_ITERATIONS);
    printf("============================================\n");

    d_list = malloc((size_t)KYBER_PERF_ITERATIONS * sizeof(*d_list));
    r_list = malloc((size_t)KYBER_PERF_ITERATIONS * sizeof(*r_list));
    m_list = malloc((size_t)KYBER_PERF_ITERATIONS * sizeof(*m_list));
    ct_list = malloc((size_t)KYBER_PERF_ITERATIONS * sizeof(*ct_list));
    if(!d_list || !r_list || !m_list || !ct_list) {
        printf("  [FAIL] Allocation failed for performance inputs\n");
        free(d_list);
        free(r_list);
        free(m_list);
        free(ct_list);
        return;
    }

    prepare_perf_inputs(d_list, r_list, m_list, mode);

    cyc_start = __rdtsc();
    ms_start = kyber_now_ms();
    for(int i = 0; i < KYBER_PERF_ITERATIONS; i++) {
        kyber_pke_keygen(ek, dk_pke, d_list[i], param);
    }
    cyc_end = __rdtsc();
    ms_end = kyber_now_ms();
    printf("  KeyGen : %12llu cycles/op | %10.4f ms/op  (total %10.2f ms)\n",
           (unsigned long long)((cyc_end - cyc_start) / KYBER_PERF_ITERATIONS),
           (ms_end - ms_start) / KYBER_PERF_ITERATIONS,
           ms_end - ms_start);

    cyc_start = __rdtsc();
    ms_start = kyber_now_ms();
    for(int i = 0; i < KYBER_PERF_ITERATIONS; i++) {
        kyber_pke_encrypt(ct, ek, m_list[i], r_list[i], param);
    }
    cyc_end = __rdtsc();
    ms_end = kyber_now_ms();
    printf("  Encrypt: %12llu cycles/op | %10.4f ms/op  (total %10.2f ms)\n",
           (unsigned long long)((cyc_end - cyc_start) / KYBER_PERF_ITERATIONS),
           (ms_end - ms_start) / KYBER_PERF_ITERATIONS,
           ms_end - ms_start);

    for(int i = 0; i < KYBER_PERF_ITERATIONS; i++) {
        kyber_pke_encrypt(ct_list[i], ek, m_list[i], r_list[i], param);
    }

    cyc_start = __rdtsc();
    ms_start = kyber_now_ms();
    for(int i = 0; i < KYBER_PERF_ITERATIONS; i++) {
        kyber_pke_decrypt(m_list[i], dk_pke, ct_list[i], param);
    }
    cyc_end = __rdtsc();
    ms_end = kyber_now_ms();
    printf("  Decrypt: %12llu cycles/op | %10.4f ms/op  (total %10.2f ms)\n",
           (unsigned long long)((cyc_end - cyc_start) / KYBER_PERF_ITERATIONS),
           (ms_end - ms_start) / KYBER_PERF_ITERATIONS,
           ms_end - ms_start);

    free(d_list);
    free(r_list);
    free(m_list);
    free(ct_list);
}

static void test_perf_all_modes(void) {
    printf("\n============================================\n");
    printf("  Test: All Kyber Modes Performance\n");
    printf("============================================\n");

    test_perf_mode(KYBER_MODE_2);
    test_perf_mode(KYBER_MODE_3);
    test_perf_mode(KYBER_MODE_5);
}

int main(void) {
    int ok = 1;

    test_init();
    printf("============================================\n");
    printf("  Test Suite: Kyber\n");
    printf("============================================\n");

    int r;

    r = test_params();
    print_case("Kyber parameter sets", r);
    ok &= r;

    r = test_byte_encode_decode();
    print_case("Kyber ByteEncode/ByteDecode", r);
    ok &= r;

    r = test_compress_decompress();
    print_case("Kyber Compress/Decompress consistency", r);
    ok &= r;

    r = test_sampling();
    print_case("Kyber SampleNTT/CBD ranges", r);
    ok &= r;

    r = test_ek_modulus_check();
    print_case("Kyber public-key modulus check", r);
    ok &= r;

    r = test_pke_all_modes();
    print_case("Kyber K-PKE roundtrip in all modes", r);
    ok &= r;

    test_perf_all_modes();

    return ok ? 0 : 1;
}
