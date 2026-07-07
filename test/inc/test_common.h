#ifndef TEST_COMMON_H
#define TEST_COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#if defined(_WIN32)
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif
#include "../../inc/types.h"
#include "../../inc/poly.h"
#include "../../inc/rand.h"

#define TEST_COEFS     256
#define TEST_Q_KYBER   3329
#define TEST_Q_DILI    8380417
#define TEST_QBITS     12

#define TEST_START(_DESCRIPTION, _TESTS)                                                    \
    printf("Testing if " _DESCRIPTION "...%*c", (int)(64 - strlen(_DESCRIPTION)), ' ');     \
    for (int _i = 0; _i < _TESTS; _i++)                                                     \
    {

#define TEST_ITEM(_CONDITION, _LABEL) \
    if (!(_CONDITION))                \
    {                                 \
        test_fail();                  \
        goto _LABEL;                  \
    }

#define TEST_FINAL \
    }              \
    test_pass();

#define BENCH_START(_LABEL, _BENCHS)                                       \
    bench_reset();                                                         \
    printf("BENCH: " _LABEL "%*c = ", (int)(32 - strlen(_LABEL)), ' ');   \
    for (int _b = 0; _b < _BENCHS; _b++)                                   \
    {

#define BENCH_FINAL(_BENCHS, _ROUNDS)  \
    }                                  \
    bench_compute(_BENCHS);            \
    bench_print(_ROUNDS);

#define BENCH_ITEM(_FUNCTION, _ROUNDS)  \
    bench_before();                     \
    for (int _r = 0; _r < _ROUNDS; _r++) \
    {                                   \
        _FUNCTION;                      \
    }                                   \
    bench_after();
uint64_t _bench_before;
static uint64_t _bench_after;
static uint64_t _bench_total;

static inline void test_fail(void) {
    printf("[FAIL]\n");
}

static inline void test_pass(void) {
    printf("[PASS]\n");
}

static inline void bench_reset(void) {
    _bench_before = 0;
    _bench_after = 0;
    _bench_total = 0;
}

static inline void bench_before(void) {
    _bench_before = __rdtsc();
}

static inline void bench_after(void) {
    _bench_after = __rdtsc();
    _bench_total += (_bench_after - _bench_before);
}

static inline void bench_compute(int benches) {
    if (benches > 0) {
        _bench_total = _bench_total / benches;
    }
}

static inline void bench_print(int rounds) {
    uint64_t cycles_per_op = _bench_total / rounds;
    
    printf("%llu cycles total | %llu cycles/op\n", 
           (unsigned long long)_bench_total, 
           (unsigned long long)cycles_per_op);
}

static inline unsigned long long bench_total(void) {
    return (unsigned long long)_bench_total;
}

static inline void print_poly(const char* label, lssi_t *a, int len) {
    if (len <= 3) {
        printf("%s: [", label);
        for (int i = 0; i < len; i++) {
            printf("%d%s", a[i], (i < len-1) ? ", " : "]\n");
        }
    } else {
        printf("%s: [%d, %d, %d, ..., %d] (%d elements)\n", 
               label, a[0], a[1], a[2], a[len-1], len);
    }
}

static inline int check_poly_range(lssi_t *a, int len, lssi_t min, lssi_t max) {
    for (int i = 0; i < len; i++) {
        if (a[i] < min || a[i] > max) {
            return 0;
        }
    }
    return 1;
}

static inline int poly_equals(const lssi_t *a, const lssi_t *b, int len) {
    return memcmp(a, b, len * sizeof(lssi_t)) == 0;
}

static inline void test_init(void) {
    uint8_t seed[] = {0x12, 0x34, 0x56, 0x78};
    rand_seed(seed, sizeof(seed));
    printf("Test environment initialized\n\n");
}

#endif // TEST_COMMON_H
