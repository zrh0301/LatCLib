#ifndef TEST_H
#define TEST_H

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
    #include <intrin.h>
#elif defined(__i386__) || defined(__x86_64__)
    #include <x86intrin.h>
#endif
#if defined(__linux__) && (defined(__i386__) || defined(__x86_64__))
    #include <sched.h>
#endif
#include "test_rng.h"
#include "types.h"
#include "poly.h"
#include "poly_sample.h"
#include "rand.h"
#include "dili.h"
#include "kyber.h"
#include "kyber_pack.h"

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

#define BENCH_START(_LABEL, _BENCHS)                                      \
    bench_reset();                                                        \
    printf("BENCH: " _LABEL "%*c = ", (int)(32 - strlen(_LABEL)), ' '); \
    for (int _b = 0; _b < _BENCHS; _b++)                                  \
    {

#define BENCH_FINAL(_BENCHS, _ROUNDS) \
    }                                 \
    bench_compute(_BENCHS);           \
    bench_print(_ROUNDS);

#define BENCH_ITEM(_FUNCTION, _ROUNDS)   \
    bench_before();                      \
    for (int _r = 0; _r < _ROUNDS; _r++) \
    {                                    \
        _FUNCTION;                       \
    }                                    \
    bench_after();

typedef struct {
    uint64_t cycles;
    struct timespec wall_time;
    unsigned int cpu_aux;
#if defined(__linux__) && (defined(__i386__) || defined(__x86_64__))
    cpu_set_t original_affinity;
    int affinity_saved;
#endif
} test_perf_mark;

typedef struct {
    uint64_t cycles;
    double elapsed_ms;
    int cycles_valid;
} test_perf_result;

static uint64_t _bench_cycles;
static double _bench_elapsed_ms;
static int _bench_cycles_valid;
static test_perf_mark _bench_before;

static inline void test_fail(void) {
    printf("[FAIL]\n");
}

static inline void test_pass(void) {
    printf("[PASS]\n");
}

static inline void bench_reset(void) {
    memset(&_bench_before, 0, sizeof(_bench_before));
    _bench_cycles = 0;
    _bench_elapsed_ms = 0.0;
    _bench_cycles_valid = 1;
}

static inline uint64_t test_cpu_cycles(unsigned int *cpu_aux) {
#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64)
    _mm_lfence();
    uint64_t cycles = __rdtscp(cpu_aux);
    _mm_lfence();
    return cycles;
#else
    *cpu_aux = 0;
    return 0;
#endif
}

static inline test_perf_mark test_perf_begin(void) {
    test_perf_mark mark;
    memset(&mark, 0, sizeof(mark));
#if defined(__linux__) && (defined(__i386__) || defined(__x86_64__))
    cpu_set_t allowed;
    if (sched_getaffinity(0, sizeof(allowed), &allowed) == 0) {
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
            if (CPU_ISSET(cpu, &allowed)) {
                cpu_set_t selected;
                CPU_ZERO(&selected);
                CPU_SET(cpu, &selected);
                mark.original_affinity = allowed;
                if (sched_setaffinity(0, sizeof(selected), &selected) == 0)
                    mark.affinity_saved = 1;
                break;
            }
        }
    }
#endif
    mark.cycles = test_cpu_cycles(&mark.cpu_aux);
    timespec_get(&mark.wall_time, TIME_UTC);
    return mark;
}

static inline test_perf_result test_perf_end(test_perf_mark start) {
    test_perf_mark end;
    memset(&end, 0, sizeof(end));
    end.cycles = test_cpu_cycles(&end.cpu_aux);
    timespec_get(&end.wall_time, TIME_UTC);
    test_perf_result result = {0, 0.0, 0};
    int64_t sec = (int64_t)end.wall_time.tv_sec - (int64_t)start.wall_time.tv_sec;
    int64_t nsec = (int64_t)end.wall_time.tv_nsec - (int64_t)start.wall_time.tv_nsec;
    if (start.cycles && end.cycles >= start.cycles && start.cpu_aux == end.cpu_aux) {
        result.cycles = end.cycles - start.cycles;
        result.cycles_valid = 1;
    }
    result.elapsed_ms = (double)sec * 1000.0 + (double)nsec / 1000000.0;
    if (result.elapsed_ms < 0.0) result.elapsed_ms = 0.0;
#if defined(__linux__) && (defined(__i386__) || defined(__x86_64__))
    if (start.affinity_saved)
        (void)sched_setaffinity(0, sizeof(start.original_affinity), &start.original_affinity);
#endif
    return result;
}

static inline void test_perf_print(const char *label, test_perf_result result,
                                   uint64_t iterations) {
    double ms_per_op = iterations ? result.elapsed_ms / (double)iterations : 0.0;
    if (result.cycles_valid) {
        printf("%-24s: %12llu cycles/op | %.6f ms/op (total %9.2f ms)\n",
               label, (unsigned long long)(iterations ? result.cycles / iterations : 0),
               ms_per_op, result.elapsed_ms);
    } else {
        printf("%-24s: %12s cycles/op | %.6f ms/op (total %9.2f ms)\n",
               label, "N/A", ms_per_op, result.elapsed_ms);
    }
}

static inline void bench_before(void) {
    _bench_before = test_perf_begin();
}

static inline void bench_after(void) {
    test_perf_result result = test_perf_end(_bench_before);
    if (result.cycles_valid) _bench_cycles += result.cycles;
    else _bench_cycles_valid = 0;
    _bench_elapsed_ms += result.elapsed_ms;
}

static inline void bench_compute(int benches) {
    if (benches > 0) {
        _bench_cycles /= (uint64_t)benches;
        _bench_elapsed_ms /= (double)benches;
    }
}

static inline void bench_print(int rounds) {
    test_perf_result result = {_bench_cycles, _bench_elapsed_ms, _bench_cycles_valid};
    test_perf_print("", result, (uint64_t)rounds);
}

static inline unsigned long long bench_total(void) {
    return (unsigned long long)_bench_cycles;
}

static inline void print_poly(const char* label, lssi_t *a, int len) {
    if (len <= 3) {
        printf("%s: [", label);
        for (int i = 0; i < len; i++) {
            printf("%d%s", a[i], (i < len-1) ? ", " : "]\n");
        }
    } else {
        printf("%s: [%d, %d, %d, ..., %d] (%d elements)\n", label, a[0], a[1], a[2], a[len-1], len);
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
    test_random_state = UINT32_C(0x6d5a56a9);
    printf("Deterministic test data initialized\n\n");
}

typedef struct {
    const char *name;
    int security_level;
    dili_param param;
    size_t sig_bytes;
    size_t rnd_bytes;
    void (*keygen)(dili_keypair_t *keypair, dili_param param);
    int (*sign)(uint8_t *sig, const uint8_t *msg, int msg_len, const dili_sk *sk, dili_param param);
    int (*sign_with_rnd)(uint8_t *sig, const uint8_t *msg, int msg_len, const dili_sk *sk, dili_param param, const uint8_t rnd[32]);
    int (*sign_with_rng)(uint8_t *sig, const uint8_t *msg, int msg_len, const dili_sk *sk, dili_param param, dili_rng_func rng, void *rng_ctx);
    int (*verify)(const uint8_t *sig, const uint8_t *msg, int msg_len, const dili_pk *pk, dili_param param);
} test_sig_alg;

extern const test_sig_alg *test_sig_algs[];
extern const size_t test_sig_alg_count;

typedef struct {
    const char *name;
    int security_level;
    kyber_param param;
    size_t pk_bytes;
    size_t sk_bytes;
    size_t ct_bytes;
    size_t msg_bytes;
    size_t seed_bytes;
    size_t coins_bytes;
    int (*check_pk)(const uint8_t *pk, kyber_param param);
    int (*keygen_with_seed)(uint8_t *pk, uint8_t *sk, const uint8_t seed[KYBER_SYMBYTES], kyber_param param);
    int (*encrypt_with_coins)(uint8_t *ct, const uint8_t *pk, const uint8_t msg[KYBER_SYMBYTES], const uint8_t coins[KYBER_SYMBYTES], kyber_param param);
    int (*decrypt)(uint8_t msg[KYBER_SYMBYTES], const uint8_t *sk, const uint8_t *ct, kyber_param param);
} test_enc_alg;

extern const test_enc_alg *test_enc_algs[];
extern const size_t test_enc_alg_count;

typedef struct {
    const uint8_t *msg;
    size_t msg_len;
    const uint8_t *digest;
} test_hash_vec;

typedef struct test_hash_alg {
    const char *name;
    size_t ctx_size;
    size_t digest_size;
    size_t block_size;
    const test_hash_vec *vecs;
    size_t vec_count;
    void (*starts)(void *ctx);
    void (*update)(void *ctx, const uint8_t *input, size_t input_len);
    void (*finish)(void *ctx, uint8_t *output);
    void (*digest)(uint8_t *output, const uint8_t *input, size_t input_len);
} test_hash_alg;

extern const test_hash_alg *test_hash_algs[];
extern const size_t test_hash_alg_count;

typedef struct {
    const char *name;
    int (*run)(void);
} test_pack_suite;

extern const test_pack_suite *test_pack_suites[];
extern const size_t test_pack_suite_count;

#endif
