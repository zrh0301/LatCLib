#ifndef PQC_TEST_RNG_H
#define PQC_TEST_RNG_H

#include <stddef.h>
#include <stdint.h>

/* Reproducible test data only. Never link this helper into the library. */
static uint32_t test_random_state = UINT32_C(0x6d5a56a9);

static inline uint32_t test_random_u32(void) {
    uint32_t x = test_random_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    test_random_state = x;
    return x;
}

static inline void test_random_bytes(uint8_t *out, size_t len) {
    for (size_t i = 0; i < len; ++i)
        out[i] = (uint8_t)test_random_u32();
}

#endif
