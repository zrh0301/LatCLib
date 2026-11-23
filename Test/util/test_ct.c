#include "util/ct.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "ct check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

static uint64_t next_random(void)
{
    static uint64_t state = UINT64_C(0x6374363374363363);
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

int main(void)
{
    uint8_t a[130] = {0}, b[130] = {0}, out[130] = {0}, x[130] = {0};

    /* Equality must agree with memcmp semantics for every supported length,
     * including the n==0 case, on identical and single-bit-different inputs. */
    for (size_t n = 0; n <= sizeof(a); ++n) {
        for (size_t i = 0; i < n; ++i) {
            a[i] = (uint8_t)next_random();
            b[i] = a[i];
        }
        CHECK(pqc_ct_equal(a, b, n) == 1);
        if (n > 0) {
            for (size_t pos = 0; pos < n; ++pos) {
                b[pos] ^= 0x01;
                CHECK(pqc_ct_equal(a, b, n) == 0);
                b[pos] ^= 0x01;
            }
        }
    }

    /* Every bit position in the accumulator must be observed, so a difference
     * confined to bit 7 of one byte still counts. */
    memset(a, 0, sizeof(a));
    memset(b, 0, sizeof(b));
    b[sizeof(a) - 1] = 0x80;
    CHECK(pqc_ct_equal(a, b, sizeof(a)) == 0);

    /* pqc_ct_select keeps out when b==0 and copies x when b==1, and an odd b
     * only uses its low bit. */
    for (size_t n = 0; n <= sizeof(out); ++n) {
        for (size_t i = 0; i < n; ++i) {
            out[i] = (uint8_t)next_random();
            x[i] = (uint8_t)next_random();
        }
        uint8_t saved[130];
        memcpy(saved, out, n);

        pqc_ct_select(out, x, n, 0);
        CHECK(memcmp(out, saved, n) == 0);

        pqc_ct_select(out, x, n, 1);
        CHECK(memcmp(out, x, n) == 0);

        memcpy(out, saved, n);
        pqc_ct_select(out, x, n, 3);
        CHECK(memcmp(out, x, n) == 0);
    }

    /* Selecting in place is the cmov pattern used by KEM implicit rejection. */
    memcpy(out, a, sizeof(out));
    memcpy(x, b, sizeof(x));
    pqc_ct_select(out, x, sizeof(out), 1);
    CHECK(memcmp(out, b, sizeof(out)) == 0);

    puts("ct: constant-time equality and select passed");
    return 0;
}
