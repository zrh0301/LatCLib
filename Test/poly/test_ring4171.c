#include "poly/ring4171.h"
#include "poly/fp32.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "ring4171 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

#define Q PQC_Q4171
#define N PQC_Q4171_N

static uint32_t next_random(void)
{
    static uint32_t state = UINT32_C(0x2545f491);
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static int congruent(int64_t a, int64_t b)
{
    int64_t d = a - b;
    d %= Q;
    return d == 0;
}

int main(void)
{
    /* Montgomery: reduce(x * R) = x for the field representative x. */
    for (int i = 0; i < 200; ++i) {
        int32_t x = (int32_t)(next_random() % (uint32_t)Q);
        int32_t r = pqc_q4171_montgomery_reduce((int64_t)x * 2208763);
        CHECK(congruent(r, x));
        CHECK(r > -Q && r < Q);
    }

    /* Barrett and general reduction agree with the modulus, in range.
     * Barrett requires |a| < 460138740; general covers the int32 range. */
    for (int i = 0; i < 200; ++i) {
        int32_t a = (int32_t)(next_random() % 460138740u) -
                    (int32_t)(460138740u / 2);
        int32_t g = (int32_t)next_random();
        int32_t rb = pqc_q4171_barrett_reduce(a);
        int32_t rg = pqc_q4171_general_reduce(g);
        CHECK(congruent(rb, a));
        CHECK(rb > -Q && rb < Q);
        CHECK(congruent(rg, g));
        CHECK(rg > -Q && rg < Q);
    }

    /* Canonical and centred representatives. */
    for (int i = 0; i < 200; ++i) {
        int32_t a = (int32_t)((next_random() % (uint32_t)(2 * Q - 1)) - (Q - 1));
        int32_t ra = pqc_q4171_amod(a);
        int32_t rc = pqc_q4171_cmod(a);
        CHECK(ra >= 0 && ra < Q);
        CHECK(congruent(ra, a));
        CHECK(rc > -Q / 2 && rc <= Q / 2);
        CHECK(congruent(rc, a));
    }

    /* NTT round trip returns R * a (mod q), as documented for the ring. */
    static int32_t a[N], b[N], work[N], conv[N], prod[N];
    for (int i = 0; i < N; ++i) {
        a[i] = (int32_t)(next_random() % (uint32_t)Q);
        b[i] = (int32_t)(next_random() % (uint32_t)Q);
    }

    memcpy(work, a, sizeof(a));
    pqc_q4171_ntt(work);
    pqc_q4171_invntt_frominvmont(work);
    for (int i = 0; i < N; ++i) {
        int32_t scaled = pqc_q4171_barrett_reduce(
            (int32_t)fp32_reduce_raw((int64_t)a[i] * 2208763, Q));
        CHECK(congruent(work[i], scaled));
    }

    /* NTT multiplication equals the schoolbook negacyclic convolution. */
    for (int k = 0; k < N; ++k) {
        int64_t sum = 0;
        for (int i = 0; i < N; ++i) {
            int j = k - i;
            if (j >= 0)
                sum += (int64_t)a[i] * b[j];
            else
                sum -= (int64_t)a[i] * b[j + N];
        }
        conv[k] = (int32_t)fp32_reduce_raw(sum, Q);
    }

    static int32_t an[N], bn[N];
    memcpy(an, a, sizeof(a));
    memcpy(bn, b, sizeof(b));
    pqc_q4171_ntt(an);
    pqc_q4171_ntt(bn);
    pqc_q4171_pointwise_invmontgomery(prod, an, bn, N);
    pqc_q4171_invntt_frominvmont(prod);
    for (int i = 0; i < N; ++i) {
        int32_t reduced = pqc_q4171_amod(prod[i]);
        CHECK(reduced == conv[i]);
    }

    /* Coefficient-wise helpers. */
    static int32_t c1[N], c2[N];
    pqc_q4171_add(c1, a, b, N);
    pqc_q4171_sub(c2, a, b, N);
    for (int i = 0; i < N; ++i) {
        CHECK(c1[i] == a[i] + b[i]);
        CHECK(c2[i] == a[i] - b[i]);
    }
    pqc_q4171_sub_scaled(c1, a, b, 7, N);
    for (int i = 0; i < N; ++i)
        CHECK(c1[i] == a[i] - 7 * b[i]);
    memcpy(c1, a, sizeof(a));
    pqc_q4171_shiftl(c1, N, 3);
    for (int i = 0; i < N; ++i)
        CHECK(c1[i] == a[i] * 8);

    /* chknorm: 1 exactly when some centred coefficient leaves the bound. */
    for (int i = 0; i < N; ++i)
        c1[i] = i % 11 - 5;
    CHECK(pqc_q4171_chknorm(c1, N, 6) == 0);
    CHECK(pqc_q4171_chknorm(c1, N, 5) == 1);
    CHECK(pqc_q4171_chknorm(c1, N, 4) == 1);

    puts("ring4171: reductions, NTT and coefficient helpers passed");
    return 0;
}
