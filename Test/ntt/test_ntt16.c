#include "poly/ntt16.h"
#include "poly/poly16.h"

#include <stdio.h>
#include <stdlib.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "ntt16 check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

#define Q 3329
#define DIM 128
#define PSI 17 /* primitive 256th root of unity modulo Q, used by the sampler NTT */

static uint64_t next_random(void)
{
    static uint64_t state = UINT64_C(0x6e747431366e7474);
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static int32_t mod_q(int64_t a)
{
    int64_t r = a % Q;
    return (int32_t)(r < 0 ? r + Q : r);
}

/* Bit-reversal of a 7-bit index, matching the bit-reversed output order. */
static unsigned bit_reverse7(unsigned x)
{
    unsigned r = 0;
    for (unsigned i = 0; i < 7; ++i) {
        r |= ((x >> i) & 1u) << (6 - i);
    }
    return r;
}

/* Direct evaluation of the negacyclic NTT: output slot j uses twiddle
 * exponent 2*bitrev(j)+1 when indexing the primitive 256th root psi. */
static void reference_ntt(const int16_t *a, int32_t *out)
{
    for (unsigned j = 0; j < DIM; ++j) {
        const int64_t e = (int64_t)(2u * bit_reverse7(j) + 1u);
        int64_t acc = 0;
        for (unsigned i = 0; i < DIM; ++i) {
            /* psi^(e*i) with e possibly larger than the group order. */
            int32_t root = 1;
            for (int64_t t = 0; t < (e * i) % 256; ++t) {
                root = mod_q((int64_t)root * PSI);
            }
            acc = mod_q(acc + (int64_t)a[i] * root);
        }
        out[j] = (int32_t)acc;
    }
}

/* Interleave the coefficient blocks as f(x) = sum x^j f_j(x^blocks). */
static int test_product(size_t blocks)
{
    enum { MAX_N = 16 * DIM };
    const size_t n = blocks * DIM;
    int16_t a[MAX_N], b[MAX_N], ta[MAX_N], tb[MAX_N], out[MAX_N];
    int64_t expected[MAX_N] = {0};
    fp16_ctx field;
    CHECK(fp16_init(&field, Q) == 0);

    for (size_t i = 0; i < n; ++i) {
        a[i] = (int16_t)(next_random() % Q);
        b[i] = (int16_t)(next_random() % Q);
        const size_t slot = (i % blocks) * DIM + i / blocks;
        ta[slot] = a[i];
        tb[slot] = b[i];
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const int64_t term = (int64_t)a[i] * b[j];
            if (i + j < n)
                expected[i + j] += term;
            else
                expected[i + j - n] -= term;
        }
    }
    pqc_ntt128_blocks(ta, blocks);
    pqc_ntt128_blocks(tb, blocks);
    CHECK(pqc_ntt128_blocks_mont_mul(out, ta, tb, blocks) == 0);
    CHECK(pqc_ntt128_blocks_mont_mul(ta, ta, tb, blocks) == 0);
    for (size_t i = 0; i < n; ++i)
        CHECK(ta[i] == out[i]);

    fp16_poly_tomont(out, n, &field);
    pqc_invntt128_blocks(out, blocks);
    for (size_t i = 0; i < n; ++i) {
        const size_t slot = (i % blocks) * DIM + i / blocks;
        CHECK(mod_q(out[slot]) == mod_q(expected[i]));
    }
    CHECK(pqc_ntt128_blocks_mont_mul(out, ta, tb, 5) == -1);
    CHECK(mod_q(out[0]) == mod_q(expected[0]));
    return 0;
}

int main(void)
{
    int16_t a[DIM], c[DIM], b[DIM];
    int32_t expect[DIM];

    /* Round trip: invntt(ntt(a)) == a, since the inverse carries 1/128. */
    for (int trial = 0; trial < 64; ++trial) {
        for (int i = 0; i < DIM; ++i) {
            a[i] = (int16_t)(next_random() % Q);
        }
        for (int i = 0; i < DIM; ++i) {
            c[i] = a[i];
        }
        pqc_ntt128(c);
        for (int i = 0; i < DIM; ++i) {
            CHECK(c[i] > -Q && c[i] < Q);
        }
        pqc_invntt128(c);
        for (int i = 0; i < DIM; ++i) {
            CHECK(mod_q(c[i]) == mod_q(a[i]));
        }
    }

    /* Forward transform against the direct definition. */
    for (int trial = 0; trial < 8; ++trial) {
        for (int i = 0; i < DIM; ++i) {
            a[i] = (int16_t)(next_random() % Q);
        }
        reference_ntt(a, expect);
        for (int i = 0; i < DIM; ++i) {
            c[i] = a[i];
        }
        pqc_ntt128(c);
        for (int i = 0; i < DIM; ++i) {
            CHECK(mod_q(c[i]) == expect[i]);
        }
    }

    /* Linearity modulo q, which also checks the block wrapper. */
    for (int i = 0; i < DIM; ++i) {
        a[i] = (int16_t)(next_random() % Q);
        b[i] = (int16_t)(next_random() % Q);
    }
    int16_t sum[DIM], ta[DIM], tb[DIM];
    for (int i = 0; i < DIM; ++i) {
        sum[i] = (int16_t)(a[i] + b[i]);
        ta[i] = a[i];
        tb[i] = b[i];
    }
    pqc_ntt128(sum);
    pqc_ntt128(ta);
    pqc_ntt128(tb);
    for (int i = 0; i < DIM; ++i) {
        CHECK(mod_q(sum[i]) == mod_q((int32_t)ta[i] + tb[i]));
    }

    /* Block wrapper: two independent blocks must not interact. */
    int16_t blocks[2 * DIM], blocks_ref[2 * DIM], blocks_orig[2 * DIM];
    for (int i = 0; i < 2 * DIM; ++i) {
        blocks[i] = (int16_t)(next_random() % Q);
        blocks_orig[i] = blocks[i];
        blocks_ref[i] = blocks[i];
    }
    pqc_ntt128_blocks(blocks, 2);
    for (int k = 0; k < 2; ++k) {
        pqc_ntt128(blocks_ref + k * DIM);
        for (int i = 0; i < DIM; ++i) {
            CHECK(blocks[k * DIM + i] == blocks_ref[k * DIM + i]);
        }
    }
    pqc_invntt128_blocks(blocks, 2);
    for (int i = 0; i < 2 * DIM; ++i) {
        CHECK(mod_q(blocks[i]) == mod_q(blocks_orig[i]));
    }

    CHECK(test_product(4) == 0);
    CHECK(test_product(8) == 0);
    CHECK(test_product(16) == 0);

    puts("ntt16: 128-point negacyclic transform matches the direct definition");
    return 0;
}
