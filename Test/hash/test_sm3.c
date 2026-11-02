/* SM3 streaming/boundary self-consistency suite.
 *
 * Absolute digests live in Test/hash/cases.c (cross-verified KATs); this
 * executable checks the properties a Monolithic KAT cannot reach:
 *   - every streaming split (1/3/7/13/31/63/64/200-byte chunks) reproduces
 *     the one-shot digest for every message length 0..270;
 *   - the MD padding branches at 55/56 and 63/64/119/120/127/128 agree;
 *   - a reused context is fully reset by sm3_init();
 *   - zero-length updates and NULL/0 inputs are no-ops.
 */

#include <stdio.h>
#include <string.h>
#include "sm3.h"

static const uint8_t sm3_abc_digest[SM3_HASH_SIZE] = {
    0x66, 0xc7, 0xf0, 0xf4, 0x62, 0xee, 0xed, 0xd9,
    0xd1, 0xf2, 0xd4, 0x6b, 0xdc, 0x10, 0xe4, 0xe2,
    0x41, 0x67, 0xc4, 0x87, 0x5c, 0xf2, 0xf7, 0xa2,
    0x29, 0x7d, 0xa0, 0x2b, 0x8f, 0x4b, 0xa8, 0xe0,
};

static int run_sweep(void) {
    /* padding: buffer spans 3 blocks + 1 byte so the absorb tail logic and
     * both padn branches run against unaligned message tails. */
    uint8_t msg[270], ref[SM3_HASH_SIZE], out[SM3_HASH_SIZE];
    sm3_ctx ctx;
    static const size_t chunk_sizes[] = {1, 3, 7, 13, 31, 63, 64, 200};
    int ok = 1;

    for (size_t i = 0; i < sizeof(msg); i++) {
        msg[i] = (uint8_t)(i * 7 + 11);
    }

    for (size_t len = 0; len <= sizeof(msg); len++) {
        sm3_hash(ref, msg, len);
        for (size_t c = 0; c < sizeof(chunk_sizes) / sizeof(chunk_sizes[0]); c++) {
            size_t step = chunk_sizes[c];
            (void)sm3_init(&ctx);
            for (size_t off = 0; off < len; off += step) {
                size_t n = (len - off < step) ? (len - off) : step;
                (void)sm3_update(&ctx, msg + off, n);
            }
            (void)sm3_update(&ctx, NULL, 0); /* zero-length update is a no-op */
            (void)sm3_final(out, &ctx);
            if (memcmp(out, ref, SM3_HASH_SIZE) != 0) {
                printf("  [FAIL] chunk=%zu len=%zu streaming digest mismatch\n",
                       step, len);
                ok = 0;
                goto done;
            }
        }
    }
done:
    return ok;
}

static int run_boundaries(void) {
    /* explicit MD padding edge lengths: last<56 (single pad block) and
     * last>=56 (pad spills to a second block) on both sides */
    static const size_t edge_lens[] = {
        55, 56, 57, 63, 64, 65, 119, 120, 121, 127, 128, 129,
    };
    uint8_t msg[270], ref[SM3_HASH_SIZE], out[SM3_HASH_SIZE];
    sm3_ctx ctx;
    int ok = 1;

    for (size_t i = 0; i < sizeof(msg); i++) {
        msg[i] = (uint8_t)(i * 0x1f + 3);
    }

    for (size_t v = 0; v < sizeof(edge_lens) / sizeof(edge_lens[0]); v++) {
        size_t len = edge_lens[v];
        sm3_hash(ref, msg, len);
        (void)sm3_init(&ctx);
        (void)sm3_update(&ctx, msg, len);
        (void)sm3_final(out, &ctx);
        if (memcmp(out, ref, SM3_HASH_SIZE) != 0) {
            printf("  [FAIL] streaming digest mismatch at length %zu\n", len);
            ok = 0;
        }
    }
    return ok;
}

static int run_reuse(void) {
    uint8_t first[SM3_HASH_SIZE], second[SM3_HASH_SIZE];
    sm3_ctx ctx;

    (void)sm3_init(&ctx);
    (void)sm3_update(&ctx, (const uint8_t *)"abc", 3);
    (void)sm3_final(first, &ctx);

    /* re-init on the same context must behave like a fresh one */
    (void)sm3_init(&ctx);
    (void)sm3_update(&ctx, (const uint8_t *)"abc", 3);
    (void)sm3_final(second, &ctx);

    return memcmp(first, sm3_abc_digest, SM3_HASH_SIZE) == 0 &&
           memcmp(second, sm3_abc_digest, SM3_HASH_SIZE) == 0;
}

static int run_zero_and_null(void) {
    uint8_t out[SM3_HASH_SIZE], ref[SM3_HASH_SIZE];
    sm3_ctx ctx;

    sm3_hash(ref, (const uint8_t *)"abc", 3);

    /* zero-length updates must not disturb the state */
    (void)sm3_init(&ctx);
    (void)sm3_update(&ctx, NULL, 0);
    (void)sm3_update(&ctx, (const uint8_t *)"abcd", 3); /* only abcd..abc: 3 bytes */
    (void)sm3_final(out, &ctx);
    if (memcmp(out, ref, SM3_HASH_SIZE) != 0) {
        printf("  [FAIL] zero-length updates changed the digest\n");
        return 0;
    }

    /* one-shot NULL/0 must equal the empty-string digest */
    sm3_hash(out, NULL, 0);
    sm3_hash(ref, (const uint8_t *)"", 0);
    return memcmp(out, ref, SM3_HASH_SIZE) == 0;
}

int main(void) {
    struct {
        const char *name;
        int (*run)(void);
    } cases[] = {
        {"sweep 0..270 all chunk splits", run_sweep},
        {"MD padding boundary lengths", run_boundaries},
        {"context reuse re-init", run_reuse},
        {"zero-length and NULL edges", run_zero_and_null},
    };
    int ok = 1;

    printf("============================================\n");
    printf("  Test Suite: SM3 Streaming Boundaries\n");
    printf("============================================\n");

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        printf("  %-34s", cases[i].name);
        if (cases[i].run()) {
            printf("[PASS]\n");
        } else {
            printf("[FAIL]\n");
            ok = 0;
        }
    }

    return ok ? 0 : 1;
}
