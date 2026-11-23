#include <stdint.h>
#include <string.h>
#include "sm3.h"

#ifndef PUT_UINT32_BE
/* Put 32-bit integer(n) to big endian bytes(b) */
#define PUT_UINT32_BE(n, b, i)                        \
    {                                                 \
        (b)[(i)] = (uint8_t)(((n) >> 24) & 0xFF);     \
        (b)[(i) + 1] = (uint8_t)(((n) >> 16) & 0xFF); \
        (b)[(i) + 2] = (uint8_t)(((n) >> 8) & 0xFF);  \
        (b)[(i) + 3] = (uint8_t)(((n)) & 0xFF);       \
    }
#endif

#ifndef GET_UINT32_BE
/* Get 32-bit integer(n) from big endian bytes(b) */
#define GET_UINT32_BE(n, b, i)                                                                                                        \
    {                                                                                                                                 \
        (n) = ((uint32_t)(b)[(i)] << 24) | ((uint32_t)(b)[(i) + 1] << 16) | ((uint32_t)(b)[(i) + 2] << 8) | ((uint32_t)(b)[(i) + 3]); \
    }
#endif

#ifndef ROTL32
/* Rotate 32-bit integer(x) to left by (n) bits */
#define ROTL32(x, n) ((((uint32_t)x) << (n)) | ((x) >> (32 - (n))))
#endif

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* Hash a single 512-bit block. This is the core of the algorithm. */
static void _sm3_process(uint32_t state[8], const uint8_t buffer[64])
{
    uint32_t a, b, c, d, e, f, g, h;
    uint32_t t1, t2;
    uint32_t W[68];

    /* Copy ctx->state[] to working vars */
    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];
    f = state[5];
    g = state[6];
    h = state[7];

    /* get W0~W1  */
    GET_UINT32_BE(W[0], buffer, 0);
    GET_UINT32_BE(W[1], buffer, 4);
    GET_UINT32_BE(W[2], buffer, 8);
    GET_UINT32_BE(W[3], buffer, 12);

#define P0(x) (x ^ ROTL32(x, 9) ^ ROTL32(x, 17))
#define P1(x) (x ^ ROTL32(x, 15) ^ ROTL32(x, 23))
#define FF1(a, b, c) (a ^ b ^ c)
#define FF2(a, b, c) ((a & b) | (a & c) | (b & c))
#define GG1(e, f, g) (e ^ f ^ g)
#define GG2(e, f, g) ((e & f) | ((~e) & g))

#define ROUND_FUNC_0(a, b, c, d, e, f, g, h, j)                       \
    {                                                                 \
        t1 = ROTL32(a, 12) + e + ROTL32(0x79cc4519, j);               \
        t1 = ROTL32(t1, 7);                                           \
        GET_UINT32_BE(W[j + 4], buffer, j * 4 + 16);                  \
        d += FF1(a, b, c) + (t1 ^ ROTL32(a, 12)) + (W[j] ^ W[j + 4]); \
        h += GG1(e, f, g) + t1 + W[j];                                \
        b = ROTL32(b, 9);                                             \
        f = ROTL32(f, 19);                                            \
        h = P0(h);                                                    \
    }

#define ROUND_FUNC_0_0(a, b, c, d, e, f, g, h, j)                     \
    {                                                                 \
        t1 = ROTL32(a, 12) + e + 0x79cc4519;                          \
        t1 = ROTL32(t1, 7);                                           \
        GET_UINT32_BE(W[j + 4], buffer, j * 4 + 16);                  \
        d += FF1(a, b, c) + (t1 ^ ROTL32(a, 12)) + (W[j] ^ W[j + 4]); \
        h += GG1(e, f, g) + t1 + W[j];                                \
        b = ROTL32(b, 9);                                             \
        f = ROTL32(f, 19);                                            \
        h = P0(h);                                                    \
    }

#define ROUND_FUNC_1(a, b, c, d, e, f, g, h, j)                       \
    {                                                                 \
        t1 = ROTL32(a, 12) + e + ROTL32(0x79cc4519, j);               \
        t1 = ROTL32(t1, 7);                                           \
        t2 = ROTL32(W[j + 1], 15) ^ W[j - 5] ^ W[j - 12];             \
        W[j + 4] = P1(t2) ^ ROTL32(W[j - 9], 7) ^ W[j - 2];           \
        d += FF1(a, b, c) + (t1 ^ ROTL32(a, 12)) + (W[j] ^ W[j + 4]); \
        h += GG1(e, f, g) + t1 + W[j];                                \
        b = ROTL32(b, 9);                                             \
        f = ROTL32(f, 19);                                            \
        h = P0(h);                                                    \
    }

#define ROUND_FUNC_2(a, b, c, d, e, f, g, h, j)                       \
    {                                                                 \
        t1 = ROTL32(a, 12) + e + ROTL32(0x7a879d8a, j & 0x1F);        \
        t1 = ROTL32(t1, 7);                                           \
        t2 = ROTL32(W[j + 1], 15) ^ W[j - 5] ^ W[j - 12];             \
        W[j + 4] = P1(t2) ^ ROTL32(W[j - 9], 7) ^ W[j - 2];           \
        d += FF2(a, b, c) + (t1 ^ ROTL32(a, 12)) + (W[j] ^ W[j + 4]); \
        h += GG2(e, f, g) + t1 + W[j];                                \
        b = ROTL32(b, 9);                                             \
        f = ROTL32(f, 19);                                            \
        h = P0(h);                                                    \
    }

#define ROUND_FUNC_2_32(a, b, c, d, e, f, g, h, j)                    \
    {                                                                 \
        t1 = ROTL32(a, 12) + e + 0x7a879d8a;                          \
        t1 = ROTL32(t1, 7);                                           \
        t2 = ROTL32(W[j + 1], 15) ^ W[j - 5] ^ W[j - 12];             \
        W[j + 4] = P1(t2) ^ ROTL32(W[j - 9], 7) ^ W[j - 2];           \
        d += FF2(a, b, c) + (t1 ^ ROTL32(a, 12)) + (W[j] ^ W[j + 4]); \
        h += GG2(e, f, g) + t1 + W[j];                                \
        b = ROTL32(b, 9);                                             \
        f = ROTL32(f, 19);                                            \
        h = P0(h);                                                    \
    }

    ROUND_FUNC_0_0(a, b, c, d, e, f, g, h, 0);
    ROUND_FUNC_0(d, a, b, c, h, e, f, g, 1);
    ROUND_FUNC_0(c, d, a, b, g, h, e, f, 2);
    ROUND_FUNC_0(b, c, d, a, f, g, h, e, 3);
    ROUND_FUNC_0(a, b, c, d, e, f, g, h, 4);
    ROUND_FUNC_0(d, a, b, c, h, e, f, g, 5);
    ROUND_FUNC_0(c, d, a, b, g, h, e, f, 6);
    ROUND_FUNC_0(b, c, d, a, f, g, h, e, 7);
    ROUND_FUNC_0(a, b, c, d, e, f, g, h, 8);
    ROUND_FUNC_0(d, a, b, c, h, e, f, g, 9);
    ROUND_FUNC_0(c, d, a, b, g, h, e, f, 10);
    ROUND_FUNC_0(b, c, d, a, f, g, h, e, 11);

    ROUND_FUNC_1(a, b, c, d, e, f, g, h, 12);
    ROUND_FUNC_1(d, a, b, c, h, e, f, g, 13);
    ROUND_FUNC_1(c, d, a, b, g, h, e, f, 14);
    ROUND_FUNC_1(b, c, d, a, f, g, h, e, 15);

    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 16);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 17);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 18);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 19);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 20);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 21);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 22);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 23);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 24);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 25);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 26);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 27);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 28);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 29);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 30);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 31);
    ROUND_FUNC_2_32(a, b, c, d, e, f, g, h, 32);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 33);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 34);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 35);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 36);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 37);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 38);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 39);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 40);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 41);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 42);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 43);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 44);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 45);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 46);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 47);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 48);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 49);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 50);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 51);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 52);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 53);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 54);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 55);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 56);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 57);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 58);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 59);
    ROUND_FUNC_2(a, b, c, d, e, f, g, h, 60);
    ROUND_FUNC_2(d, a, b, c, h, e, f, g, 61);
    ROUND_FUNC_2(c, d, a, b, g, h, e, f, 62);
    ROUND_FUNC_2(b, c, d, a, f, g, h, e, 63);

    /* Add the working vars back into ctx.state[] */
    state[0] ^= a;
    state[1] ^= b;
    state[2] ^= c;
    state[3] ^= d;
    state[4] ^= e;
    state[5] ^= f;
    state[6] ^= g;
    state[7] ^= h;
    /* Wipe variables */
    //   a = b = c = d = e = f = g = h =0;

#undef ROUND_FUNC_0
#undef ROUND_FUNC_0_0
#undef ROUND_FUNC_1
#undef ROUND_FUNC_2
#undef ROUND_FUNC_2_32

#undef P0
#undef P1
#undef FF1
#undef FF2
#undef GG1
#undef GG2
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

int sm3_init(sm3_ctx *md)
{
    /* Initialization constants */
    md->bytes[0] = 0;
    md->bytes[1] = 0;
    md->state[0] = 0x7380166f;
    md->state[1] = 0x4914b2b9;
    md->state[2] = 0x172442d7;
    md->state[3] = 0xda8a0600;
    md->state[4] = 0xa96f30bc;
    md->state[5] = 0x163138aa;
    md->state[6] = 0xe38dee4d;
    md->state[7] = 0xb0fb0e4e;
    // memset(md->buffer, 0, 64);
    return 1;
}

int sm3_update(sm3_ctx *md, const void *input, size_t ilen)
{
    const uint8_t *in = (const uint8_t *)input;
    uint32_t left, fill;

    if (in == NULL || ilen == 0)
        return 1;

    left = md->bytes[0] & 0x3F;
    fill = 64 - left;

    md->bytes[0] += (uint32_t)ilen;
    if (md->bytes[0] < (uint32_t)ilen)
        md->bytes[1]++;

    if (left && ilen >= fill)
    {
        memcpy(md->buffer + left, in, fill);
        _sm3_process(md->state, md->buffer);
        in += fill;
        ilen -= fill;
        left = 0;
    }
    while (ilen >= 64)
    {
        _sm3_process(md->state, in);
        in += 64;
        ilen -= 64;
    }
    if (ilen > 0)
    {
        memcpy(md->buffer + left, in, ilen);
    }

    return 1;
}

static const uint8_t sm3_padding[64] = {
    0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int sm3_final(uint8_t md[SM3_HASH_SIZE], sm3_ctx *ctx)
{
    uint32_t last, padn;
    uint32_t lenh, lenl;
    uint8_t msglen[8];
    size_t i;

    lenh = (ctx->bytes[0] >> 29) | (ctx->bytes[1] << 3);
    lenl = (ctx->bytes[0] << 3);
    PUT_UINT32_BE(lenh, msglen, 0);
    PUT_UINT32_BE(lenl, msglen, 4);

    last = ctx->bytes[0] & 0x3F;
    padn = (last < 56) ? (56 - last) : (120 - last);

    (void)sm3_update(ctx, sm3_padding, padn);
    (void)sm3_update(ctx, msglen, 8);

    for (i = 0; i < 8; i++)
        PUT_UINT32_BE(ctx->state[i], md, i * 4);

    return 1;
}

/**
 * @brief          Output = SM3( input buffer )
 *
 * @param  md      SM3 hash result
 * @param  input   buffer holding the input data
 * @param  ilen    length of the input data in bytes
 */
void sm3_hash(uint8_t md[SM3_HASH_SIZE], const uint8_t *input, size_t ilen)
{
    sm3_ctx ctx;
    (void)sm3_init(&ctx);
    (void)sm3_update(&ctx, input, ilen);
    (void)sm3_final(md, &ctx);
}
