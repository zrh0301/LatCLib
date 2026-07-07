#include "../inc/keccak.h"
#include <assert.h>
#include <stdio.h>

/* ================= 配置与宏定义 ================= */

#define KECCAK_2X

#define SHA3_BLOCKSIZE(bitlen) (1600 - (bitlen) * 2) / 8
#define OSSL_NELEM(x) (sizeof(x)/sizeof((x)[0]))

// 循环左移
static uint64_t ROL64(uint64_t val, int offset) {
    return (offset == 0) ? val : ((val << offset) | (val >> (64 - offset)));
}

// 常量表
static const unsigned char rhotates[5][5] = {
    { 0, 1, 62, 28, 27 }, { 36, 44, 6, 55, 20 },
    { 3, 10, 43, 25, 39 }, { 41, 45, 15, 21, 8 },
    { 18, 2, 61, 56, 14 }
};

static const uint64_t iotas[] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
    0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
    0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
    0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
    0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

/* ================= 核心置换函数 (KECCAK_2X) ================= */

static void Round(uint64_t R[5][5], uint64_t A[5][5], size_t i) {
    uint64_t C[5], D[5];

    C[0] = A[0][0] ^ A[1][0] ^ A[2][0] ^ A[3][0] ^ A[4][0];
    C[1] = A[0][1] ^ A[1][1] ^ A[2][1] ^ A[3][1] ^ A[4][1];
    C[2] = A[0][2] ^ A[1][2] ^ A[2][2] ^ A[3][2] ^ A[4][2];
    C[3] = A[0][3] ^ A[1][3] ^ A[2][3] ^ A[3][3] ^ A[4][3];
    C[4] = A[0][4] ^ A[1][4] ^ A[2][4] ^ A[3][4] ^ A[4][4];

    D[0] = ROL64(C[1], 1) ^ C[4];
    D[1] = ROL64(C[2], 1) ^ C[0];
    D[2] = ROL64(C[3], 1) ^ C[1];
    D[3] = ROL64(C[4], 1) ^ C[2];
    D[4] = ROL64(C[0], 1) ^ C[3];

    C[0] = A[0][0] ^ D[0];
    C[1] = ROL64(A[1][1] ^ D[1], rhotates[1][1]);
    C[2] = ROL64(A[2][2] ^ D[2], rhotates[2][2]);
    C[3] = ROL64(A[3][3] ^ D[3], rhotates[3][3]);
    C[4] = ROL64(A[4][4] ^ D[4], rhotates[4][4]);

    R[0][0] = C[0] ^ (~C[1] & C[2]) ^ iotas[i];
    R[0][1] = C[1] ^ (~C[2] & C[3]);
    R[0][2] = C[2] ^ (~C[3] & C[4]);
    R[0][3] = C[3] ^ (~C[4] & C[0]);
    R[0][4] = C[4] ^ (~C[0] & C[1]);

    C[0] = ROL64(A[0][3] ^ D[3], rhotates[0][3]);
    C[1] = ROL64(A[1][4] ^ D[4], rhotates[1][4]);
    C[2] = ROL64(A[2][0] ^ D[0], rhotates[2][0]);
    C[3] = ROL64(A[3][1] ^ D[1], rhotates[3][1]);
    C[4] = ROL64(A[4][2] ^ D[2], rhotates[4][2]);

    R[1][0] = C[0] ^ (~C[1] & C[2]);
    R[1][1] = C[1] ^ (~C[2] & C[3]);
    R[1][2] = C[2] ^ (~C[3] & C[4]);
    R[1][3] = C[3] ^ (~C[4] & C[0]);
    R[1][4] = C[4] ^ (~C[0] & C[1]);

    C[0] = ROL64(A[0][1] ^ D[1], rhotates[0][1]);
    C[1] = ROL64(A[1][2] ^ D[2], rhotates[1][2]);
    C[2] = ROL64(A[2][3] ^ D[3], rhotates[2][3]);
    C[3] = ROL64(A[3][4] ^ D[4], rhotates[3][4]);
    C[4] = ROL64(A[4][0] ^ D[0], rhotates[4][0]);

    R[2][0] = C[0] ^ (~C[1] & C[2]);
    R[2][1] = C[1] ^ (~C[2] & C[3]);
    R[2][2] = C[2] ^ (~C[3] & C[4]);
    R[2][3] = C[3] ^ (~C[4] & C[0]);
    R[2][4] = C[4] ^ (~C[0] & C[1]);

    C[0] = ROL64(A[0][4] ^ D[4], rhotates[0][4]);
    C[1] = ROL64(A[1][0] ^ D[0], rhotates[1][0]);
    C[2] = ROL64(A[2][1] ^ D[1], rhotates[2][1]);
    C[3] = ROL64(A[3][2] ^ D[2], rhotates[3][2]);
    C[4] = ROL64(A[4][3] ^ D[3], rhotates[4][3]);

    R[3][0] = C[0] ^ (~C[1] & C[2]);
    R[3][1] = C[1] ^ (~C[2] & C[3]);
    R[3][2] = C[2] ^ (~C[3] & C[4]);
    R[3][3] = C[3] ^ (~C[4] & C[0]);
    R[3][4] = C[4] ^ (~C[0] & C[1]);

    C[0] = ROL64(A[0][2] ^ D[2], rhotates[0][2]);
    C[1] = ROL64(A[1][3] ^ D[3], rhotates[1][3]);
    C[2] = ROL64(A[2][4] ^ D[4], rhotates[2][4]);
    C[3] = ROL64(A[3][0] ^ D[0], rhotates[3][0]);
    C[4] = ROL64(A[4][1] ^ D[1], rhotates[4][1]);

    R[4][0] = C[0] ^ (~C[1] & C[2]);
    R[4][1] = C[1] ^ (~C[2] & C[3]);
    R[4][2] = C[2] ^ (~C[3] & C[4]);
    R[4][3] = C[3] ^ (~C[4] & C[0]);
    R[4][4] = C[4] ^ (~C[0] & C[1]);
}

static void KeccakF1600(uint64_t A[5][5]) {
    uint64_t T[5][5];
    size_t i;
    for (i = 0; i < 24; i += 2) {
        Round(T, A, i);
        Round(A, T, i + 1);
    }
}

/* ================= 吸收与挤压 ================= */

static size_t SHA3_absorb(uint64_t A[5][5], const unsigned char *inp, size_t len, size_t r) {
    uint64_t *A_flat = (uint64_t *)A;
    size_t i, w = r / 8;

    while (len >= r) {
        for (i = 0; i < w; i++) {
            uint64_t Ai = (uint64_t)inp[0] | (uint64_t)inp[1] << 8 | 
                          (uint64_t)inp[2] << 16 | (uint64_t)inp[3] << 24 | 
                          (uint64_t)inp[4] << 32 | (uint64_t)inp[5] << 40 | 
                          (uint64_t)inp[6] << 48 | (uint64_t)inp[7] << 56;
            inp += 8;
            A_flat[i] ^= Ai; // 假设大端/小端不需要互换(BitInterleave)
        }
        KeccakF1600(A);
        len -= r;
    }
    return len;
}

static void internal_SHA3_squeeze(uint64_t A[5][5], unsigned char *out, size_t len, size_t r) {
    uint64_t *A_flat = (uint64_t *)A;
    size_t i, w = r / 8;
    int next = 0; // 第一次调用不置换，因为 absorb 结束时已置换过

    while (len != 0) {
        if (next) KeccakF1600(A);
        next = 1;
        for (i = 0; i < w && len != 0; i++) {
            uint64_t Ai = A_flat[i];
            if (len < 8) {
                for (size_t j = 0; j < len; j++) {
                    *out++ = (unsigned char)Ai;
                    Ai >>= 8;
                }
                return;
            }
            out[0] = (unsigned char)(Ai);
            out[1] = (unsigned char)(Ai >> 8);
            out[2] = (unsigned char)(Ai >> 16);
            out[3] = (unsigned char)(Ai >> 24);
            out[4] = (unsigned char)(Ai >> 32);
            out[5] = (unsigned char)(Ai >> 40);
            out[6] = (unsigned char)(Ai >> 48);
            out[7] = (unsigned char)(Ai >> 56);
            out += 8;
            len -= 8;
        }
    }
}

/* ================= 外部接口实现 ================= */

void sha3_reset(KECCAK1600_CTX *ctx) {
    memset(ctx->A, 0, sizeof(ctx->A));
    ctx->bufsz = 0;
    ctx->squeeze_off = 0;
}

int sha3_init(KECCAK1600_CTX *ctx, unsigned char pad, size_t bitlen) {
    size_t bsz = SHA3_BLOCKSIZE(bitlen);
    if (bsz <= sizeof(ctx->buf)) {
        sha3_reset(ctx);
        ctx->block_size = bsz;
        ctx->md_size = bitlen / 8;
        ctx->pad = pad;
        return 1;
    }
    return 0;
}

int sha3_update(KECCAK1600_CTX *ctx, const void *_inp, size_t len) {
    const unsigned char *inp = _inp;
    size_t bsz = ctx->block_size;
    size_t num, rem;

    if (len == 0) return 1;

    if ((num = ctx->bufsz) != 0) {
        rem = bsz - num;
        if (len < rem) {
            memcpy(ctx->buf + num, inp, len);
            ctx->bufsz += len;
            return 1;
        }
        memcpy(ctx->buf + num, inp, rem);
        inp += rem; len -= rem;
        (void)SHA3_absorb(ctx->A, ctx->buf, bsz, bsz);
        ctx->bufsz = 0;
    }

    if (len >= bsz)
        rem = SHA3_absorb(ctx->A, inp, len, bsz);
    else
        rem = len;

    if (rem) {
        memcpy(ctx->buf, inp + len - rem, rem);
        ctx->bufsz = rem;
    }
    return 1;
}

int sha3_final(unsigned char *md, KECCAK1600_CTX *ctx) {
    size_t bsz = ctx->block_size;
    size_t num = ctx->bufsz;

    if (ctx->md_size == 0) return 1;

    memset(ctx->buf + num, 0, bsz - num);
    ctx->buf[num] = ctx->pad;
    ctx->buf[bsz - 1] |= 0x80;

    (void)SHA3_absorb(ctx->A, ctx->buf, bsz, bsz);
    internal_SHA3_squeeze(ctx->A, md, ctx->md_size, bsz);
    return 1;
}

/* SHAKE / PQC 专用接口 */
void shake128_init(KECCAK1600_CTX *ctx) { sha3_init(ctx, 0x1F, 128); }
void shake256_init(KECCAK1600_CTX *ctx) { sha3_init(ctx, 0x1F, 256); }

void shake_xof_final(KECCAK1600_CTX *ctx) {
    size_t bsz = ctx->block_size;
    size_t num = ctx->bufsz;
    memset(ctx->buf + num, 0, bsz - num);
    ctx->buf[num] = ctx->pad;
    ctx->buf[bsz - 1] |= 0x80;
    (void)SHA3_absorb(ctx->A, ctx->buf, bsz, bsz);
    ctx->bufsz = 0;
    ctx->squeeze_off = 0;   // 状态已置换, 从块首开始挤压
}

// 连续 XOF 挤压: 跨多次调用维持流位置 (x86 小端下直接按字节读状态)
void shake_squeeze(KECCAK1600_CTX *ctx, unsigned char *out, size_t len) {
    const size_t r = ctx->block_size;
    const unsigned char *A8 = (const unsigned char *)ctx->A;
    while (len) {
        if (ctx->squeeze_off == r) {
            KeccakF1600(ctx->A);
            ctx->squeeze_off = 0;
        }
        size_t take = r - ctx->squeeze_off;
        if (take > len) take = len;
        memcpy(out, A8 + ctx->squeeze_off, take);
        out += take;
        len -= take;
        ctx->squeeze_off += take;
    }
}

void shake128(unsigned char *output, size_t outlen, const unsigned char *input, size_t inlen) {
    KECCAK1600_CTX ctx;
    shake128_init(&ctx);
    sha3_update(&ctx, input, inlen);
    shake_xof_final(&ctx);
    shake_squeeze(&ctx, output, outlen);
}

void shake256(unsigned char *output, size_t outlen, const unsigned char *input, size_t inlen) {
    KECCAK1600_CTX ctx;
    shake256_init(&ctx);
    sha3_update(&ctx, input, inlen);
    shake_xof_final(&ctx);
    shake_squeeze(&ctx, output, outlen);
}