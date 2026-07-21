#ifndef PQC_KECCAK_H
#define PQC_KECCAK_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
    uint64_t A[5][5];       // Keccak 状态矩阵
    size_t block_size;      // 率 (Rate)
    size_t md_size;         // 输出长度 (对于 SHAKE 可忽略)
    size_t bufsz;           // 缓冲区当前使用量
    size_t squeeze_off;     // XOF 挤压偏移: 当前已从置换态输出的字节数
    unsigned char pad;      // 填充字符 (0x06 或 0x1F)
    unsigned char buf[200]; // 临时缓冲区 (Keccak 最大 rate 为 1600/8 = 200)
} KECCAK1600_CTX;

/* 基础 API */
void sha3_reset(KECCAK1600_CTX *ctx);
int sha3_init(KECCAK1600_CTX *ctx, unsigned char pad, size_t bitlen);
int sha3_update(KECCAK1600_CTX *ctx, const void *_inp, size_t len);
int sha3_final(unsigned char *md, KECCAK1600_CTX *ctx);

/* PQC 常用 API (SHAKE XOF) */
void shake128_init(KECCAK1600_CTX *ctx);
void shake256_init(KECCAK1600_CTX *ctx);
void shake_xof_final(KECCAK1600_CTX *ctx);
void shake_squeeze(KECCAK1600_CTX *ctx, unsigned char *out, size_t len);

/* 一次性调用简便接口 */
void shake128(unsigned char *output, size_t outlen, const unsigned char *input, size_t inlen);
void shake256(unsigned char *output, size_t outlen, const unsigned char *input, size_t inlen);

#ifdef __cplusplus
}
#endif

#endif // PQC_KECCAK_H