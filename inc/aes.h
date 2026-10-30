#ifndef PQC_AES_H
#define PQC_AES_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* AES 密钥与分组长度 (字节) */
#define PQC_AES128_KEYBYTES 16
#define PQC_AES256_KEYBYTES 32
#define PQC_AES_BLOCKBYTES 16

/* AES-128 ECB 加密: 对 input 起始的 num 字节 (16 的整数倍) 用 key (16 字节)
 * 逐块独立加密, 结果写入 output。 */
void AES128_ECB(const unsigned char *key, const uint8_t *input, unsigned char *output, size_t num);

/* AES-256 ECB 加密单个 16 字节块: 将 key (32 字节) 下的明文 ctr (16 字节)
 * 加密后写入 buffer。NIST DRBG 与 seed expander 以此按计数器生成数据流。 */
void AES256_ECB(const unsigned char *key, const unsigned char *ctr, unsigned char *buffer);

/* AES-128 CTR 密钥流: 以 input (inputByteLen 字节) 为密钥, 全零 16 字节 IV,
 * 块计数器从 0 开始递增, 生成的密钥流直接覆盖 output 的前 outputByteLen
 * 字节 (与明文无关)。返回值恒为 outputByteLen。 */
int AES128_CTR(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen);

#ifdef __cplusplus
}
#endif

#endif /* PQC_AES_H */
