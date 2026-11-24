/* Scloud+ 参考实现要求的 AES-128-CTR 行展开接口。
 *
 * 参考包自带一份可移植 AES（ref/aes_reference.c）。本库已经有 src/sym/aes.c
 * 的 AES-128（纯 C 或 AES-NI 后端），这里用它的 AES128_ECB 实现同一语义：
 * 计数器块为 16 字节，前 4 字节是小端计数器，其余为零。AES-128 ECB 是标准
 * 运算，因此生成的公开矩阵与参考实现逐字节一致。
 *
 * 参考代码把 16*11 字节的缓冲区当作“密钥调度表”传递，这里只在前 16 字节
 * 保存原始密钥，其余字节清零。 */

#include "aes.h"
#include "sym/aes.h"

#include <stdint.h>
#include <string.h>

void AES128_load_schedule(const uint8_t *key, uint8_t *schedule)
{
    memcpy(schedule, key, 16U);
    memset(schedule + 16U, 0, 16U * 11U - 16U);
}

void AES128_free_schedule(uint8_t *schedule)
{
    memset(schedule, 0, 16U * 11U);
}

void AES128_CTR_zero_sch(uint32_t counter_start, size_t block_count,
                         const uint8_t *schedule, uint8_t *ciphertext)
{
    uint8_t blocks[block_count * 16U];

    for (size_t i = 0; i < block_count; i++) {
        uint8_t *block = blocks + 16U * i;
        uint32_t counter = counter_start + (uint32_t)i;

        memset(block, 0, 16U);
        block[0] = (uint8_t)counter;
        block[1] = (uint8_t)(counter >> 8);
        block[2] = (uint8_t)(counter >> 16);
        block[3] = (uint8_t)(counter >> 24);
    }
    AES128_ECB(schedule, blocks, ciphertext, block_count * 16U);
    memset(blocks, 0, sizeof(blocks));
}
