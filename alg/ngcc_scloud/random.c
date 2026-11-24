/* Scloud+ 参考实现的 randombytes 入口。
 *
 * 参考代码在 pke_keygen / KEM 封装内部调用 randombytes 取随机种子。这里把它
 * 接到 NGCC 的确定性 DRNG 上：常规调用由 adapter 在每次操作前用系统 RNG
 * 播种；KAT 测试直接对同一上下文播种，从而复现官方测试向量。 */

#include "random.h"
#include "ngcc/drng.h"

extern DRNG_ctx pqc_ngcc_scloud_drng;

int randombytes(unsigned char *buffer, unsigned int size)
{
    if (size == 0U)
        return 0;
    if (buffer == 0)
        return -1;
    return get_random_number(&pqc_ngcc_scloud_drng, buffer,
                             (unsigned long long)size * 8ULL);
}
