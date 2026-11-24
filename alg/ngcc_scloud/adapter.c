/* Scloud+（NGCC 第一轮 KEM 候选，格基非结构 LWE）到本库 KEM API 的适配层。
 *
 * 参考实现：ICCS NGCC 第一轮提交包 Reference_Implementation/Scloudplus-128，
 * 主密钥流为 AES 族（矩阵 A 用 AES-128-CTR 展开，F/G/H/K 用带标签的
 * SHAKE256/SHA3-512）。这里只接入提交默认的 128 位安全级、AES 族。
 *
 * 本地改动：
 *   - 参考实现的 randombytes 改为由 NGCC DRNG 提供（见 random.c），KAT 测试
 *     直接对 pqc_ngcc_scloud_drng 播种即可逐字节复现官方向量；
 *   - AES-128-CTR 行展开改由库内 src/sym/aes.c 统一实现（见 aes_shim.c）；
 *   - 所有对外符号统一加 pqc_ngcc_scloud_ 前缀（namespace.h）。
 *
 * 评审报告 kem-35-1 只影响 256/384/512 等级（BD6/BD12 拒绝采样泄漏解密
 * 消息类别）；128 级使用 BD4/BD2 非拒绝采样，属于报告中的对照实现，因此
 * 默认的 128 级不触发该问题。 */

#include "kem.h"
#include "core/pqc_kem.h"
#include "util/rand.h"
#include "ngcc/drng.h"

/* 参考实现的 randombytes（pke_keygen / kem 封装）与 KAT 测试共用此上下文。 */
DRNG_ctx pqc_ngcc_scloud_drng;

static int seed_drng(void)
{
    uint8_t seed[32];

    if (pqc_randombytes(seed, sizeof(seed), NULL, NULL) != 0)
        return -1;
    if (init_random_number(&pqc_ngcc_scloud_drng, seed, sizeof(seed)) != 0)
        return -1;
    return 0;
}

static int keypair(uint8_t *pk, uint8_t *sk)
{
    if (seed_drng() != 0)
        return -1;
    return scloud_kemkeygen(pk, sk);
}

static int encaps(uint8_t *ct, uint8_t *ss, const uint8_t *pk)
{
    if (seed_drng() != 0)
        return -1;
    return scloud_kemencaps(pk, ct, ss);
}

static int decaps(uint8_t *ss, const uint8_t *ct, const uint8_t *sk)
{
    return scloud_kemdecaps(sk, ct, ss);
}

static const pqc_kem algorithm = {
    SYSTEM_NAME,
    SCLOUDPLUS_PUBLIC_KEY_BYTES,
    SCLOUDPLUS_SECRET_KEY_BYTES,
    SCLOUDPLUS_CIPHERTEXT_BYTES,
    SCLOUDPLUS_SHARED_SECRET_BYTES,
    keypair,
    encaps,
    decaps,
};

const pqc_kem *pqc_ngcc_scloud_algorithm(void)
{
    return &algorithm;
}
