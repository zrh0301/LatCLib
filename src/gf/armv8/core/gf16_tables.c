#include "gf16.h"

#include <string.h>

#if defined(_MSC_VER)
#define GF16_ALIGN(N) __declspec(align(N))
#else
#define GF16_ALIGN(N) __attribute__((aligned(N)))
#endif

uint8_t pqc_kernel_gf16_mul_table[16][16];
uint8_t pqc_kernel_gf16_log_table[16];
uint8_t pqc_kernel_gf16_exp_table[15];
pqc_kernel_gf16_vec pqc_kernel_gf16_mtk2_16[256];
pqc_kernel_gf16_vec pqc_kernel_gf16_mtk_uniform[16];
pqc_kernel_gf16_vec pqc_kernel_gf16_log_vec;
pqc_kernel_gf16_vec pqc_kernel_gf16_exp_vec;
pqc_kernel_gf16_vec pqc_kernel_gf16_mul2_vec;
pqc_kernel_gf16_vec pqc_kernel_gf16_mul4_vec;
pqc_kernel_gf16_vec pqc_kernel_gf16_mul8_vec;

/* GF16 参考乘法。
 * 输入：a、b 的低 4 bit 为 GF16 元素。
 * 输出：单个 GF16 乘积的低 4 bit。
 * 特点：按多项式基做移位异或，是最朴素但最可靠的正确性来源。
 */
uint8_t pqc_kernel_gf16_mul_ref(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    a &= 0x0f;
    b &= 0x0f;
    for (int i = 0; i < 4; ++i) {
        if (b & 1u) {
            p ^= a;
        }
        uint8_t carry = (uint8_t)(a & 0x08u);
        a <<= 1;
        if (carry) {
            a ^= 0x13u;
        }
        b >>= 1;
    }
    return (uint8_t)(p & 0x0f);
}

/* 把 16 项查表转换成一个 SIMD shuffle 常量向量。
 * 输入：tab16[16] 是一个 16 元素的 GF16 映射表。
 * 输出：一个向量常量，后续可直接喂给 shuffle 指令。
 * 场景：把“表”提前变成“向量常量”，避免热循环里重复装配。
 */
static void make_vec_table(pqc_kernel_gf16_vec *dst, const uint8_t tab16[16]) {
    GF16_ALIGN(32) uint8_t tmp[GF16_VEC_BYTES];
    for (size_t i = 0; i < GF16_VEC_BYTES; ++i) {
        tmp[i] = tab16[i & 0x0f];
    }
    *dst = pqc_kernel_gf16_vec_load(tmp);
}

/* 初始化所有 GF16 运行时常量。
 * 输入：无。
 * 输出：填好标量表、log/exp 表、双输出 mtk2 表、uniform shuffle 表、mul2/mul4/mul8 常量表。
 * 场景：程序启动时调用一次，后续所有算法族共享。
 */
void pqc_kernel_gf16_tables_init(void) {
    static int done = 0;
    if (done) {
        return;
    }

    memset(pqc_kernel_gf16_mul_table, 0, sizeof(pqc_kernel_gf16_mul_table));
    for (int a = 0; a < 16; ++a) {
        for (int b = 0; b < 16; ++b) {
            pqc_kernel_gf16_mul_table[a][b] = pqc_kernel_gf16_mul_ref((uint8_t)a, (uint8_t)b);
        }
    }

    static const uint8_t exp_init[15] = {
        1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9
    };
    memcpy(pqc_kernel_gf16_exp_table, exp_init, sizeof(pqc_kernel_gf16_exp_table));
    memset(pqc_kernel_gf16_log_table, 0, sizeof(pqc_kernel_gf16_log_table));
    for (int i = 0; i < 15; ++i) {
        pqc_kernel_gf16_log_table[pqc_kernel_gf16_exp_table[i]] = (uint8_t)i;
    }

    for (int low_k = 0; low_k < 16; ++low_k) {
        for (int high_k = 0; high_k < 16; ++high_k) {
            uint8_t tab16[16];
            for (int x = 0; x < 16; ++x) {
                uint8_t lo = pqc_kernel_gf16_mul_table[x & 0x0f][low_k];
                uint8_t hi = pqc_kernel_gf16_mul_table[x & 0x0f][high_k];
                tab16[x] = (uint8_t)(lo | (hi << 4));
            }
            make_vec_table(&pqc_kernel_gf16_mtk2_16[low_k | (high_k << 4)], tab16);
        }
    }

    for (int k = 0; k < 16; ++k) {
        uint8_t tab16[16];
        for (int x = 0; x < 16; ++x) {
            tab16[x] = pqc_kernel_gf16_mul_table[x][k];
        }
        make_vec_table(&pqc_kernel_gf16_mtk_uniform[k], tab16);
    }

    uint8_t log16[16];
    uint8_t exp16[16];
    uint8_t mul2[16];
    uint8_t mul4[16];
    uint8_t mul8[16];
    for (int i = 0; i < 16; ++i) {
        log16[i] = pqc_kernel_gf16_log_table[i];
        exp16[i] = (i < 15) ? pqc_kernel_gf16_exp_table[i] : 0;
        mul2[i] = pqc_kernel_gf16_mul_table[i][2];
        mul4[i] = pqc_kernel_gf16_mul_table[i][4];
        mul8[i] = pqc_kernel_gf16_mul_table[i][8];
    }
    make_vec_table(&pqc_kernel_gf16_log_vec, log16);
    make_vec_table(&pqc_kernel_gf16_exp_vec, exp16);
    make_vec_table(&pqc_kernel_gf16_mul2_vec, mul2);
    make_vec_table(&pqc_kernel_gf16_mul4_vec, mul4);
    make_vec_table(&pqc_kernel_gf16_mul8_vec, mul8);
    make_vec_table(&pqc_kernel_gf16_inv_vec, pqc_kernel_gf16_inv16_table);

    done = 1;
}

/* 16 项逆元表（正文 3.9 节）。约定 INV[0] = 0，使零元素自动落到零。 */
const uint8_t pqc_kernel_gf16_inv16_table[16] = {
    0x00, 0x01, 0x09, 0x0e, 0x0d, 0x0b, 0x07, 0x06,
    0x0f, 0x02, 0x0c, 0x05, 0x0a, 0x04, 0x03, 0x08
};
pqc_kernel_gf16_vec pqc_kernel_gf16_inv_vec;
