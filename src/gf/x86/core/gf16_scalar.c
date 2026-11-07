#include "gf16.h"

/* 标量表乘。
 * 输入：a、b 的低 4 bit 为 GF16 元素。
 * 输出：a*b 的低 4 bit。
 * 特点：直接查 16x16 表，简单稳妥。
 */
uint8_t pqc_kernel_gf16_mul_table_scalar(uint8_t a, uint8_t b) {
    return pqc_kernel_gf16_mul_table[a & 0x0f][b & 0x0f];
}

/* 标量 log/exp 乘。
 * 输入：a、b 的低 4 bit 为 GF16 元素。
 * 输出：a*b 的低 4 bit。
 * 特点：把非零元素映射成指数相加，再做模 15 归约；适合研究 log/exp 代价。
 */
uint8_t pqc_kernel_gf16_mul_modexp_scalar(uint8_t a, uint8_t b) {
    a &= 0x0f;
    b &= 0x0f;
    if ((a | b) == 0 || a == 0 || b == 0) {
        return 0;
    }
    uint8_t e = (uint8_t)(pqc_kernel_gf16_log_table[a] + pqc_kernel_gf16_log_table[b]);
    if (e >= 15) {
        e = (uint8_t)(e - 15);
    }
    return pqc_kernel_gf16_exp_table[e];
}

/* 标量数组表乘。
 * 输入：a/b 为普通 byte GF16 数组；n 为元素数。
 * 输出：out[i] = a[i] * b[i]。
 * 场景：scalar baseline、正确性参考。
 */
void pqc_kernel_gf16_mul_table_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_table[a[i] & 0x0f][b[i] & 0x0f];
    }
}

/* 标量数组 log/exp 乘。
 * 输入：a/b 为普通 byte GF16 数组。
 * 输出：逐元素乘积。
 * 场景：和表乘做标量口径比较。
 */
void pqc_kernel_gf16_mul_modexp_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_modexp_scalar(a[i], b[i]);
    }
}

/* 标量常数乘表。
 * 输入：a 是数组，k 是常数系数。
 * 输出：out[i] = a[i] * k。
 * 场景：固定系数矩阵作用、后续可对照 SIMD shuffle。
 */
void pqc_kernel_gf16_mul_const_table_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    const uint8_t *tab = pqc_kernel_gf16_mul_table[k & 0x0f];
    for (size_t i = 0; i < n; ++i) {
        out[i] = tab[a[i] & 0x0f];
    }
}

/* 标量常数 log/exp 乘。
 * 输入：a 是数组，k 是常数系数。
 * 输出：out[i] = a[i] * k。
 * 场景：研究 modexp 在常数乘上的成本。
 */
void pqc_kernel_gf16_mul_const_modexp_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_modexp_scalar(a[i], k);
    }
}

/* GF(16) 求逆（正文 3.9 节）：16 项逆元表、复用 log/exp 表、费马幂 a^14。 */
uint8_t pqc_kernel_gf16_inv_table_scalar(uint8_t a) {
    return pqc_kernel_gf16_inv16_table[a & 0x0fu];
}

uint8_t pqc_kernel_gf16_inv_logexp_scalar(uint8_t a) {
    const uint8_t v = (uint8_t)(a & 0x0fu);
    if (!v) return 0;
    return pqc_kernel_gf16_exp_table[(15u - pqc_kernel_gf16_log_table[v]) % 15u];
}

uint8_t pqc_kernel_gf16_inv_fermat_scalar(uint8_t a) {
    const uint8_t v = (uint8_t)(a & 0x0fu);
    const uint8_t a2 = pqc_kernel_gf16_mul_ref(v, v);
    const uint8_t a4 = pqc_kernel_gf16_mul_ref(a2, a2);
    const uint8_t a8 = pqc_kernel_gf16_mul_ref(a4, a4);
    return pqc_kernel_gf16_mul_ref(pqc_kernel_gf16_mul_ref(a2, a4), a8);
}

void pqc_kernel_gf16_inv_fermat_array(uint8_t *out, const uint8_t *a, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = pqc_kernel_gf16_inv_fermat_scalar(a[i]);
}