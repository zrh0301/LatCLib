#include "gf16.h"

/* 常数乘 shuffle 内核。
 * 输入：a 为普通 byte GF16 向量；k 为常数系数。
 * 输出：每个 lane 都乘以 k。
 * 场景：固定系数矩阵乘法，是 shuffle 家族的最基础形态。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_shuffle(pqc_kernel_gf16_vec a, uint8_t k) {
    return pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mtk_uniform[k & 0x0f], pqc_kernel_gf16_vec_mask_u8(a, 0x0f));
}

/* 双输出 shuffle 内核（mtk2）。
 * 输入：a 为普通 byte GF16 向量；low_k/high_k 是两个常数系数。
 * 输出：每个 lane 同时产生两个乘积，结果打进一个 packed byte。
 * 场景：同一个输入要同时供两路常数系数消费，例如 rank-pair 的双行更新。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_dual_coeff_shuffle(pqc_kernel_gf16_vec a, uint8_t low_k, uint8_t high_k) {
    return pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mtk2_16[(low_k & 0x0f) | ((high_k & 0x0f) << 4)], pqc_kernel_gf16_vec_mask_u8(a, 0x0f));
}

/* packed pair shuffle。
 * 输入：packed byte 的低/高 nibble 分别代表两个 GF16 元素；low_k/high_k 是两个常数系数。
 * 输出：低 nibble 乘 low_k，高 nibble 乘 high_k，再 pack 回一个 byte。
 * 场景：输入本身就是 packed 表示时，保留布局并避免中间展开。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_pair_shuffle(pqc_kernel_gf16_vec packed, uint8_t low_k, uint8_t high_k) {
    pqc_kernel_gf16_vec lo = pqc_kernel_gf16_vec_mask_u8(packed, 0x0f);
    pqc_kernel_gf16_vec hi = pqc_kernel_gf16_vec_srli4(packed);
    pqc_kernel_gf16_vec low_res = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mtk_uniform[low_k & 0x0f], lo);
    pqc_kernel_gf16_vec high_res = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mtk_uniform[high_k & 0x0f], hi);
    return pqc_kernel_gf16_vec_or(low_res, pqc_kernel_gf16_vec_slli4(high_res));
}

/* 向量化常数乘数组。
 * 输入：a 为 byte GF16 数组，k 为常数系数。
 * 输出：out[i] = a[i] * k。
 * 场景：低层 shuffle 基线，和标量常数乘对照。
 */
void pqc_kernel_gf16_mul_const_shuffle_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec v = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_const_shuffle(v, k);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    const uint8_t *tab = pqc_kernel_gf16_mul_table[k & 0x0f];
    for (; i < n; ++i) {
        out[i] = tab[a[i] & 0x0f];
    }
}

/* 向量化双常数系数乘数组。
 * 输入：a 为普通 byte GF16 数组；low_k/high_k 是两个常数系数。
 * 输出：一个 packed byte 同时存两路结果。
 * 场景：mtk2 双输出，适合两行或两累加器共享同一个输入。
 */
void pqc_kernel_gf16_mul_dual_coeff_shuffle_array(uint8_t *out, const uint8_t *a, uint8_t low_k, uint8_t high_k, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec v = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_dual_coeff_shuffle(v, low_k, high_k);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        uint8_t x = a[i] & 0x0f;
        out[i] = (uint8_t)(pqc_kernel_gf16_mul_table[x][low_k & 0x0f] |
                           (pqc_kernel_gf16_mul_table[x][high_k & 0x0f] << 4));
    }
}

/* packed pair shuffle 数组。
 * 输入：packed 的低/高 nibble 是两个 GF16 元素。
 * 输出：两个 nibble 分别乘以不同常数后再 pack。
 * 场景：packed L/R 之类已经打包好的数据通路。
 */
void pqc_kernel_gf16_mul_packed_pair_shuffle_array(uint8_t *out, const uint8_t *packed, uint8_t low_k, uint8_t high_k, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec v = pqc_kernel_gf16_vec_load(packed + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_packed_pair_shuffle(v, low_k, high_k);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        uint8_t lo = packed[i] & 0x0f;
        uint8_t hi = (packed[i] >> 4) & 0x0f;
        out[i] = (uint8_t)(pqc_kernel_gf16_mul_table[lo][low_k & 0x0f] |
                           (pqc_kernel_gf16_mul_table[hi][high_k & 0x0f] << 4));
    }
}

/* GF(16) 求逆的 SIMD 版本：16 项表正好落在一次 shuffle 的索引空间内。
 * 求逆路径比乘法路径更短：没有加法、没有 mod 15（log/exp 版除外）。 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_inv_table_shuffle(pqc_kernel_gf16_vec a) {
    return pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_inv_vec, pqc_kernel_gf16_vec_mask_u8(a, 0x0f));
}

/* 复用 3.2 节的 log/exp 表：inv(a) = exp[(15 - log a) mod 15]。 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_inv_logexp(pqc_kernel_gf16_vec a) {
    const pqc_kernel_gf16_vec idx = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);
    const pqc_kernel_gf16_vec lg = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, idx);
    pqc_kernel_gf16_vec t = pqc_kernel_gf16_vec_sub_u8(pqc_kernel_gf16_vec_set1(15u), lg);
    t = pqc_kernel_gf16_vec_and(t, pqc_kernel_gf16_vec_not(pqc_kernel_gf16_vec_eq(t, pqc_kernel_gf16_vec_set1(15u))));
    const pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, t);
    const pqc_kernel_gf16_vec nz = pqc_kernel_gf16_vec_not(pqc_kernel_gf16_vec_eq(idx, pqc_kernel_gf16_vec_zero()));
    return pqc_kernel_gf16_vec_and(r, nz);
}

void pqc_kernel_gf16_inv_table_shuffle_array(uint8_t *out, const uint8_t *a, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_vec_inv_table_shuffle(pqc_kernel_gf16_vec_load(a + i)));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf16_inv16_table[a[i] & 0x0fu];
}

void pqc_kernel_gf16_inv_logexp_vec_array(uint8_t *out, const uint8_t *a, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_vec_inv_logexp(pqc_kernel_gf16_vec_load(a + i)));
    }
    for (; i < n; ++i) out[i] = pqc_kernel_gf16_inv_logexp_scalar(a[i]);
}