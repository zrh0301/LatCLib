#include "gf16.h"

/* 向量常数乘 log/exp 版本。
 * 输入：a 为普通 byte GF16 向量；k 为常数系数。
 * 输出：a*k。
 * 特点：先查 log，再做指数加法和模 15 归约，再查 exp。
 * 场景：研究 log/exp 在常数乘上的 SIMD 成本。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_modexp(pqc_kernel_gf16_vec a, uint8_t k) {
    a = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);
    if ((k & 0x0f) == 0) {
        return pqc_kernel_gf16_vec_zero();
    }
    pqc_kernel_gf16_vec zero = pqc_kernel_gf16_vec_zero();
    pqc_kernel_gf16_vec loga = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, a);
    pqc_kernel_gf16_vec e = pqc_kernel_gf16_vec_add_u8(loga, pqc_kernel_gf16_vec_set1(pqc_kernel_gf16_log_table[k & 0x0f]));
    pqc_kernel_gf16_vec ge15 = pqc_kernel_gf16_vec_gt_u8_small(e, pqc_kernel_gf16_vec_set1(14));
    e = pqc_kernel_gf16_vec_sub_u8(e, pqc_kernel_gf16_vec_and(ge15, pqc_kernel_gf16_vec_set1(15)));
    pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, e);
    pqc_kernel_gf16_vec nz = pqc_kernel_gf16_vec_eq(a, zero);
    return pqc_kernel_gf16_vec_and(r, pqc_kernel_gf16_vec_xor(nz, pqc_kernel_gf16_vec_set1(0xff)));
}

/* 向量变量乘变量 log/exp 版本。
 * 输入：a、b 为普通 byte GF16 向量。
 * 输出：逐 lane 乘积。
 * 特点：对称、统一、布局干净，但通常比 shuffle 更重。
 * 场景：变量乘变量、以及希望保留普通 byte 布局的通用场景。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_modexp(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) {
    a = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);
    b = pqc_kernel_gf16_vec_mask_u8(b, 0x0f);
    pqc_kernel_gf16_vec zero = pqc_kernel_gf16_vec_zero();
    pqc_kernel_gf16_vec loga = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, a);
    pqc_kernel_gf16_vec logb = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, b);
    pqc_kernel_gf16_vec e = pqc_kernel_gf16_vec_add_u8(loga, logb);
    pqc_kernel_gf16_vec ge15 = pqc_kernel_gf16_vec_gt_u8_small(e, pqc_kernel_gf16_vec_set1(14));
    e = pqc_kernel_gf16_vec_sub_u8(e, pqc_kernel_gf16_vec_and(ge15, pqc_kernel_gf16_vec_set1(15)));
    pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, e);
    pqc_kernel_gf16_vec az = pqc_kernel_gf16_vec_eq(a, zero);
    pqc_kernel_gf16_vec bz = pqc_kernel_gf16_vec_eq(b, zero);
    pqc_kernel_gf16_vec z = pqc_kernel_gf16_vec_or(az, bz);
    return pqc_kernel_gf16_vec_and(r, pqc_kernel_gf16_vec_xor(z, pqc_kernel_gf16_vec_set1(0xff)));
}

/* 向量常数乘数组 log/exp 版本。
 * 输入：a 为数组，k 为常数。
 * 输出：out[i] = a[i] * k。
 */
void pqc_kernel_gf16_mul_const_modexp_vec_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec v = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_const_modexp(v, k);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_modexp_scalar(a[i], k);
    }
}

/* 向量变量乘变量数组 log/exp 版本。
 * 输入：a/b 为普通 byte GF16 数组。
 * 输出：out[i] = a[i] * b[i]。
 */
void pqc_kernel_gf16_mul_modexp_vec_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec va = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec vb = pqc_kernel_gf16_vec_load(b + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_modexp(va, vb);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_modexp_scalar(a[i], b[i]);
    }
}
