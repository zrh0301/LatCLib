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

/* 4 路 ILP 版 log/exp 乘法：四个独立链交错，隐藏 vpshufb 与 ALU 延迟。 */
void pqc_kernel_gf16_vec_mul_modexp_ilp4(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                              pqc_kernel_gf16_vec b0, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b3,
                              pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1, pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3) {
    const pqc_kernel_gf16_vec zero = pqc_kernel_gf16_vec_zero();
    const pqc_kernel_gf16_vec v14  = pqc_kernel_gf16_vec_set1(14);
    const pqc_kernel_gf16_vec v15  = pqc_kernel_gf16_vec_set1(15);
    const pqc_kernel_gf16_vec all  = pqc_kernel_gf16_vec_set1(0xff);

    a0 = pqc_kernel_gf16_vec_mask_u8(a0, 0x0f); a1 = pqc_kernel_gf16_vec_mask_u8(a1, 0x0f);
    a2 = pqc_kernel_gf16_vec_mask_u8(a2, 0x0f); a3 = pqc_kernel_gf16_vec_mask_u8(a3, 0x0f);
    b0 = pqc_kernel_gf16_vec_mask_u8(b0, 0x0f); b1 = pqc_kernel_gf16_vec_mask_u8(b1, 0x0f);
    b2 = pqc_kernel_gf16_vec_mask_u8(b2, 0x0f); b3 = pqc_kernel_gf16_vec_mask_u8(b3, 0x0f);

    /* 阶段 1：四路 log 查表（八次 shuffle 互不依赖） */
    const pqc_kernel_gf16_vec la0 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, a0);
    const pqc_kernel_gf16_vec la1 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, a1);
    const pqc_kernel_gf16_vec la2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, a2);
    const pqc_kernel_gf16_vec la3 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, a3);
    const pqc_kernel_gf16_vec lb0 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, b0);
    const pqc_kernel_gf16_vec lb1 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, b1);
    const pqc_kernel_gf16_vec lb2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, b2);
    const pqc_kernel_gf16_vec lb3 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_log_vec, b3);

    /* 阶段 2：四路指数相加 + mod 15（四路互不依赖） */
    pqc_kernel_gf16_vec e0 = pqc_kernel_gf16_vec_add_u8(la0, lb0);
    pqc_kernel_gf16_vec e1 = pqc_kernel_gf16_vec_add_u8(la1, lb1);
    pqc_kernel_gf16_vec e2 = pqc_kernel_gf16_vec_add_u8(la2, lb2);
    pqc_kernel_gf16_vec e3 = pqc_kernel_gf16_vec_add_u8(la3, lb3);
    e0 = pqc_kernel_gf16_vec_sub_u8(e0, pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec_gt_u8_small(e0, v14), v15));
    e1 = pqc_kernel_gf16_vec_sub_u8(e1, pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec_gt_u8_small(e1, v14), v15));
    e2 = pqc_kernel_gf16_vec_sub_u8(e2, pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec_gt_u8_small(e2, v14), v15));
    e3 = pqc_kernel_gf16_vec_sub_u8(e3, pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec_gt_u8_small(e3, v14), v15));

    /* 阶段 3：四路 exp 查表 */
    const pqc_kernel_gf16_vec x0 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, e0);
    const pqc_kernel_gf16_vec x1 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, e1);
    const pqc_kernel_gf16_vec x2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, e2);
    const pqc_kernel_gf16_vec x3 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_exp_vec, e3);

    /* 阶段 4：零掩码（独立于上述所有阶段，可并行） */
    const pqc_kernel_gf16_vec z0 = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_eq(a0, zero), pqc_kernel_gf16_vec_eq(b0, zero));
    const pqc_kernel_gf16_vec z1 = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_eq(a1, zero), pqc_kernel_gf16_vec_eq(b1, zero));
    const pqc_kernel_gf16_vec z2 = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_eq(a2, zero), pqc_kernel_gf16_vec_eq(b2, zero));
    const pqc_kernel_gf16_vec z3 = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_eq(a3, zero), pqc_kernel_gf16_vec_eq(b3, zero));

    *r0 = pqc_kernel_gf16_vec_and(x0, pqc_kernel_gf16_vec_xor(z0, all));
    *r1 = pqc_kernel_gf16_vec_and(x1, pqc_kernel_gf16_vec_xor(z1, all));
    *r2 = pqc_kernel_gf16_vec_and(x2, pqc_kernel_gf16_vec_xor(z2, all));
    *r3 = pqc_kernel_gf16_vec_and(x3, pqc_kernel_gf16_vec_xor(z3, all));
}

/* 4 路 ILP 数组版：每次取 4 个向量同时算。 */
void pqc_kernel_gf16_mul_modexp_vec_array_ilp(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + 4u * GF16_VEC_BYTES <= n; i += 4u * GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec a0 = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec a1 = pqc_kernel_gf16_vec_load(a + i + GF16_VEC_BYTES);
        pqc_kernel_gf16_vec a2 = pqc_kernel_gf16_vec_load(a + i + 2u * GF16_VEC_BYTES);
        pqc_kernel_gf16_vec a3 = pqc_kernel_gf16_vec_load(a + i + 3u * GF16_VEC_BYTES);
        pqc_kernel_gf16_vec b0 = pqc_kernel_gf16_vec_load(b + i);
        pqc_kernel_gf16_vec b1 = pqc_kernel_gf16_vec_load(b + i + GF16_VEC_BYTES);
        pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_load(b + i + 2u * GF16_VEC_BYTES);
        pqc_kernel_gf16_vec b3 = pqc_kernel_gf16_vec_load(b + i + 3u * GF16_VEC_BYTES);
        pqc_kernel_gf16_vec r0, r1, r2, r3;
        pqc_kernel_gf16_vec_mul_modexp_ilp4(a0, a1, a2, a3, b0, b1, b2, b3, &r0, &r1, &r2, &r3);
        pqc_kernel_gf16_vec_store(out + i, r0);
        pqc_kernel_gf16_vec_store(out + i + GF16_VEC_BYTES, r1);
        pqc_kernel_gf16_vec_store(out + i + 2u * GF16_VEC_BYTES, r2);
        pqc_kernel_gf16_vec_store(out + i + 3u * GF16_VEC_BYTES, r3);
    }
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_modexp(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i));
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_modexp_scalar(a[i], b[i]);
    }
}
