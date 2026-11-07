#include "gf16.h"

/* packed byte 的常数乘辅助。
 * 输入：packed 的低/高 nibble 各表示一个 GF16 元素；k 是常数系数。
 * 输出：低 nibble 与高 nibble 分别乘以 k，再 pack 回一个 byte。
 * 场景：packed pair / packed-a 的局部步骤。
 */
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_k_native(pqc_kernel_gf16_vec packed, uint8_t k) {
    pqc_kernel_gf16_vec lo = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mtk_uniform[k & 0x0f], pqc_kernel_gf16_vec_mask_u8(packed, 0x0f));
    pqc_kernel_gf16_vec hi = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mtk_uniform[k & 0x0f], pqc_kernel_gf16_vec_srli4(packed));
    return pqc_kernel_gf16_vec_or(lo, pqc_kernel_gf16_vec_slli4(hi));
}

/* bit-plane 乘法的“右侧预计算”辅助。
 * 输入：a 是左操作数，b1/b2/b4/b8 是右操作数的四个 bit-plane 选择结果。
 * 输出：a*b。
 * 特点：先看 a 的 bit，再 XOR 选择对应的右侧平面。
 */
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_precomputed_b_native(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8) {
    a = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);
    pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_zero();
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b1, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(1)), pqc_kernel_gf16_vec_set1(1))));
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b2, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(2)), pqc_kernel_gf16_vec_set1(2))));
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b4, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(4)), pqc_kernel_gf16_vec_set1(4))));
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b8, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(8)), pqc_kernel_gf16_vec_set1(8))));
    return r;
}

static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8) {
    const pqc_kernel_gf16_vec v1 = pqc_kernel_gf16_vec_set1(1);
    const pqc_kernel_gf16_vec v2 = pqc_kernel_gf16_vec_set1(2);
    const pqc_kernel_gf16_vec v4 = pqc_kernel_gf16_vec_set1(4);
    const pqc_kernel_gf16_vec v8 = pqc_kernel_gf16_vec_set1(8);

    a = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);

    const pqc_kernel_gf16_vec m1 = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, v1), v1);
    const pqc_kernel_gf16_vec m2 = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, v2), v2);
    const pqc_kernel_gf16_vec m4 = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, v4), v4);
    const pqc_kernel_gf16_vec m8 = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, v8), v8);

    const pqc_kernel_gf16_vec t1 = pqc_kernel_gf16_vec_and(b1, m1);
    const pqc_kernel_gf16_vec t2 = pqc_kernel_gf16_vec_and(b2, m2);
    const pqc_kernel_gf16_vec t4 = pqc_kernel_gf16_vec_and(b4, m4);
    const pqc_kernel_gf16_vec t8 = pqc_kernel_gf16_vec_and(b8, m8);

    const pqc_kernel_gf16_vec x12 = pqc_kernel_gf16_vec_xor(t1, t2);
    const pqc_kernel_gf16_vec x48 = pqc_kernel_gf16_vec_xor(t4, t8);
    return pqc_kernel_gf16_vec_xor(x12, x48);
}

/* 普通 bit-plane 变量乘变量。
 * 输入：a、b 为普通 byte GF16 向量。
 * 输出：逐 lane 乘积。
 * 场景：变量乘变量，尤其是右操作数可被多次复用时。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) {
    a = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);
    b = pqc_kernel_gf16_vec_mask_u8(b, 0x0f);
    pqc_kernel_gf16_vec b1 = b;
    pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b);
    pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b);
    pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b);
    pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_zero();
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b1, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(1)), pqc_kernel_gf16_vec_set1(1))));
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b2, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(2)), pqc_kernel_gf16_vec_set1(2))));
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b4, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(4)), pqc_kernel_gf16_vec_set1(4))));
    r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b8, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(8)), pqc_kernel_gf16_vec_set1(8))));
    return r;
}

/* ILP 版普通 bit-plane 变量乘变量。
 * 输入/输出与上一个函数相同。
 * 特点：先集中生成四个 bit-plane，再用树形 XOR 归约。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_bitplane_ilp(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) {
    a = pqc_kernel_gf16_vec_mask_u8(a, 0x0f);
    b = pqc_kernel_gf16_vec_mask_u8(b, 0x0f);
    pqc_kernel_gf16_vec b1 = b;
    pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b);
    pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b);
    pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b);
    return pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a, b1, b2, b4, b8);
}

void pqc_kernel_gf16_vec_mul_same_b4_ilp(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                              pqc_kernel_gf16_vec b, pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1,
                              pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3) {
    b = pqc_kernel_gf16_vec_mask_u8(b, 0x0f);
    const pqc_kernel_gf16_vec b1 = b;
    const pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b);
    const pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b);
    const pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b);
    *r0 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(pqc_kernel_gf16_vec_mask_u8(a0, 0x0f), b1, b2, b4, b8);
    *r1 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(pqc_kernel_gf16_vec_mask_u8(a1, 0x0f), b1, b2, b4, b8);
    *r2 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(pqc_kernel_gf16_vec_mask_u8(a2, 0x0f), b1, b2, b4, b8);
    *r3 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(pqc_kernel_gf16_vec_mask_u8(a3, 0x0f), b1, b2, b4, b8);
}

/* bit-plane 数组版本。
 * 输入：a/b 为普通 byte GF16 数组。
 * 输出：out[i] = a[i] * b[i]。
 * 场景：和 table / modexp 做标量口径对照。
 */
void pqc_kernel_gf16_mul_bitplane_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec va = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec vb = pqc_kernel_gf16_vec_load(b + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_bitplane(va, vb);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_table[a[i] & 0x0f][b[i] & 0x0f];
    }
}

/* ILP 版 bit-plane 数组版本。
 * 输入/输出与上一个函数相同。
 * 场景：对照普通 bit-plane，检查依赖链展开是否有收益。
 */
void pqc_kernel_gf16_mul_bitplane_array_ilp(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec va = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec vb = pqc_kernel_gf16_vec_load(b + i);
        pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_bitplane_ilp(va, vb);
        pqc_kernel_gf16_vec_store(out + i, r);
    }
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_mul_table[a[i] & 0x0f][b[i] & 0x0f];
    }
}

/* 同一个 b 更新 4 个累加器。
 * 输入：a0..a3 是四路左操作数；b 是同一个右操作数。
 * 输出：s0..s3 分别累加 a_j * b。
 * 场景：SNOVA 里同一右操作数被多行/多块复用的情况。
 */
void pqc_kernel_gf16_mul_bitplane_same_b4_array(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                     const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                     const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec vb = pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(b + i), 0x0f);
        pqc_kernel_gf16_vec b1 = vb;
        pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, vb);
        pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, vb);
        pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, vb);
        pqc_kernel_gf16_vec av[4] = {
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a0 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a1 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a2 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a3 + i), 0x0f)
        };
        uint8_t *sum[4] = {s0, s1, s2, s3};
        for (int j = 0; j < 4; ++j) {
            pqc_kernel_gf16_vec a = av[j];
            pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_zero();
            r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b1, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(1)), pqc_kernel_gf16_vec_set1(1))));
            r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b2, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(2)), pqc_kernel_gf16_vec_set1(2))));
            r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b4, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(4)), pqc_kernel_gf16_vec_set1(4))));
            r = pqc_kernel_gf16_vec_xor(r, pqc_kernel_gf16_vec_and(b8, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a, pqc_kernel_gf16_vec_set1(8)), pqc_kernel_gf16_vec_set1(8))));
            pqc_kernel_gf16_vec old = pqc_kernel_gf16_vec_load(sum[j] + i);
            pqc_kernel_gf16_vec_store(sum[j] + i, pqc_kernel_gf16_vec_xor(old, r));
        }
    }
    for (; i < n; ++i) {
        s0[i] ^= pqc_kernel_gf16_mul_table[a0[i] & 0x0f][b[i] & 0x0f];
        s1[i] ^= pqc_kernel_gf16_mul_table[a1[i] & 0x0f][b[i] & 0x0f];
        s2[i] ^= pqc_kernel_gf16_mul_table[a2[i] & 0x0f][b[i] & 0x0f];
        s3[i] ^= pqc_kernel_gf16_mul_table[a3[i] & 0x0f][b[i] & 0x0f];
    }
}

/* ILP 版同一个 b 更新 4 个累加器。
 * 输入/输出与上一个函数相同。
 * 特点：复用 ILP 版 bit-plane 内核，降低单链条 XOR 依赖。
 */
void pqc_kernel_gf16_mul_bitplane_same_b4_array_ilp(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                         const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                         const uint8_t *b, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec vb = pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(b + i), 0x0f);
        pqc_kernel_gf16_vec b1 = vb;
        pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, vb);
        pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, vb);
        pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, vb);
        pqc_kernel_gf16_vec av[4] = {
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a0 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a1 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a2 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a3 + i), 0x0f)
        };
        uint8_t *sum[4] = {s0, s1, s2, s3};
        for (int j = 0; j < 4; ++j) {
            pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(av[j], b1, b2, b4, b8);
            pqc_kernel_gf16_vec old = pqc_kernel_gf16_vec_load(sum[j] + i);
            pqc_kernel_gf16_vec_store(sum[j] + i, pqc_kernel_gf16_vec_xor(old, r));
        }
    }
    for (; i < n; ++i) {
        s0[i] ^= pqc_kernel_gf16_mul_table[a0[i] & 0x0f][b[i] & 0x0f];
        s1[i] ^= pqc_kernel_gf16_mul_table[a1[i] & 0x0f][b[i] & 0x0f];
        s2[i] ^= pqc_kernel_gf16_mul_table[a2[i] & 0x0f][b[i] & 0x0f];
        s3[i] ^= pqc_kernel_gf16_mul_table[a3[i] & 0x0f][b[i] & 0x0f];
    }
}

/* 同一个 packed b 更新 4 个累加器。
 * 输入：a0..a3 是四路左操作数；b_packed 的低/高 nibble 分别代表两路右操作数。
 * 输出：每个累加器都得到对应的 packed 结果。
 * 场景：packed right operand 的复用路径。
 */
void pqc_kernel_gf16_mul_bitplane_same_packed_b4_array(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                            const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                            const uint8_t *b_packed, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec b1 = pqc_kernel_gf16_vec_load(b_packed + i);
        pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_mul_packed_k_native(b1, 2);
        pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_mul_packed_k_native(b1, 4);
        pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_mul_packed_k_native(b1, 8);
        pqc_kernel_gf16_vec r0 = pqc_kernel_gf16_vec_mul_precomputed_b_native(pqc_kernel_gf16_vec_load(a0 + i), b1, b2, b4, b8);
        pqc_kernel_gf16_vec r1 = pqc_kernel_gf16_vec_mul_precomputed_b_native(pqc_kernel_gf16_vec_load(a1 + i), b1, b2, b4, b8);
        pqc_kernel_gf16_vec r2 = pqc_kernel_gf16_vec_mul_precomputed_b_native(pqc_kernel_gf16_vec_load(a2 + i), b1, b2, b4, b8);
        pqc_kernel_gf16_vec r3 = pqc_kernel_gf16_vec_mul_precomputed_b_native(pqc_kernel_gf16_vec_load(a3 + i), b1, b2, b4, b8);
        pqc_kernel_gf16_vec_store(s0 + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(s0 + i), r0));
        pqc_kernel_gf16_vec_store(s1 + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(s1 + i), r1));
        pqc_kernel_gf16_vec_store(s2 + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(s2 + i), r2));
        pqc_kernel_gf16_vec_store(s3 + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(s3 + i), r3));
    }
    for (; i < n; ++i) {
        uint8_t lo = b_packed[i] & 0x0f;
        uint8_t hi = (b_packed[i] >> 4) & 0x0f;
        s0[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a0[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a0[i] & 0x0f][hi] << 4));
        s1[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a1[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a1[i] & 0x0f][hi] << 4));
        s2[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a2[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a2[i] & 0x0f][hi] << 4));
        s3[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a3[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a3[i] & 0x0f][hi] << 4));
    }
}

void pqc_kernel_gf16_mul_bitplane_same_packed_b4_ilp_array(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                                const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                                const uint8_t *b_packed, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec vb = pqc_kernel_gf16_vec_load(b_packed + i);
        pqc_kernel_gf16_vec b1 = vb;
        pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_mul_packed_k_native(vb, 2);
        pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_mul_packed_k_native(vb, 4);
        pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_mul_packed_k_native(vb, 8);
        pqc_kernel_gf16_vec av[4] = {
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a0 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a1 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a2 + i), 0x0f),
            pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec_load(a3 + i), 0x0f)
        };
        uint8_t *sum[4] = {s0, s1, s2, s3};
        for (int j = 0; j < 4; ++j) {
            pqc_kernel_gf16_vec r = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(av[j], b1, b2, b4, b8);
            pqc_kernel_gf16_vec old = pqc_kernel_gf16_vec_load(sum[j] + i);
            pqc_kernel_gf16_vec_store(sum[j] + i, pqc_kernel_gf16_vec_xor(old, r));
        }
    }
    for (; i < n; ++i) {
        uint8_t lo = b_packed[i] & 0x0f;
        uint8_t hi = (b_packed[i] >> 4) & 0x0f;
        s0[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a0[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a0[i] & 0x0f][hi] << 4));
        s1[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a1[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a1[i] & 0x0f][hi] << 4));
        s2[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a2[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a2[i] & 0x0f][hi] << 4));
        s3[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[a3[i] & 0x0f][lo] | (pqc_kernel_gf16_mul_table[a3[i] & 0x0f][hi] << 4));
    }
}

/* packed-a 的 split 双输出。
 * 输入：a_packed 的低/高 nibble 是两个不同 GF16 元素；b1/b2/b4/b8 是右侧预计算结果。
 * 输出：out_lo = a_lo*b，out_hi = a_hi*b。
 * 场景：packed 左输入，且希望长期保持 packed 输出的一条路径。
 */
void pqc_kernel_gf16_vec_mul_packed_a_split_dual(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8,
                                      pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi) {
    pqc_kernel_gf16_vec mask = pqc_kernel_gf16_vec_set1(0x0f);
    pqc_kernel_gf16_vec a_lo = pqc_kernel_gf16_vec_and(a_packed, mask);
    pqc_kernel_gf16_vec a_hi = pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec_srli4(a_packed), mask);

    pqc_kernel_gf16_vec r_lo = pqc_kernel_gf16_vec_zero();
    r_lo = pqc_kernel_gf16_vec_xor(r_lo, pqc_kernel_gf16_vec_and(b1, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_lo, pqc_kernel_gf16_vec_set1(1)), pqc_kernel_gf16_vec_set1(1))));
    r_lo = pqc_kernel_gf16_vec_xor(r_lo, pqc_kernel_gf16_vec_and(b2, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_lo, pqc_kernel_gf16_vec_set1(2)), pqc_kernel_gf16_vec_set1(2))));
    r_lo = pqc_kernel_gf16_vec_xor(r_lo, pqc_kernel_gf16_vec_and(b4, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_lo, pqc_kernel_gf16_vec_set1(4)), pqc_kernel_gf16_vec_set1(4))));
    r_lo = pqc_kernel_gf16_vec_xor(r_lo, pqc_kernel_gf16_vec_and(b8, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_lo, pqc_kernel_gf16_vec_set1(8)), pqc_kernel_gf16_vec_set1(8))));

    pqc_kernel_gf16_vec r_hi = pqc_kernel_gf16_vec_zero();
    r_hi = pqc_kernel_gf16_vec_xor(r_hi, pqc_kernel_gf16_vec_and(b1, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_hi, pqc_kernel_gf16_vec_set1(1)), pqc_kernel_gf16_vec_set1(1))));
    r_hi = pqc_kernel_gf16_vec_xor(r_hi, pqc_kernel_gf16_vec_and(b2, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_hi, pqc_kernel_gf16_vec_set1(2)), pqc_kernel_gf16_vec_set1(2))));
    r_hi = pqc_kernel_gf16_vec_xor(r_hi, pqc_kernel_gf16_vec_and(b4, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_hi, pqc_kernel_gf16_vec_set1(4)), pqc_kernel_gf16_vec_set1(4))));
    r_hi = pqc_kernel_gf16_vec_xor(r_hi, pqc_kernel_gf16_vec_and(b8, pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_hi, pqc_kernel_gf16_vec_set1(8)), pqc_kernel_gf16_vec_set1(8))));

    *out_lo = r_lo;
    *out_hi = r_hi;
}

/* packed-a split 的 ILP 双输出。
 * 输入/输出与上一个函数相同。
 * 特点：把 a 的高低 nibble 拆分后，直接复用 ILP 版 bit-plane 核心。
 */
static inline void pqc_kernel_gf16_vec_mul_packed_a_split_dual_ilp(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8,
                                                        pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi) {
    pqc_kernel_gf16_vec mask = pqc_kernel_gf16_vec_set1(0x0f);
    pqc_kernel_gf16_vec a_lo = pqc_kernel_gf16_vec_and(a_packed, mask);
    pqc_kernel_gf16_vec a_hi = pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec_srli4(a_packed), mask);
    *out_lo = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a_lo, b1, b2, b4, b8);
    *out_hi = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a_hi, b1, b2, b4, b8);
}

/* packed-a 的 compare-mask 双输出。
 * 输入/输出与 split 相同。
 * 特点：直接对 packed byte 的高/低比特位做比较生成掩码。
 * 场景：作为 split 的对照实现，用来比较代价。
 */
void pqc_kernel_gf16_vec_mul_packed_a_cmp_dual(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8,
                                    pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi) {
    pqc_kernel_gf16_vec m1_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x01)), pqc_kernel_gf16_vec_set1(0x01));
    pqc_kernel_gf16_vec m2_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x02)), pqc_kernel_gf16_vec_set1(0x02));
    pqc_kernel_gf16_vec m4_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x04)), pqc_kernel_gf16_vec_set1(0x04));
    pqc_kernel_gf16_vec m8_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x08)), pqc_kernel_gf16_vec_set1(0x08));

    pqc_kernel_gf16_vec m1_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x10)), pqc_kernel_gf16_vec_set1(0x10));
    pqc_kernel_gf16_vec m2_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x20)), pqc_kernel_gf16_vec_set1(0x20));
    pqc_kernel_gf16_vec m4_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x40)), pqc_kernel_gf16_vec_set1(0x40));
    pqc_kernel_gf16_vec m8_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, pqc_kernel_gf16_vec_set1(0x80)), pqc_kernel_gf16_vec_set1(0x80));

    pqc_kernel_gf16_vec lo = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_and(b1, m1_lo), pqc_kernel_gf16_vec_and(b2, m2_lo)),
                                            pqc_kernel_gf16_vec_and(b4, m4_lo)),
                               pqc_kernel_gf16_vec_and(b8, m8_lo));
    pqc_kernel_gf16_vec hi = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_and(b1, m1_hi), pqc_kernel_gf16_vec_and(b2, m2_hi)),
                                            pqc_kernel_gf16_vec_and(b4, m4_hi)),
                               pqc_kernel_gf16_vec_and(b8, m8_hi));
    *out_lo = lo;
    *out_hi = hi;
}

/* packed-a compare-mask 的 ILP 双输出。
 * 输入/输出与上一个函数相同。
 * 特点：先集中生成掩码，再树形归约 low/high 两路结果。
 */
static inline void pqc_kernel_gf16_vec_mul_packed_a_cmp_dual_ilp(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8,
                                                     pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi) {
    const pqc_kernel_gf16_vec v1 = pqc_kernel_gf16_vec_set1(0x01);
    const pqc_kernel_gf16_vec v2 = pqc_kernel_gf16_vec_set1(0x02);
    const pqc_kernel_gf16_vec v4 = pqc_kernel_gf16_vec_set1(0x04);
    const pqc_kernel_gf16_vec v8 = pqc_kernel_gf16_vec_set1(0x08);
    const pqc_kernel_gf16_vec v1h = pqc_kernel_gf16_vec_set1(0x10);
    const pqc_kernel_gf16_vec v2h = pqc_kernel_gf16_vec_set1(0x20);
    const pqc_kernel_gf16_vec v4h = pqc_kernel_gf16_vec_set1(0x40);
    const pqc_kernel_gf16_vec v8h = pqc_kernel_gf16_vec_set1(0x80);

    pqc_kernel_gf16_vec m1_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v1), v1);
    pqc_kernel_gf16_vec m2_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v2), v2);
    pqc_kernel_gf16_vec m4_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v4), v4);
    pqc_kernel_gf16_vec m8_lo = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v8), v8);

    pqc_kernel_gf16_vec m1_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v1h), v1h);
    pqc_kernel_gf16_vec m2_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v2h), v2h);
    pqc_kernel_gf16_vec m4_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v4h), v4h);
    pqc_kernel_gf16_vec m8_hi = pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec_and(a_packed, v8h), v8h);

    pqc_kernel_gf16_vec lo1 = pqc_kernel_gf16_vec_and(b1, m1_lo);
    pqc_kernel_gf16_vec lo2 = pqc_kernel_gf16_vec_and(b2, m2_lo);
    pqc_kernel_gf16_vec lo3 = pqc_kernel_gf16_vec_and(b4, m4_lo);
    pqc_kernel_gf16_vec lo4 = pqc_kernel_gf16_vec_and(b8, m8_lo);
    pqc_kernel_gf16_vec hi1 = pqc_kernel_gf16_vec_and(b1, m1_hi);
    pqc_kernel_gf16_vec hi2 = pqc_kernel_gf16_vec_and(b2, m2_hi);
    pqc_kernel_gf16_vec hi3 = pqc_kernel_gf16_vec_and(b4, m4_hi);
    pqc_kernel_gf16_vec hi4 = pqc_kernel_gf16_vec_and(b8, m8_hi);

    *out_lo = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_xor(lo1, lo2), pqc_kernel_gf16_vec_xor(lo3, lo4));
    *out_hi = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_xor(hi1, hi2), pqc_kernel_gf16_vec_xor(hi3, hi4));
}

/* packed-a/packed-b 变量乘法。
 * 输入两个 packed byte，输出 low*low 与 high*high 两路普通向量。
 * 输出不重新 pack，供 verify 收缩循环在寄存器/普通累加器中继续使用。
 */
void pqc_kernel_gf16_vec_mul_packed_ab_split_dual(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed,
                                       pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi) {
    const pqc_kernel_gf16_vec a_lo = pqc_kernel_gf16_vec_mask_u8(a_packed, 0x0f);
    const pqc_kernel_gf16_vec a_hi = pqc_kernel_gf16_vec_srli4(a_packed);
    const pqc_kernel_gf16_vec b_lo = pqc_kernel_gf16_vec_mask_u8(b_packed, 0x0f);
    const pqc_kernel_gf16_vec b_hi = pqc_kernel_gf16_vec_srli4(b_packed);
    const pqc_kernel_gf16_vec b_lo_2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b_lo);
    const pqc_kernel_gf16_vec b_lo_4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b_lo);
    const pqc_kernel_gf16_vec b_lo_8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b_lo);
    const pqc_kernel_gf16_vec b_hi_2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b_hi);
    const pqc_kernel_gf16_vec b_hi_4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b_hi);
    const pqc_kernel_gf16_vec b_hi_8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b_hi);
    *out_lo = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a_lo, b_lo, b_lo_2, b_lo_4, b_lo_8);
    *out_hi = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a_hi, b_hi, b_hi_2, b_hi_4, b_hi_8);
}

/* 全 packed 变量乘（目标形态）：A_p × B_p，一次 pass 产出两条 packed 向量。
 * out_straight = [a_hi*b_hi | a_lo*b_lo]
 * out_cross    = [a_hi*b_lo | a_lo*b_hi]
 * 全程不拆包、不重打包；4 个叉积一次 pass 全部得到。
 *
 * 实现手法：先把 B_p 一次性生成 packed 四平面 b1,b2,b4,b8（两个 nibble 同时乘
 * 2/4/8），再用 cmp 风格掩码按 A_p 的 8 个 bit 选择，得到两个“半成品” packed：
 *   lo_all = [a_lo*b_hi | a_lo*b_lo]
 *   hi_all = [a_hi*b_hi | a_hi*b_lo]
 * 最后只做 nibble 合并（AND+OR / 移位），不涉及任何拆包重算。
 */
void pqc_kernel_gf16_vec_mul_packed_ab_dual_packed(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed,
                                        pqc_kernel_gf16_vec *out_straight, pqc_kernel_gf16_vec *out_cross) {
    const pqc_kernel_gf16_vec b1 = b_packed;
    const pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 2);
    const pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 4);
    const pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 8);
    pqc_kernel_gf16_vec lo_all, hi_all;
    pqc_kernel_gf16_vec_mul_packed_a_cmp_dual(a_packed, b1, b2, b4, b8, &lo_all, &hi_all);
    *out_straight = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_mask_u8(lo_all, 0x0f), pqc_kernel_gf16_vec_mask_u8(hi_all, 0xf0));
    *out_cross    = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_slli4(pqc_kernel_gf16_vec_mask_u8(hi_all, 0x0f)), pqc_kernel_gf16_vec_srli4(lo_all));
}

/* 全 packed 变量乘，只输出直配 packed 结果（verify S^T*T 的常用形态）。 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_ab_straight(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed) {
    const pqc_kernel_gf16_vec b1 = b_packed;
    const pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 2);
    const pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 4);
    const pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 8);
    pqc_kernel_gf16_vec lo_all, hi_all;
    pqc_kernel_gf16_vec_mul_packed_a_cmp_dual(a_packed, b1, b2, b4, b8, &lo_all, &hi_all);
    return pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_mask_u8(lo_all, 0x0f), pqc_kernel_gf16_vec_mask_u8(hi_all, 0xf0));
}

/* 全 packed 变量乘（ILP 版，只输出直配 packed 结果）。 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_ab_straight_ilp(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed) {
    const pqc_kernel_gf16_vec b1 = b_packed;
    const pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 2);
    const pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 4);
    const pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 8);
    pqc_kernel_gf16_vec lo_all, hi_all;
    pqc_kernel_gf16_vec_mul_packed_a_cmp_dual_ilp(a_packed, b1, b2, b4, b8, &lo_all, &hi_all);
    return pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_mask_u8(lo_all, 0x0f), pqc_kernel_gf16_vec_mask_u8(hi_all, 0xf0));
}

/* 全 packed 变量乘（split 版，只输出直配 packed 结果）：
 * 拆成低/高两路各自 bit-plane 乘，再拼回一个 packed byte。 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_ab_straight_split(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed) {
    pqc_kernel_gf16_vec lo, hi;
    pqc_kernel_gf16_vec_mul_packed_ab_split_dual(a_packed, b_packed, &lo, &hi);
    return pqc_kernel_gf16_vec_or(lo, pqc_kernel_gf16_vec_slli4(hi));
}

/* 全 packed 变量乘（ILP 版，直配 + 交叉）。 */
void pqc_kernel_gf16_vec_mul_packed_ab_dual_packed_ilp(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed,
                                            pqc_kernel_gf16_vec *out_straight, pqc_kernel_gf16_vec *out_cross) {
    const pqc_kernel_gf16_vec b1 = b_packed;
    const pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 2);
    const pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 4);
    const pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_mul_packed_k_native(b_packed, 8);
    pqc_kernel_gf16_vec lo_all, hi_all;
    pqc_kernel_gf16_vec_mul_packed_a_cmp_dual_ilp(a_packed, b1, b2, b4, b8, &lo_all, &hi_all);
    *out_straight = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_mask_u8(lo_all, 0x0f), pqc_kernel_gf16_vec_mask_u8(hi_all, 0xf0));
    *out_cross    = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_slli4(pqc_kernel_gf16_vec_mask_u8(hi_all, 0x0f)), pqc_kernel_gf16_vec_srli4(lo_all));
}

/* 全 packed 变量乘数组（直配）。输入/输出都是 packed，长度 n 是 packed byte 数。 */
void pqc_kernel_gf16_mul_packed_ab_straight_array(uint8_t *out, const uint8_t *a_packed,
                                       const uint8_t *b_packed, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_vec_mul_packed_ab_straight(
            pqc_kernel_gf16_vec_load(a_packed + i), pqc_kernel_gf16_vec_load(b_packed + i)));
    }
    for (; i < n; ++i) {
        const uint8_t a_lo = a_packed[i] & 0x0f;
        const uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        const uint8_t b_lo = b_packed[i] & 0x0f;
        const uint8_t b_hi = (b_packed[i] >> 4) & 0x0f;
        out[i] = (uint8_t)(pqc_kernel_gf16_mul_table[a_lo][b_lo] | (pqc_kernel_gf16_mul_table[a_hi][b_hi] << 4));
    }
}

/* 全 packed 变量乘数组（ILP，直配）。 */
void pqc_kernel_gf16_mul_packed_ab_straight_array_ilp(uint8_t *out, const uint8_t *a_packed,
                                           const uint8_t *b_packed, size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i, pqc_kernel_gf16_vec_mul_packed_ab_straight_ilp(
            pqc_kernel_gf16_vec_load(a_packed + i), pqc_kernel_gf16_vec_load(b_packed + i)));
    }
    for (; i < n; ++i) {
        const uint8_t a_lo = a_packed[i] & 0x0f;
        const uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        const uint8_t b_lo = b_packed[i] & 0x0f;
        const uint8_t b_hi = (b_packed[i] >> 4) & 0x0f;
        out[i] = (uint8_t)(pqc_kernel_gf16_mul_table[a_lo][b_lo] | (pqc_kernel_gf16_mul_table[a_hi][b_hi] << 4));
    }
}

/* 全 packed 变量乘数组（直配 + 交叉，4 个乘积一次 pass）。 */
void pqc_kernel_gf16_mul_packed_ab_dual_packed_array(uint8_t *out_straight, uint8_t *out_cross,
                                          const uint8_t *a_packed, const uint8_t *b_packed,
                                          size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec s, c;
        pqc_kernel_gf16_vec_mul_packed_ab_dual_packed(pqc_kernel_gf16_vec_load(a_packed + i),
                                           pqc_kernel_gf16_vec_load(b_packed + i), &s, &c);
        pqc_kernel_gf16_vec_store(out_straight + i, s);
        pqc_kernel_gf16_vec_store(out_cross + i, c);
    }
    for (; i < n; ++i) {
        const uint8_t a_lo = a_packed[i] & 0x0f;
        const uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        const uint8_t b_lo = b_packed[i] & 0x0f;
        const uint8_t b_hi = (b_packed[i] >> 4) & 0x0f;
        out_straight[i] = (uint8_t)(pqc_kernel_gf16_mul_table[a_lo][b_lo] | (pqc_kernel_gf16_mul_table[a_hi][b_hi] << 4));
        out_cross[i]    = (uint8_t)(pqc_kernel_gf16_mul_table[a_lo][b_hi] | (pqc_kernel_gf16_mul_table[a_hi][b_lo] << 4));
    }
}

/* 双流位平面 × 共享右操作数（4 行复用，与 same-b4 同构）：
 * right packed 的低/高 nibble 各生成一套四平面（只算一次），
 * 四行 packed 左输入各自按 bit 选择，得到四个 packed 直配结果。 */
void pqc_kernel_gf16_vec_mul_packed_ab_same_b4(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                                    pqc_kernel_gf16_vec b_packed,
                                    pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1, pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3) {
    const pqc_kernel_gf16_vec b_lo = pqc_kernel_gf16_vec_mask_u8(b_packed, 0x0f);
    const pqc_kernel_gf16_vec b_hi = pqc_kernel_gf16_vec_srli4(b_packed);
    const pqc_kernel_gf16_vec b_lo_2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b_lo);
    const pqc_kernel_gf16_vec b_lo_4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b_lo);
    const pqc_kernel_gf16_vec b_lo_8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b_lo);
    const pqc_kernel_gf16_vec b_hi_2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b_hi);
    const pqc_kernel_gf16_vec b_hi_4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b_hi);
    const pqc_kernel_gf16_vec b_hi_8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b_hi);

    const pqc_kernel_gf16_vec av[4] = {a0, a1, a2, a3};
    pqc_kernel_gf16_vec *rv[4] = {r0, r1, r2, r3};
    for (int j = 0; j < 4; ++j) {
        const pqc_kernel_gf16_vec a_lo = pqc_kernel_gf16_vec_mask_u8(av[j], 0x0f);
        const pqc_kernel_gf16_vec a_hi = pqc_kernel_gf16_vec_srli4(av[j]);
        const pqc_kernel_gf16_vec lo = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a_lo, b_lo, b_lo_2, b_lo_4, b_lo_8);
        const pqc_kernel_gf16_vec hi = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(a_hi, b_hi, b_hi_2, b_hi_4, b_hi_8);
        *rv[j] = pqc_kernel_gf16_vec_or(lo, pqc_kernel_gf16_vec_slli4(hi));
    }
}

/* ILP 版 packed-A/packed-B x 共享 packed-B（4 路输出）。 */
void pqc_kernel_gf16_vec_mul_packed_ab_same_b4_ilp(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                                        pqc_kernel_gf16_vec b_packed,
                                        pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1,
                                        pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3) {
    const pqc_kernel_gf16_vec b_lo = pqc_kernel_gf16_vec_mask_u8(b_packed, 0x0f);
    const pqc_kernel_gf16_vec b_hi = pqc_kernel_gf16_vec_srli4(b_packed);
    const pqc_kernel_gf16_vec b_lo_2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b_lo);
    const pqc_kernel_gf16_vec b_lo_4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b_lo);
    const pqc_kernel_gf16_vec b_lo_8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b_lo);
    const pqc_kernel_gf16_vec b_hi_2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b_hi);
    const pqc_kernel_gf16_vec b_hi_4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b_hi);
    const pqc_kernel_gf16_vec b_hi_8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b_hi);
    const pqc_kernel_gf16_vec a0_lo = pqc_kernel_gf16_vec_mask_u8(a0, 0x0f);
    const pqc_kernel_gf16_vec a1_lo = pqc_kernel_gf16_vec_mask_u8(a1, 0x0f);
    const pqc_kernel_gf16_vec a2_lo = pqc_kernel_gf16_vec_mask_u8(a2, 0x0f);
    const pqc_kernel_gf16_vec a3_lo = pqc_kernel_gf16_vec_mask_u8(a3, 0x0f);
    const pqc_kernel_gf16_vec a0_hi = pqc_kernel_gf16_vec_srli4(a0);
    const pqc_kernel_gf16_vec a1_hi = pqc_kernel_gf16_vec_srli4(a1);
    const pqc_kernel_gf16_vec a2_hi = pqc_kernel_gf16_vec_srli4(a2);
    const pqc_kernel_gf16_vec a3_hi = pqc_kernel_gf16_vec_srli4(a3);
    const pqc_kernel_gf16_vec lo0 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a0_lo, b_lo, b_lo_2, b_lo_4, b_lo_8);
    const pqc_kernel_gf16_vec lo1 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a1_lo, b_lo, b_lo_2, b_lo_4, b_lo_8);
    const pqc_kernel_gf16_vec lo2 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a2_lo, b_lo, b_lo_2, b_lo_4, b_lo_8);
    const pqc_kernel_gf16_vec lo3 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a3_lo, b_lo, b_lo_2, b_lo_4, b_lo_8);
    const pqc_kernel_gf16_vec hi0 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a0_hi, b_hi, b_hi_2, b_hi_4, b_hi_8);
    const pqc_kernel_gf16_vec hi1 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a1_hi, b_hi, b_hi_2, b_hi_4, b_hi_8);
    const pqc_kernel_gf16_vec hi2 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a2_hi, b_hi, b_hi_2, b_hi_4, b_hi_8);
    const pqc_kernel_gf16_vec hi3 = pqc_kernel_gf16_vec_mul_precomputed_b_native_ilp(
        a3_hi, b_hi, b_hi_2, b_hi_4, b_hi_8);
    *r0 = pqc_kernel_gf16_vec_or(lo0, pqc_kernel_gf16_vec_slli4(hi0));
    *r1 = pqc_kernel_gf16_vec_or(lo1, pqc_kernel_gf16_vec_slli4(hi1));
    *r2 = pqc_kernel_gf16_vec_or(lo2, pqc_kernel_gf16_vec_slli4(hi2));
    *r3 = pqc_kernel_gf16_vec_or(lo3, pqc_kernel_gf16_vec_slli4(hi3));
}

/* packed-a split 双输出数组版本。
 * 输入：a_packed 为 packed 左输入；b1/b2/b4/b8 是右侧预计算表。
 * 输出：out_lo/out_hi 两路结果。
 */
void pqc_kernel_gf16_mul_packed_a_split_dual_array(uint8_t *out_lo, uint8_t *out_hi,
                                        const uint8_t *a_packed,
                                        const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                        size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec lo, hi;
        pqc_kernel_gf16_vec_mul_packed_a_split_dual(pqc_kernel_gf16_vec_load(a_packed + i),
                                         pqc_kernel_gf16_vec_load(b1 + i),
                                         pqc_kernel_gf16_vec_load(b2 + i),
                                         pqc_kernel_gf16_vec_load(b4 + i),
                                         pqc_kernel_gf16_vec_load(b8 + i),
                                         &lo, &hi);
        pqc_kernel_gf16_vec_store(out_lo + i, lo);
        pqc_kernel_gf16_vec_store(out_hi + i, hi);
    }
    for (; i < n; ++i) {
        uint8_t a_lo = a_packed[i] & 0x0f;
        uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        uint8_t base = b1[i] & 0x0f;
        out_lo[i] = pqc_kernel_gf16_mul_table[a_lo][base];
        out_hi[i] = pqc_kernel_gf16_mul_table[a_hi][base];
    }
}

void pqc_kernel_gf16_mul_packed_a_split_dual_array_ilp(uint8_t *out_lo, uint8_t *out_hi,
                                            const uint8_t *a_packed,
                                            const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                            size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec lo, hi;
        pqc_kernel_gf16_vec_mul_packed_a_split_dual_ilp(pqc_kernel_gf16_vec_load(a_packed + i),
                                             pqc_kernel_gf16_vec_load(b1 + i),
                                             pqc_kernel_gf16_vec_load(b2 + i),
                                             pqc_kernel_gf16_vec_load(b4 + i),
                                             pqc_kernel_gf16_vec_load(b8 + i),
                                             &lo, &hi);
        pqc_kernel_gf16_vec_store(out_lo + i, lo);
        pqc_kernel_gf16_vec_store(out_hi + i, hi);
    }
    for (; i < n; ++i) {
        uint8_t a_lo = a_packed[i] & 0x0f;
        uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        uint8_t base = b1[i] & 0x0f;
        out_lo[i] = pqc_kernel_gf16_mul_table[a_lo][base];
        out_hi[i] = pqc_kernel_gf16_mul_table[a_hi][base];
    }
}

/* packed-a cmp 双输出数组版本。
 * 输入/输出与 split 相同。
 * 场景：packed-a split 的对照实现。
 */
void pqc_kernel_gf16_mul_packed_a_cmp_dual_array(uint8_t *out_lo, uint8_t *out_hi,
                                      const uint8_t *a_packed,
                                      const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                      size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec lo, hi;
        pqc_kernel_gf16_vec_mul_packed_a_cmp_dual(pqc_kernel_gf16_vec_load(a_packed + i),
                                       pqc_kernel_gf16_vec_load(b1 + i),
                                       pqc_kernel_gf16_vec_load(b2 + i),
                                       pqc_kernel_gf16_vec_load(b4 + i),
                                       pqc_kernel_gf16_vec_load(b8 + i),
                                       &lo, &hi);
        pqc_kernel_gf16_vec_store(out_lo + i, lo);
        pqc_kernel_gf16_vec_store(out_hi + i, hi);
    }
    for (; i < n; ++i) {
        uint8_t a_lo = a_packed[i] & 0x0f;
        uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        uint8_t base = b1[i] & 0x0f;
        out_lo[i] = pqc_kernel_gf16_mul_table[a_lo][base];
        out_hi[i] = pqc_kernel_gf16_mul_table[a_hi][base];
    }
}

void pqc_kernel_gf16_mul_packed_a_cmp_dual_array_ilp(uint8_t *out_lo, uint8_t *out_hi,
                                          const uint8_t *a_packed,
                                          const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                          size_t n) {
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec lo, hi;
        pqc_kernel_gf16_vec_mul_packed_a_cmp_dual_ilp(pqc_kernel_gf16_vec_load(a_packed + i),
                                           pqc_kernel_gf16_vec_load(b1 + i),
                                           pqc_kernel_gf16_vec_load(b2 + i),
                                           pqc_kernel_gf16_vec_load(b4 + i),
                                           pqc_kernel_gf16_vec_load(b8 + i),
                                           &lo, &hi);
        pqc_kernel_gf16_vec_store(out_lo + i, lo);
        pqc_kernel_gf16_vec_store(out_hi + i, hi);
    }
    for (; i < n; ++i) {
        uint8_t a_lo = a_packed[i] & 0x0f;
        uint8_t a_hi = (a_packed[i] >> 4) & 0x0f;
        uint8_t base = b1[i] & 0x0f;
        out_lo[i] = pqc_kernel_gf16_mul_table[a_lo][base];
        out_hi[i] = pqc_kernel_gf16_mul_table[a_hi][base];
    }
}
