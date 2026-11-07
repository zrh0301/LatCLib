// SPDX-License-Identifier: CC0 OR Apache-2.0
//
// pqc_kernel_gf16_gfni_pack.c -- GFNI 方法的 packed 变体（GFNI_pack_affine）。
//
// 与 pqc_kernel_gf16_gfni.c 的区别：
//   pqc_kernel_gf16_gfni.c 走"GF(16) 嵌入 GF(2^8) 子域 + GF2P8MULB"，输入输出是
//   一字节一个元素（unpacked）；GF2P8MULB 在 packed 布局下不可用，因为
//   GF(256) 乘法的进位会跨 nibble 串扰。
//
//   本文件走"块对角 8x8 GF(2) 线性映射 + GF2P8AFFINEQB"，输入输出全程
//   packed（一字节两个 GF(16) 元素），一条 affine 指令处理 32 字节 = 64 个元素。
//
// 数学：k 固定时 v -> v*k 是 GF(2) 线性的，且低 nibble 只依赖低 nibble、
// 高 nibble 只依赖高 nibble，所以整字节映射是两块 4x4 的 M_k 做块对角：
//       A_k = [[M_k, 0], [0, M_k]]
// GFNI 编码约定（Intel SDM, GF2P8AFFINEQB）：
//       dst.bit[i] = parity(src AND A.byte[7-i]) XOR b_i
// 即 A.byte[m] 存"输出 bit (7-m)"那一行。
//
// 参考常量表与构造算法见：
//   GF_cal/tools/pqc_kernel_gf16_pack_affine_matrix.py   （生成 + 逐项校验）

#include "gf16.h"

#if defined(GF16_ENABLE_GFNI)

#include <immintrin.h>
#include <stdint.h>
#include <string.h>

/*
 * 16 个 A_k，k = 0..15。表体积 128 字节，L1 常驻。
 * 由 pqc_kernel_gf16_pack_affine_matrix.py 生成；对全部 16 x 256 组 (k, v) 与直接
 * GF(16) 乘法比对通过。
 * 其中 k = 1,2,4,8 四项与pqov src/gfni/pqc_kernel_gf16_gfni.h 的 __gf16_mulbit
 * 的 4 个 qword 逐字节一致，可作为交叉验证锚点。
 */
static const uint64_t pqc_kernel_gf16_pack_affine_table[16] = {
    UINT64_C(0x0000000000000000), /* k =  0 */
    UINT64_C(0x0102040810204080), /* k =  1 */
    UINT64_C(0x0809020480902040), /* k =  2 */
    UINT64_C(0x090b060c90b060c0), /* k =  3 */
    UINT64_C(0x040c090240c09020), /* k =  4 */
    UINT64_C(0x050e0d0a50e0d0a0), /* k =  5 */
    UINT64_C(0x0c050b06c050b060), /* k =  6 */
    UINT64_C(0x0d070f0ed070f0e0), /* k =  7 */
    UINT64_C(0x02060c092060c090), /* k =  8 */
    UINT64_C(0x0304080130408010), /* k =  9 */
    UINT64_C(0x0a0f0e0da0f0e0d0), /* k = 10 */
    UINT64_C(0x0b0d0a05b0d0a050), /* k = 11 */
    UINT64_C(0x060a050b60a050b0), /* k = 12 */
    UINT64_C(0x0708010370801030), /* k = 13 */
    UINT64_C(0x0e03070fe03070f0), /* k = 14 */
    UINT64_C(0x0f010307f0103070), /* k = 15 */
};

/* 4 个基矩阵 A_1, A_2, A_4, A_8（"M 对 k 线性"现拼用） */
static const uint64_t pqc_kernel_gf16_pack_affine_basis[4] = {
    UINT64_C(0x0102040810204080), /* A_1 */
    UINT64_C(0x0809020480902040), /* A_2 */
    UINT64_C(0x040c090240c09020), /* A_4 */
    UINT64_C(0x02060c092060c090), /* A_8 */
};

/* 算法 1：由 k 现构造 A_k。
 * 用于校验静态表，以及量化"现构造 vs 查表"的成本差。 */
static uint64_t pqc_kernel_gf16_pack_affine_build(uint8_t k) {
    uint8_t img[4];
    uint8_t row[8];
    for (int j = 0; j < 4; ++j) {
        img[j] = pqc_kernel_gf16_mul_table[(uint8_t)(1u << j)][k & 0x0f];
    }
    for (int i = 0; i < 4; ++i) {
        uint8_t m = 0;
        for (int j = 0; j < 4; ++j) {
            if ((img[j] >> i) & 1u) m |= (uint8_t)(1u << j);
        }
        row[i] = m;
        row[i + 4] = (uint8_t)(m << 4);
    }
    uint64_t a = 0;
    for (int m = 0; m < 8; ++m) {
        a |= (uint64_t)row[7 - m] << (8 * m);
    }
    return a;
}

/* 基矩阵线性组合构造：只存 4 个基矩阵，用 k 的位把它们 XOR 起来。
 * 因为 M_k 对 k 是 GF(2)-线性的，A_k = XOR_{bit t of k} A_{2^t}。 */
static uint64_t pqc_kernel_gf16_pack_affine_build_from_basis(uint8_t k) {
    uint64_t a = 0;
    for (int t = 0; t < 4; ++t) {
        if ((k >> t) & 1u) a ^= pqc_kernel_gf16_pack_affine_basis[t];
    }
    return a;
}

/*
 * 自检：静态表 == 算法 1 构造 == 基矩阵现拼，且 affine 结果 == 直接 GF(16) 乘法。
 * 返回 0 表示全部通过，否则返回第一个失败处的编码。
 */
static int pqc_kernel_gf16_gfni_pack_selftest_scalar(uint64_t x, uint64_t a) {
    /* 按 SDM 语义在标量上模拟 GF2P8AFFINEQB（b = 0） */
    uint8_t bytes[8];
    uint8_t out = 0;
    for (int m = 0; m < 8; ++m) bytes[m] = (uint8_t)((a >> (8 * m)) & 0xff);
    for (int i = 0; i < 8; ++i) {
        uint8_t t = (uint8_t)(x & bytes[7 - i]);
        t ^= (uint8_t)(t >> 4);
        t ^= (uint8_t)(t >> 2);
        t ^= (uint8_t)(t >> 1);
        out |= (uint8_t)((t & 1u) << i);
    }
    return out;
}

int pqc_kernel_gf16_gfni_pack_selftest(void) {
    for (int k = 0; k < 16; ++k) {
        const uint64_t a0 = pqc_kernel_gf16_pack_affine_table[k];
        if (a0 != pqc_kernel_gf16_pack_affine_build((uint8_t)k)) return 1;
        if (a0 != pqc_kernel_gf16_pack_affine_build_from_basis((uint8_t)k)) return 2;
        for (int v = 0; v < 256; ++v) {
            /* packed 语义：低 nibble 出低 nibble，高 nibble 出高 nibble */
            const uint8_t want = (uint8_t)(
                pqc_kernel_gf16_mul_table[v & 0x0f][k] |
                (pqc_kernel_gf16_mul_table[(v >> 4) & 0x0f][k] << 4));
            const uint8_t got = pqc_kernel_gf16_gfni_pack_selftest_scalar(
                (uint64_t)(unsigned)v, a0);
            if (got != want) return 3;
        }
    }
    return 0;
}

/* 算法 2：packed 向量 x 常数。1 条 affine / 32 字节 / 64 元素。 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_gfni_pack(pqc_kernel_gf16_vec a, uint8_t k) {
    const __m256i mat = _mm256_set1_epi64x((long long)pqc_kernel_gf16_pack_affine_table[k & 0x0f]);
    return _mm256_gf2p8affine_epi64_epi8(a, mat, 0);
}

/* 算法 3：packed 乘加 acc ^= a x k */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_madd_gfni_pack(pqc_kernel_gf16_vec acc, pqc_kernel_gf16_vec a, uint8_t k) {
    return _mm256_xor_si256(acc, pqc_kernel_gf16_vec_mul_const_gfni_pack(a, k));
}

/* 数组形式，含标量尾部 */
void pqc_kernel_gf16_mul_const_gfni_pack_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    const uint8_t kk = (uint8_t)(k & 0x0f);
    const __m256i mat = _mm256_set1_epi64x((long long)pqc_kernel_gf16_pack_affine_table[kk]);
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec_store(out + i,
                       _mm256_gf2p8affine_epi64_epi8(pqc_kernel_gf16_vec_load(a + i), mat, 0));
    }
    for (; i < n; ++i) {
        const uint8_t v = a[i];
        out[i] = (uint8_t)(pqc_kernel_gf16_mul_table[v & 0x0f][kk] |
                           (pqc_kernel_gf16_mul_table[(v >> 4) & 0x0f][kk] << 4));
    }
}

/* 数组形式，乘加：acc ^= a x k（packed） */
void pqc_kernel_gf16_mulacc_const_gfni_pack_array(uint8_t *acc, const uint8_t *a,
                                       uint8_t k, size_t n) {
    const uint8_t kk = (uint8_t)(k & 0x0f);
    const __m256i mat = _mm256_set1_epi64x((long long)pqc_kernel_gf16_pack_affine_table[kk]);
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf16_vec p =
            _mm256_gf2p8affine_epi64_epi8(pqc_kernel_gf16_vec_load(a + i), mat, 0);
        pqc_kernel_gf16_vec_store(acc + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(acc + i), p));
    }
    for (; i < n; ++i) {
        const uint8_t v = a[i];
        acc[i] ^= (uint8_t)(pqc_kernel_gf16_mul_table[v & 0x0f][kk] |
                            (pqc_kernel_gf16_mul_table[(v >> 4) & 0x0f][kk] << 4));
    }
}

/* 算法 4：packed 常数收缩 P*S。
 * coeff[r*cols+c] 是低 nibble 标量；sig 与 out 全程 packed。 */
void pqc_kernel_gf16_snova_ps_gfni_pack(uint8_t *out, const uint8_t *coeff,
                             const uint8_t *sig, size_t rows, size_t cols,
                             size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_madd_gfni_pack(
                    acc, pqc_kernel_gf16_vec_load(sig + c * lane_count + i),
                    (uint8_t)(coeff[r * cols + c] & 0x0f));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t v = sig[c * lane_count + i];
                const uint8_t kk = (uint8_t)(coeff[r * cols + c] & 0x0f);
                acc ^= (uint8_t)(pqc_kernel_gf16_mul_table[v & 0x0f][kk] |
                                 (pqc_kernel_gf16_mul_table[(v >> 4) & 0x0f][kk] << 4));
            }
            dst[i] = acc;
        }
    }
}

/* 供基准量化"表构造方式"差异的三个入口 */
uint64_t pqc_kernel_gf16_gfni_pack_matrix_static(uint8_t k) {
    return pqc_kernel_gf16_pack_affine_table[k & 0x0f];
}

uint64_t pqc_kernel_gf16_gfni_pack_matrix_build(uint8_t k) {
    return pqc_kernel_gf16_pack_affine_build((uint8_t)(k & 0x0f));
}

uint64_t pqc_kernel_gf16_gfni_pack_matrix_from_basis(uint8_t k) {
    return pqc_kernel_gf16_pack_affine_build_from_basis((uint8_t)(k & 0x0f));
}

const uint64_t *pqc_kernel_gf16_gfni_pack_matrix_table(void) {
    return pqc_kernel_gf16_pack_affine_table;
}

#else /* !GF16_ENABLE_GFNI */

#include <string.h>

int pqc_kernel_gf16_gfni_pack_selftest(void) { return -1; }

pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_gfni_pack(pqc_kernel_gf16_vec a, uint8_t k) {
    (void)k;
    return a;
}

pqc_kernel_gf16_vec pqc_kernel_gf16_vec_madd_gfni_pack(pqc_kernel_gf16_vec acc, pqc_kernel_gf16_vec a, uint8_t k) {
    (void)k;
    return pqc_kernel_gf16_vec_xor(acc, a);
}

void pqc_kernel_gf16_mul_const_gfni_pack_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    (void)k;
    memcpy(out, a, n);
}

void pqc_kernel_gf16_mulacc_const_gfni_pack_array(uint8_t *acc, const uint8_t *a, uint8_t k, size_t n) {
    (void)k;
    for (size_t i = 0; i < n; ++i) acc[i] ^= a[i];
}

void pqc_kernel_gf16_snova_ps_gfni_pack(uint8_t *out, const uint8_t *coeff,
                             const uint8_t *sig, size_t rows, size_t cols,
                             size_t lane_count) {
    (void)coeff;
    (void)sig;
    (void)cols;
    memset(out, 0, rows * lane_count);
}

uint64_t pqc_kernel_gf16_gfni_pack_matrix_static(uint8_t k) { (void)k; return 0; }
uint64_t pqc_kernel_gf16_gfni_pack_matrix_build(uint8_t k) { (void)k; return 0; }
uint64_t pqc_kernel_gf16_gfni_pack_matrix_from_basis(uint8_t k) { (void)k; return 0; }
const uint64_t *pqc_kernel_gf16_gfni_pack_matrix_table(void) { return 0; }

#endif /* GF16_ENABLE_GFNI */
