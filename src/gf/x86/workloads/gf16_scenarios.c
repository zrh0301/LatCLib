#include "gf16.h"

#include <stdlib.h>
#include <string.h>

/* SNOVA-like P*S。
 * 输入：coeff 是 rows x cols 的常数系数；sig 是 cols 组输入向量。
 * 输出：out 是 rows 组收缩结果。
 * 场景：常数矩阵乘普通向量，是固定系数路径的核心模型。
 */
void pqc_kernel_gf16_snova_ps_table(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                         size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            const uint8_t *src = sig + c * lane_count;
            for (size_t i = 0; i < lane_count; ++i) {
                dst[i] ^= pqc_kernel_gf16_mul_table[src[i] & 0x0f][k];
            }
        }
    }
}

/* P*S 的标量 log/exp 版本。 */
void pqc_kernel_gf16_snova_ps_modexp_scalar(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                                 size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            const uint8_t *src = sig + c * lane_count;
            for (size_t i = 0; i < lane_count; ++i) {
                dst[i] ^= pqc_kernel_gf16_mul_modexp_scalar(src[i], k);
            }
        }
    }
}

/* P*S 的 SIMD shuffle 版本。
 * 输入/输出与 table 版相同。
 * 特点：固定系数可直接做 shuffle，通常是常数乘场景的主力。
 */
void pqc_kernel_gf16_snova_ps_shuffle(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                           size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t k = coeff[r * cols + c] & 0x0f;
                acc = pqc_kernel_gf16_vec_xor(acc,
                    pqc_kernel_gf16_vec_mul_const_shuffle(pqc_kernel_gf16_vec_load(sig + c * lane_count + i), k));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[sig[c * lane_count + i] & 0x0f]
                                      [coeff[r * cols + c] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

/* P*S 的 SIMD log/exp 版本。
 * 输入/输出与 table 版相同。
 * 特点：保留普通 byte 布局，便于和变量乘变量统一。
 */
void pqc_kernel_gf16_snova_ps_modexp_vec(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                              size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t k = coeff[r * cols + c] & 0x0f;
                acc = pqc_kernel_gf16_vec_xor(acc,
                    pqc_kernel_gf16_vec_mul_const_modexp(pqc_kernel_gf16_vec_load(sig + c * lane_count + i), k));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_modexp_scalar(sig[c * lane_count + i],
                                              coeff[r * cols + c] & 0x0f);
            }
            dst[i] = acc;
        }
    }
}

/* P*S 的 bit-plane 版本。
 * 输入/输出与 table 版相同。
 * 特点：按 b,2b,4b,8b 选项做乘加，适合和多累加器复用结合。
 */
void pqc_kernel_gf16_snova_ps_bitplane(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                            size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const pqc_kernel_gf16_vec kvec = pqc_kernel_gf16_vec_set1(coeff[r * cols + c] & 0x0f);
                acc = pqc_kernel_gf16_vec_xor(acc,
                    pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec_load(sig + c * lane_count + i), kvec));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[sig[c * lane_count + i] & 0x0f]
                                      [coeff[r * cols + c] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

void pqc_kernel_gf16_snova_ps_bitplane_ilp(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                                size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const pqc_kernel_gf16_vec kvec = pqc_kernel_gf16_vec_set1(coeff[r * cols + c] & 0x0f);
                acc = pqc_kernel_gf16_vec_xor(acc,
                    pqc_kernel_gf16_vec_mul_bitplane_ilp(pqc_kernel_gf16_vec_load(sig + c * lane_count + i), kvec));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[sig[c * lane_count + i] & 0x0f]
                                      [coeff[r * cols + c] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

/* P*S 的 mtk2 双输出版本。
 * 输入：两个连续 row 共享同一个 sig 切片。
 * 输出：一趟 shuffle 同时得到两行结果。
 * 场景：rank-pair / 双行更新。
 */
void pqc_kernel_gf16_snova_ps_mtk2_pair(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                             size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    size_t r = 0;
    for (; r + 1 < rows; r += 2) {
        uint8_t *dst0 = out + r * lane_count;
        uint8_t *dst1 = out + (r + 1) * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc0 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec acc1 = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t k0 = coeff[r * cols + c] & 0x0f;
                const uint8_t k1 = coeff[(r + 1) * cols + c] & 0x0f;
                const pqc_kernel_gf16_vec packed_prod = pqc_kernel_gf16_vec_mul_dual_coeff_shuffle(
                    pqc_kernel_gf16_vec_load(sig + c * lane_count + i), k0, k1);
                acc0 = pqc_kernel_gf16_vec_xor(acc0, pqc_kernel_gf16_vec_mask_u8(packed_prod, 0x0f));
                acc1 = pqc_kernel_gf16_vec_xor(acc1, pqc_kernel_gf16_vec_srli4(packed_prod));
            }
            pqc_kernel_gf16_vec_store(dst0 + i, acc0);
            pqc_kernel_gf16_vec_store(dst1 + i, acc1);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc0 = 0;
            uint8_t acc1 = 0;
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t x = sig[c * lane_count + i] & 0x0f;
                acc0 ^= pqc_kernel_gf16_mul_table[x][coeff[r * cols + c] & 0x0f];
                acc1 ^= pqc_kernel_gf16_mul_table[x][coeff[(r + 1) * cols + c] & 0x0f];
            }
            dst0[i] = acc0;
            dst1[i] = acc1;
        }
    }
    if (r < rows) {
        pqc_kernel_gf16_snova_ps_shuffle(out + r * lane_count, coeff + r * cols, sig, 1, cols, lane_count);
    }
}

/* P*S 的expanded 对照。
 * 输入：coeff / sig 与其它版本相同。
 * 输出：普通 byte GF16 结果。
 * 特点：每个系数都把 sig 展开成 uint16_t，再做延迟压缩式累加。
 */
void pqc_kernel_gf16_snova_ps_official(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                            size_t rows, size_t cols, size_t lane_count) {
    uint16_t *expanded_sig = (uint16_t *)malloc(lane_count * sizeof(uint16_t));
    uint16_t *acc = (uint16_t *)malloc(lane_count * sizeof(uint16_t));
    if (!expanded_sig || !acc) {
        free(expanded_sig);
        free(acc);
        return;
    }
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        memset(acc, 0, lane_count * sizeof(uint16_t));
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *src = sig + c * lane_count;
            uint8_t k = coeff[r * cols + c] & 0x0f;
            pqc_kernel_gf16_official_expand_array(expanded_sig, src, lane_count);
            pqc_kernel_gf16_official_const_expanded_muladd(acc, k, expanded_sig, lane_count);
        }
        pqc_kernel_gf16_official_compress_array(out + r * lane_count, acc, lane_count);
    }
    free(expanded_sig);
    free(acc);
}

/* S^T*T 的标量 table 版本。
 * 输入：left 是 rows x cols 的左矩阵块，每个块内部按 lane_count 排列；right 是 cols 组右向量。
 * 输出：out 是 rows 组收缩结果。
 * 场景：变量乘变量的普通基线。
 */
void pqc_kernel_gf16_snova_sts_table(uint8_t *out, const uint8_t *left, const uint8_t *right,
                          size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        const uint8_t *lrow = left + r * cols * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *a = lrow + c * lane_count;
            const uint8_t *b = right + c * lane_count;
            for (size_t i = 0; i < lane_count; ++i) {
                dst[i] ^= pqc_kernel_gf16_mul_table[a[i] & 0x0f][b[i] & 0x0f];
            }
        }
    }
}

/* S^T*T 的标量 log/exp 版本。 */
void pqc_kernel_gf16_snova_sts_modexp_scalar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                                  size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        const uint8_t *lrow = left + r * cols * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *a = lrow + c * lane_count;
            const uint8_t *b = right + c * lane_count;
            for (size_t i = 0; i < lane_count; ++i) {
                dst[i] ^= pqc_kernel_gf16_mul_modexp_scalar(a[i], b[i]);
            }
        }
    }
}

/* S^T*T 的 bit-plane 版本。
 * 特点：右侧向量按 b,2b,4b,8b 预计算后再与左侧逐 lane 组合。
 */
void pqc_kernel_gf16_snova_sts_bitplane(uint8_t *out, const uint8_t *left, const uint8_t *right,
                             size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        const uint8_t *lrow = left + r * cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc, pqc_kernel_gf16_vec_mul_bitplane(
                    pqc_kernel_gf16_vec_load(lrow + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right + c * lane_count + i)));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[lrow[c * lane_count + i] & 0x0f]
                                      [right[c * lane_count + i] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

void pqc_kernel_gf16_snova_sts_bitplane_ilp(uint8_t *out, const uint8_t *left, const uint8_t *right,
                                 size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        const uint8_t *lrow = left + r * cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc, pqc_kernel_gf16_vec_mul_bitplane_ilp(
                    pqc_kernel_gf16_vec_load(lrow + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right + c * lane_count + i)));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_table[lrow[c * lane_count + i] & 0x0f]
                                      [right[c * lane_count + i] & 0x0f];
            }
            dst[i] = acc;
        }
    }
}

/* S^T*T 的 same-b 四行路径。
 * 对同一个右向量只生成一次 b,2b,4b,8b，再同时更新四个输出行。
 * 向量累加器按 lane 保持到所有列完成，避免每列读写输出矩阵。
 */
void pqc_kernel_gf16_snova_sts_same_b4(uint8_t *out, const uint8_t *left, const uint8_t *right,
                            size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    size_t r = 0;
    for (; r + 4u <= rows; r += 4u) {
        const uint8_t *l0 = left + r * cols * lane_count;
        const uint8_t *l1 = l0 + cols * lane_count;
        const uint8_t *l2 = l1 + cols * lane_count;
        const uint8_t *l3 = l2 + cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec a0 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec a1 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec a2 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec a3 = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                pqc_kernel_gf16_vec p0, p1, p2, p3;
                pqc_kernel_gf16_vec_mul_same_b4_ilp(
                    pqc_kernel_gf16_vec_load(l0 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(l1 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(l2 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(l3 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right + c * lane_count + i),
                    &p0, &p1, &p2, &p3);
                a0 = pqc_kernel_gf16_vec_xor(a0, p0);
                a1 = pqc_kernel_gf16_vec_xor(a1, p1);
                a2 = pqc_kernel_gf16_vec_xor(a2, p2);
                a3 = pqc_kernel_gf16_vec_xor(a3, p3);
            }
            pqc_kernel_gf16_vec_store(out + r * lane_count + i, a0);
            pqc_kernel_gf16_vec_store(out + (r + 1u) * lane_count + i, a1);
            pqc_kernel_gf16_vec_store(out + (r + 2u) * lane_count + i, a2);
            pqc_kernel_gf16_vec_store(out + (r + 3u) * lane_count + i, a3);
        }
        for (; i < lane_count; ++i) {
            uint8_t a0 = 0, a1 = 0, a2 = 0, a3 = 0;
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t b = right[c * lane_count + i] & 0x0f;
                a0 ^= pqc_kernel_gf16_mul_table[l0[c * lane_count + i] & 0x0f][b];
                a1 ^= pqc_kernel_gf16_mul_table[l1[c * lane_count + i] & 0x0f][b];
                a2 ^= pqc_kernel_gf16_mul_table[l2[c * lane_count + i] & 0x0f][b];
                a3 ^= pqc_kernel_gf16_mul_table[l3[c * lane_count + i] & 0x0f][b];
            }
            out[r * lane_count + i] = a0;
            out[(r + 1u) * lane_count + i] = a1;
            out[(r + 2u) * lane_count + i] = a2;
            out[(r + 3u) * lane_count + i] = a3;
        }
    }
    if (r < rows) {
        pqc_kernel_gf16_snova_sts_bitplane_ilp(out + r * lane_count,
                                    left + r * cols * lane_count,
                                    right, rows - r, cols, lane_count);
    }
}

/* S^T*T 的全 packed 路径（矩阵场景）：left/right 都是 packed，输出 packed。
 * 只算直配 [a_hi*b_hi | a_lo*b_lo]，不产生交叉项；中间用 split 拆两路各自
 * bit-plane 乘，再拼回 packed，全程不拆包。实测比 cmp（含交叉）快约 2.3x、
 * 比非 packed 快约 1.6x（见文档 10.6）。
 */
void pqc_kernel_gf16_snova_sts_packed_ab_packed(uint8_t *out_packed, const uint8_t *left_packed,
                                     const uint8_t *right_packed,
                                     size_t rows, size_t cols, size_t lane_count) {
    memset(out_packed, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out_packed + r * lane_count;
        const uint8_t *lrow = left_packed + r * cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc, pqc_kernel_gf16_vec_mul_packed_ab_straight_split(
                    pqc_kernel_gf16_vec_load(lrow + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right_packed + c * lane_count + i)));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t lo = 0, hi = 0;
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t a = lrow[c * lane_count + i];
                const uint8_t b = right_packed[c * lane_count + i];
                lo ^= pqc_kernel_gf16_mul_table[a & 0x0f][b & 0x0f];
                hi ^= pqc_kernel_gf16_mul_table[(a >> 4) & 0x0f][(b >> 4) & 0x0f];
            }
            dst[i] = (uint8_t)((lo & 0x0f) | ((hi & 0x0f) << 4));
        }
    }
}

/* S^T*T 的全 packed 路径，一次 pass 同时产出直配 + 交叉两条 packed 结果。 */
void pqc_kernel_gf16_snova_sts_packed_ab_dual(uint8_t *out_straight_packed, uint8_t *out_cross_packed,
                                   const uint8_t *left_packed, const uint8_t *right_packed,
                                   size_t rows, size_t cols, size_t lane_count) {
    memset(out_straight_packed, 0, rows * lane_count);
    memset(out_cross_packed, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst_s = out_straight_packed + r * lane_count;
        uint8_t *dst_c = out_cross_packed + r * lane_count;
        const uint8_t *lrow = left_packed + r * cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc_s = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec acc_c = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                pqc_kernel_gf16_vec s, c2;
                pqc_kernel_gf16_vec_mul_packed_ab_dual_packed(
                    pqc_kernel_gf16_vec_load(lrow + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right_packed + c * lane_count + i), &s, &c2);
                acc_s = pqc_kernel_gf16_vec_xor(acc_s, s);
                acc_c = pqc_kernel_gf16_vec_xor(acc_c, c2);
            }
            pqc_kernel_gf16_vec_store(dst_s + i, acc_s);
            pqc_kernel_gf16_vec_store(dst_c + i, acc_c);
        }
        for (; i < lane_count; ++i) {
            uint8_t slo = 0, shi = 0, clo = 0, chi = 0;
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t a = lrow[c * lane_count + i];
                const uint8_t b = right_packed[c * lane_count + i];
                const uint8_t a_lo = a & 0x0f, a_hi = (a >> 4) & 0x0f;
                const uint8_t b_lo = b & 0x0f, b_hi = (b >> 4) & 0x0f;
                slo ^= pqc_kernel_gf16_mul_table[a_lo][b_lo];
                shi ^= pqc_kernel_gf16_mul_table[a_hi][b_hi];
                clo ^= pqc_kernel_gf16_mul_table[a_lo][b_hi];
                chi ^= pqc_kernel_gf16_mul_table[a_hi][b_lo];
            }
            dst_s[i] = (uint8_t)((slo & 0x0f) | ((shi & 0x0f) << 4));
            dst_c[i] = (uint8_t)((clo & 0x0f) | ((chi & 0x0f) << 4));
        }
    }
}

/* S^T*T 全 packed + 4 行复用（same-b4 的双元素版本）：
 * right packed 的平面跨 4 行只算一次，四行 packed 左输入各自选择。 */
void pqc_kernel_gf16_snova_sts_packed_ab_same_b4(uint8_t *out_packed, const uint8_t *left_packed,
                                      const uint8_t *right_packed,
                                      size_t rows, size_t cols, size_t lane_count) {
    memset(out_packed, 0, rows * lane_count);
    size_t r = 0;
    for (; r + 4u <= rows; r += 4u) {
        const uint8_t *l0 = left_packed + r * cols * lane_count;
        const uint8_t *l1 = l0 + cols * lane_count;
        const uint8_t *l2 = l1 + cols * lane_count;
        const uint8_t *l3 = l2 + cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc0 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec acc1 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec acc2 = pqc_kernel_gf16_vec_zero();
            pqc_kernel_gf16_vec acc3 = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                pqc_kernel_gf16_vec p0, p1, p2, p3;
                pqc_kernel_gf16_vec_mul_packed_ab_same_b4(
                    pqc_kernel_gf16_vec_load(l0 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(l1 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(l2 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(l3 + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right_packed + c * lane_count + i),
                    &p0, &p1, &p2, &p3);
                acc0 = pqc_kernel_gf16_vec_xor(acc0, p0);
                acc1 = pqc_kernel_gf16_vec_xor(acc1, p1);
                acc2 = pqc_kernel_gf16_vec_xor(acc2, p2);
                acc3 = pqc_kernel_gf16_vec_xor(acc3, p3);
            }
            pqc_kernel_gf16_vec_store(out_packed + r * lane_count + i, acc0);
            pqc_kernel_gf16_vec_store(out_packed + (r + 1u) * lane_count + i, acc1);
            pqc_kernel_gf16_vec_store(out_packed + (r + 2u) * lane_count + i, acc2);
            pqc_kernel_gf16_vec_store(out_packed + (r + 3u) * lane_count + i, acc3);
        }
        for (; i < lane_count; ++i) {
            uint8_t lo[4] = {0, 0, 0, 0};
            uint8_t hi[4] = {0, 0, 0, 0};
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t b_lo = right_packed[c * lane_count + i] & 0x0f;
                const uint8_t b_hi = (right_packed[c * lane_count + i] >> 4) & 0x0f;
                const uint8_t *lrows[4] = {l0, l1, l2, l3};
                for (int j = 0; j < 4; ++j) {
                    const uint8_t a = lrows[j][c * lane_count + i];
                    lo[j] ^= pqc_kernel_gf16_mul_table[a & 0x0f][b_lo];
                    hi[j] ^= pqc_kernel_gf16_mul_table[(a >> 4) & 0x0f][b_hi];
                }
            }
            for (int j = 0; j < 4; ++j) {
                out_packed[(r + (size_t)j) * lane_count + i] =
                    (uint8_t)((lo[j] & 0x0f) | ((hi[j] & 0x0f) << 4));
            }
        }
    }
    for (; r < rows; ++r) {
        uint8_t *dst = out_packed + r * lane_count;
        const uint8_t *lrow = left_packed + r * cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc, pqc_kernel_gf16_vec_mul_packed_ab_straight_split(
                    pqc_kernel_gf16_vec_load(lrow + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right_packed + c * lane_count + i)));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t lo = 0, hi = 0;
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t a = lrow[c * lane_count + i];
                const uint8_t b = right_packed[c * lane_count + i];
                lo ^= pqc_kernel_gf16_mul_table[a & 0x0f][b & 0x0f];
                hi ^= pqc_kernel_gf16_mul_table[(a >> 4) & 0x0f][(b >> 4) & 0x0f];
            }
            dst[i] = (uint8_t)((lo & 0x0f) | ((hi & 0x0f) << 4));
        }
    }
}

/* S^T*T 的 SIMD log/exp 版本。 */
void pqc_kernel_gf16_snova_sts_modexp_vec(uint8_t *out, const uint8_t *left, const uint8_t *right,
                               size_t rows, size_t cols, size_t lane_count) {
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out + r * lane_count;
        const uint8_t *lrow = left + r * cols * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                acc = pqc_kernel_gf16_vec_xor(acc, pqc_kernel_gf16_vec_mul_modexp(
                    pqc_kernel_gf16_vec_load(lrow + c * lane_count + i),
                    pqc_kernel_gf16_vec_load(right + c * lane_count + i)));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf16_mul_modexp_scalar(lrow[c * lane_count + i],
                                              right[c * lane_count + i]);
            }
            dst[i] = acc;
        }
    }
}

/* S^T*T 的expanded 对照。
 * 输入：right 先扩展成 uint16_t。
 * 输出：普通 byte GF16。
 * 场景：和verify 的延迟压缩方式对拍。
 */
void pqc_kernel_gf16_snova_sts_official(uint8_t *out, const uint8_t *left, const uint8_t *right,
                             size_t rows, size_t cols, size_t lane_count) {
    uint16_t *expanded = (uint16_t *)malloc(lane_count * sizeof(uint16_t));
    uint16_t *acc = (uint16_t *)malloc(lane_count * sizeof(uint16_t));
    if (!expanded || !acc) {
        free(expanded);
        free(acc);
        return;
    }
    memset(out, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        memset(acc, 0, lane_count * sizeof(uint16_t));
        const uint8_t *lrow = left + r * cols * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *a = lrow + c * lane_count;
            const uint8_t *b = right + c * lane_count;
            pqc_kernel_gf16_official_expand_array(expanded, b, lane_count);
            pqc_kernel_gf16_official_raw_expanded_muladd(acc, a, expanded, lane_count);
        }
        pqc_kernel_gf16_official_compress_array(out + r * lane_count, acc, lane_count);
    }
    free(expanded);
    free(acc);
}

void pqc_kernel_gf16_snova_packed_pair_contract(uint8_t *out_packed, const uint8_t *coeff_pair,
                                     const uint8_t *in_packed, size_t rows, size_t cols,
                                     size_t lane_count) {
    memset(out_packed, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst = out_packed + r * lane_count;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const size_t ci = 2u * (r * cols + c);
                acc = pqc_kernel_gf16_vec_xor(acc, pqc_kernel_gf16_vec_mul_packed_pair_shuffle(
                    pqc_kernel_gf16_vec_load(in_packed + c * lane_count + i),
                    coeff_pair[ci] & 0x0f, coeff_pair[ci + 1u] & 0x0f));
            }
            pqc_kernel_gf16_vec_store(dst + i, acc);
        }
        for (; i < lane_count; ++i) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                const size_t ci = 2u * (r * cols + c);
                const uint8_t x = in_packed[c * lane_count + i];
                acc ^= (uint8_t)(pqc_kernel_gf16_mul_table[x & 0x0f][coeff_pair[ci] & 0x0f] |
                                 (pqc_kernel_gf16_mul_table[(x >> 4) & 0x0f][coeff_pair[ci + 1u] & 0x0f] << 4));
            }
            dst[i] = acc;
        }
    }
}

void pqc_kernel_gf16_snova_packed_pair_chain(uint8_t *out_packed,
                                   const uint8_t *coeff_stage1,
                                   const uint8_t *coeff_stage2,
                                   const uint8_t *in_packed,
                                   uint8_t *scratch_packed,
                                   size_t rows, size_t cols, size_t lane_count) {
    pqc_kernel_gf16_snova_packed_pair_contract(scratch_packed, coeff_stage1, in_packed,
                                    rows, cols, lane_count);
    pqc_kernel_gf16_snova_packed_pair_contract(out_packed, coeff_stage2, scratch_packed,
                                    rows, rows, lane_count);
}

static void pqc_kernel_gf16_expanded_pair_contract(uint16_t *out_lo, uint16_t *out_hi,
                                        const uint8_t *coeff_pair,
                                        const uint16_t *in_lo, const uint16_t *in_hi,
                                        size_t rows, size_t cols, size_t lane_count) {
    memset(out_lo, 0, rows * lane_count * sizeof(uint16_t));
    memset(out_hi, 0, rows * lane_count * sizeof(uint16_t));
    for (size_t r = 0; r < rows; ++r) {
        uint16_t *dst_lo = out_lo + r * lane_count;
        uint16_t *dst_hi = out_hi + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            const size_t ci = 2u * (r * cols + c);
            pqc_kernel_gf16_official_const_expanded_muladd(dst_lo, coeff_pair[ci],
                                                in_lo + c * lane_count, lane_count);
            pqc_kernel_gf16_official_const_expanded_muladd(dst_hi, coeff_pair[ci + 1],
                                                in_hi + c * lane_count, lane_count);
        }
    }
}

void pqc_kernel_gf16_snova_expanded_pair_chain_preexpanded(uint16_t *out_lo, uint16_t *out_hi,
                                                const uint8_t *coeff_stage1,
                                                const uint8_t *coeff_stage2,
                                                const uint16_t *in_lo, const uint16_t *in_hi,
                                                uint16_t *scratch_lo, uint16_t *scratch_hi,
                                                uint8_t *raw_mid_lo, uint8_t *raw_mid_hi,
                                                size_t rows, size_t cols, size_t lane_count) {
    pqc_kernel_gf16_expanded_pair_contract(scratch_lo, scratch_hi, coeff_stage1,
                                in_lo, in_hi, rows, cols, lane_count);
    pqc_kernel_gf16_official_compress_array(raw_mid_lo, scratch_lo, rows * lane_count);
    pqc_kernel_gf16_official_compress_array(raw_mid_hi, scratch_hi, rows * lane_count);
    pqc_kernel_gf16_official_expand_array(scratch_lo, raw_mid_lo, rows * lane_count);
    pqc_kernel_gf16_official_expand_array(scratch_hi, raw_mid_hi, rows * lane_count);
    pqc_kernel_gf16_expanded_pair_contract(out_lo, out_hi, coeff_stage2,
                                scratch_lo, scratch_hi, rows, rows, lane_count);
}

void pqc_kernel_gf16_snova_expanded_pair_chain(uint16_t *out_lo, uint16_t *out_hi,
                                    const uint8_t *coeff_stage1,
                                    const uint8_t *coeff_stage2,
                                    const uint8_t *in_lo, const uint8_t *in_hi,
                                    uint16_t *expanded_in_lo, uint16_t *expanded_in_hi,
                                    uint16_t *scratch_lo, uint16_t *scratch_hi,
                                    uint8_t *raw_mid_lo, uint8_t *raw_mid_hi,
                                    size_t rows, size_t cols, size_t lane_count) {
    pqc_kernel_gf16_official_expand_array(expanded_in_lo, in_lo, cols * lane_count);
    pqc_kernel_gf16_official_expand_array(expanded_in_hi, in_hi, cols * lane_count);
    pqc_kernel_gf16_snova_expanded_pair_chain_preexpanded(out_lo, out_hi,
                                               coeff_stage1, coeff_stage2,
                                               expanded_in_lo, expanded_in_hi,
                                               scratch_lo, scratch_hi, raw_mid_lo, raw_mid_hi,
                                               rows, cols, lane_count);
}

void pqc_kernel_gf16_snova_packeda_split_ilp(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                                  size_t rows, size_t cols, size_t lane_count) {
    memset(out_lo, 0, rows * lane_count);
    memset(out_hi, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst_lo = out_lo + r * lane_count;
        uint8_t *dst_hi = out_hi + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            pqc_kernel_gf16_vec b1 = pqc_kernel_gf16_vec_set1(k);
            pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b1);
            pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b1);
            pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b1);
            uint8_t b1_bytes[GF16_VEC_BYTES];
            uint8_t b2_bytes[GF16_VEC_BYTES];
            uint8_t b4_bytes[GF16_VEC_BYTES];
            uint8_t b8_bytes[GF16_VEC_BYTES];
            pqc_kernel_gf16_vec_store(b1_bytes, b1);
            pqc_kernel_gf16_vec_store(b2_bytes, b2);
            pqc_kernel_gf16_vec_store(b4_bytes, b4);
            pqc_kernel_gf16_vec_store(b8_bytes, b8);
            const uint8_t *src = sig_packed + c * lane_count;
            size_t i = 0;
            for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
                uint8_t lo[GF16_VEC_BYTES];
                uint8_t hi[GF16_VEC_BYTES];
                pqc_kernel_gf16_mul_packed_a_split_dual_array_ilp(lo, hi, src + i, b1_bytes, b2_bytes, b4_bytes, b8_bytes,
                                                       GF16_VEC_BYTES);
                pqc_kernel_gf16_vec_store(dst_lo + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_lo + i), pqc_kernel_gf16_vec_load(lo)));
                pqc_kernel_gf16_vec_store(dst_hi + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_hi + i), pqc_kernel_gf16_vec_load(hi)));
            }
            for (; i < lane_count; ++i) {
                uint8_t x = src[i];
                dst_lo[i] ^= pqc_kernel_gf16_mul_table[x & 0x0f][k];
                dst_hi[i] ^= pqc_kernel_gf16_mul_table[(x >> 4) & 0x0f][k];
            }
        }
    }
}

void pqc_kernel_gf16_snova_packeda_cmp_ilp(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                                size_t rows, size_t cols, size_t lane_count) {
    memset(out_lo, 0, rows * lane_count);
    memset(out_hi, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst_lo = out_lo + r * lane_count;
        uint8_t *dst_hi = out_hi + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            pqc_kernel_gf16_vec b1 = pqc_kernel_gf16_vec_set1(k);
            pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b1);
            pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b1);
            pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b1);
            uint8_t b1_bytes[GF16_VEC_BYTES];
            uint8_t b2_bytes[GF16_VEC_BYTES];
            uint8_t b4_bytes[GF16_VEC_BYTES];
            uint8_t b8_bytes[GF16_VEC_BYTES];
            pqc_kernel_gf16_vec_store(b1_bytes, b1);
            pqc_kernel_gf16_vec_store(b2_bytes, b2);
            pqc_kernel_gf16_vec_store(b4_bytes, b4);
            pqc_kernel_gf16_vec_store(b8_bytes, b8);
            const uint8_t *src = sig_packed + c * lane_count;
            size_t i = 0;
            for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
                uint8_t lo[GF16_VEC_BYTES];
                uint8_t hi[GF16_VEC_BYTES];
                pqc_kernel_gf16_mul_packed_a_cmp_dual_array_ilp(lo, hi, src + i, b1_bytes, b2_bytes, b4_bytes, b8_bytes,
                                                     GF16_VEC_BYTES);
                pqc_kernel_gf16_vec_store(dst_lo + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_lo + i), pqc_kernel_gf16_vec_load(lo)));
                pqc_kernel_gf16_vec_store(dst_hi + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_hi + i), pqc_kernel_gf16_vec_load(hi)));
            }
            for (; i < lane_count; ++i) {
                uint8_t x = src[i];
                dst_lo[i] ^= pqc_kernel_gf16_mul_table[x & 0x0f][k];
                dst_hi[i] ^= pqc_kernel_gf16_mul_table[(x >> 4) & 0x0f][k];
            }
        }
    }
}

/* packed-a 常数系数收缩。
 * 输入：coeff 是常数系数矩阵；sig_packed 的低/高 nibble 是两路输入。
 * 输出：out_lo / out_hi 分别存两路结果。
 * 场景：packed 左输入的常数乘路径。
 */
void pqc_kernel_gf16_snova_packeda_table(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                              size_t rows, size_t cols, size_t lane_count) {
    memset(out_lo, 0, rows * lane_count);
    memset(out_hi, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst_lo = out_lo + r * lane_count;
        uint8_t *dst_hi = out_hi + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            const uint8_t *src = sig_packed + c * lane_count;
            for (size_t i = 0; i < lane_count; ++i) {
                uint8_t x = src[i];
                dst_lo[i] ^= pqc_kernel_gf16_mul_table[x & 0x0f][k];
                dst_hi[i] ^= pqc_kernel_gf16_mul_table[(x >> 4) & 0x0f][k];
            }
        }
    }
}

/* packed-a 常数系数收缩的 split 版本。
 * 特点：左输入按 low/high 两路分开做 bit-plane 选择，再分别累加。
 */
void pqc_kernel_gf16_snova_packeda_split(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                              size_t rows, size_t cols, size_t lane_count) {
    memset(out_lo, 0, rows * lane_count);
    memset(out_hi, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst_lo = out_lo + r * lane_count;
        uint8_t *dst_hi = out_hi + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            pqc_kernel_gf16_vec b1 = pqc_kernel_gf16_vec_set1(k);
            pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b1);
            pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b1);
            pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b1);
            const uint8_t *src = sig_packed + c * lane_count;
            size_t i = 0;
            for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
                pqc_kernel_gf16_vec lo, hi;
                pqc_kernel_gf16_vec_mul_packed_a_split_dual(pqc_kernel_gf16_vec_load(src + i), b1, b2, b4, b8, &lo, &hi);
                pqc_kernel_gf16_vec_store(dst_lo + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_lo + i), lo));
                pqc_kernel_gf16_vec_store(dst_hi + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_hi + i), hi));
            }
            for (; i < lane_count; ++i) {
                uint8_t x = src[i];
                dst_lo[i] ^= pqc_kernel_gf16_mul_table[x & 0x0f][k];
                dst_hi[i] ^= pqc_kernel_gf16_mul_table[(x >> 4) & 0x0f][k];
            }
        }
    }
}

/* packed-a 常数系数收缩的 cmp 版本。
 * 特点：用 compare-mask 直接生成选择掩码，作为 split 的对照实现。
 */
void pqc_kernel_gf16_snova_packeda_cmp(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                            size_t rows, size_t cols, size_t lane_count) {
    memset(out_lo, 0, rows * lane_count);
    memset(out_hi, 0, rows * lane_count);
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *dst_lo = out_lo + r * lane_count;
        uint8_t *dst_hi = out_hi + r * lane_count;
        for (size_t c = 0; c < cols; ++c) {
            uint8_t k = coeff[r * cols + c] & 0x0f;
            pqc_kernel_gf16_vec b1 = pqc_kernel_gf16_vec_set1(k);
            pqc_kernel_gf16_vec b2 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul2_vec, b1);
            pqc_kernel_gf16_vec b4 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul4_vec, b1);
            pqc_kernel_gf16_vec b8 = pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_mul8_vec, b1);
            const uint8_t *src = sig_packed + c * lane_count;
            size_t i = 0;
            for (; i + GF16_VEC_BYTES <= lane_count; i += GF16_VEC_BYTES) {
                pqc_kernel_gf16_vec lo, hi;
                pqc_kernel_gf16_vec_mul_packed_a_cmp_dual(pqc_kernel_gf16_vec_load(src + i), b1, b2, b4, b8, &lo, &hi);
                pqc_kernel_gf16_vec_store(dst_lo + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_lo + i), lo));
                pqc_kernel_gf16_vec_store(dst_hi + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(dst_hi + i), hi));
            }
            for (; i < lane_count; ++i) {
                uint8_t x = src[i];
                dst_lo[i] ^= pqc_kernel_gf16_mul_table[x & 0x0f][k];
                dst_hi[i] ^= pqc_kernel_gf16_mul_table[(x >> 4) & 0x0f][k];
            }
        }
    }
}
