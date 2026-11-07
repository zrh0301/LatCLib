#ifndef GF16_H
#define GF16_H

#include "gf16_platform.h"

#include <stdint.h>
#include <stddef.h>

extern uint8_t pqc_kernel_gf16_mul_table[16][16];
extern uint8_t pqc_kernel_gf16_log_table[16];
extern uint8_t pqc_kernel_gf16_exp_table[15];
extern pqc_kernel_gf16_vec pqc_kernel_gf16_mtk2_16[256];
extern pqc_kernel_gf16_vec pqc_kernel_gf16_mtk_uniform[16];
extern pqc_kernel_gf16_vec pqc_kernel_gf16_log_vec;
extern pqc_kernel_gf16_vec pqc_kernel_gf16_exp_vec;
extern pqc_kernel_gf16_vec pqc_kernel_gf16_mul2_vec;
extern pqc_kernel_gf16_vec pqc_kernel_gf16_mul4_vec;
extern pqc_kernel_gf16_vec pqc_kernel_gf16_mul8_vec;

/* 初始化所有 GF16 查表和向量常量。
 * 输入：无。
 * 输出：填充普通 16x16 乘法表、log/exp 表，以及 shuffle/bit-plane 预计算向量。
 * 场景：程序启动时调用一次，之后所有算法族共享这些常量。
 */
void pqc_kernel_gf16_tables_init(void);

/* 标量 GF16 乘法。
 * 输入：a、b 的低 4 bit 是 GF16 元素。
 * 输出：单个 GF16 结果的低 4 bit。
 * 场景：正确性参考、标量回退、小数据控制流。
 */
uint8_t pqc_kernel_gf16_mul_ref(uint8_t a, uint8_t b);
uint8_t pqc_kernel_gf16_mul_table_scalar(uint8_t a, uint8_t b);
uint8_t pqc_kernel_gf16_mul_modexp_scalar(uint8_t a, uint8_t b);

/* 标量数组版本。
 * 输入：a/b 为普通 byte GF16 数组；n 为元素个数。
 * 输出：out[i] = a[i] * b[i]。
 * 场景：最朴素基线、与 SIMD 原语做同口径对照。
 */
void pqc_kernel_gf16_mul_table_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf16_mul_modexp_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf16_mul_const_table_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf16_mul_const_modexp_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);

/* SIMD shuffle 常数乘。
 * 输入：a 为 byte GF16 向量；k 为常数系数。
 * 输出：每个 lane 都乘以 k。
 * 场景：固定系数矩阵乘法、SNOVA 的 P*S、A/B/q 固定块作用。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_shuffle(pqc_kernel_gf16_vec a, uint8_t k);
/* 一条向量同时输出两个常数系数乘积。
 * 输入：a 为普通 byte GF16 向量；low_k/high_k 为两个常数系数。
 * 输出：每个 lane 的低/高半字节分别承载两个结果。
 * 场景：rank-pair / 双行复用。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_dual_coeff_shuffle(pqc_kernel_gf16_vec a, uint8_t low_k, uint8_t high_k);
/* packed pair 版本。
 * 输入：packed byte 内低/高 nibble 是两个不同 GF16 元素。
 * 输出：低 nibble 乘 low_k，高 nibble 乘 high_k。
 * 场景：已经 packed 的输入需要保留 packed 输出。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_pair_shuffle(pqc_kernel_gf16_vec packed, uint8_t low_k, uint8_t high_k);
void pqc_kernel_gf16_mul_const_shuffle_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf16_mul_dual_coeff_shuffle_array(uint8_t *out, const uint8_t *a, uint8_t low_k, uint8_t high_k, size_t n);
void pqc_kernel_gf16_mul_packed_pair_shuffle_array(uint8_t *out, const uint8_t *packed, uint8_t low_k, uint8_t high_k, size_t n);

/* SIMD log/exp 乘法。
 * 输入：a、b 为 byte GF16 向量。
 * 输出：逐 lane 乘积。
 * 场景：变量乘变量、希望保留普通 byte 布局的场景。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_const_modexp(pqc_kernel_gf16_vec a, uint8_t k);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_modexp(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b);
void pqc_kernel_gf16_mul_const_modexp_vec_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf16_mul_modexp_vec_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);

/* bit-plane 乘法。
 * 输入：a、b 为普通 byte GF16 向量；内部会预计算 b,2b,4b,8b。
 * 输出：逐 lane 乘积。
 * 场景：变量乘变量、同一 b 复用多次、same_b4 / same_packed_b4。
 */
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_bitplane(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_bitplane_ilp(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b);
void pqc_kernel_gf16_vec_mul_same_b4_ilp(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                              pqc_kernel_gf16_vec b, pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1,
                              pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3);
void pqc_kernel_gf16_mul_bitplane_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf16_mul_bitplane_array_ilp(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf16_mul_bitplane_same_b4_array(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                     const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                     const uint8_t *b, size_t n);
void pqc_kernel_gf16_mul_bitplane_same_b4_array_ilp(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                         const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                         const uint8_t *b, size_t n);
/* packed-a 双输出。
 * 输入：a_packed 的低/高 nibble 是两个不同 GF16 元素；b1/b2/b4/b8 是右侧预计算表。
 * 输出：out_lo = a_lo * b，out_hi = a_hi * b。
 * 场景：左侧长期保持 packed 的路径；split 是默认候选，cmp 是对照路径。
 */
void pqc_kernel_gf16_vec_mul_packed_a_split_dual(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8,
                                      pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi);
void pqc_kernel_gf16_vec_mul_packed_a_cmp_dual(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b1, pqc_kernel_gf16_vec b2, pqc_kernel_gf16_vec b4, pqc_kernel_gf16_vec b8,
                                    pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi);
/* packed-a/packed-b 变量乘法。
 * 输入两个 packed byte：各自的 low/high nibble 是独立 GF(16) 元素。
 * 输出两个普通向量，分别为 low*low 与 high*high。
 * 这是 verify 长链中保持两路数据 packed 的核心变量乘法。
 */
void pqc_kernel_gf16_vec_mul_packed_ab_split_dual(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed,
                                       pqc_kernel_gf16_vec *out_lo, pqc_kernel_gf16_vec *out_hi);
/* 全 packed 变量乘：A_p × B_p -> packed 输出（不拆包、不重打包）。 */
void pqc_kernel_gf16_vec_mul_packed_ab_dual_packed(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed,
                                        pqc_kernel_gf16_vec *out_straight, pqc_kernel_gf16_vec *out_cross);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_ab_straight(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_ab_straight_ilp(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mul_packed_ab_straight_split(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed);
void pqc_kernel_gf16_vec_mul_packed_ab_dual_packed_ilp(pqc_kernel_gf16_vec a_packed, pqc_kernel_gf16_vec b_packed,
                                            pqc_kernel_gf16_vec *out_straight, pqc_kernel_gf16_vec *out_cross);
void pqc_kernel_gf16_mul_packed_ab_straight_array(uint8_t *out, const uint8_t *a_packed,
                                       const uint8_t *b_packed, size_t n);
void pqc_kernel_gf16_mul_packed_ab_straight_array_ilp(uint8_t *out, const uint8_t *a_packed,
                                           const uint8_t *b_packed, size_t n);
void pqc_kernel_gf16_mul_packed_ab_dual_packed_array(uint8_t *out_straight, uint8_t *out_cross,
                                          const uint8_t *a_packed, const uint8_t *b_packed,
                                          size_t n);
void pqc_kernel_gf16_vec_mul_packed_ab_same_b4(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                                    pqc_kernel_gf16_vec b_packed,
                                    pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1, pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3);
void pqc_kernel_gf16_vec_mul_packed_ab_same_b4_ilp(pqc_kernel_gf16_vec a0, pqc_kernel_gf16_vec a1, pqc_kernel_gf16_vec a2, pqc_kernel_gf16_vec a3,
                                        pqc_kernel_gf16_vec b_packed,
                                        pqc_kernel_gf16_vec *r0, pqc_kernel_gf16_vec *r1, pqc_kernel_gf16_vec *r2, pqc_kernel_gf16_vec *r3);
void pqc_kernel_gf16_mul_packed_a_split_dual_array(uint8_t *out_lo, uint8_t *out_hi,
                                        const uint8_t *a_packed,
                                        const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                        size_t n);
void pqc_kernel_gf16_mul_packed_a_split_dual_array_ilp(uint8_t *out_lo, uint8_t *out_hi,
                                            const uint8_t *a_packed,
                                            const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                            size_t n);
void pqc_kernel_gf16_mul_packed_a_cmp_dual_array(uint8_t *out_lo, uint8_t *out_hi,
                                      const uint8_t *a_packed,
                                      const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                      size_t n);
void pqc_kernel_gf16_mul_packed_a_cmp_dual_array_ilp(uint8_t *out_lo, uint8_t *out_hi,
                                          const uint8_t *a_packed,
                                          const uint8_t *b1, const uint8_t *b2, const uint8_t *b4, const uint8_t *b8,
                                          size_t n);
void pqc_kernel_gf16_mul_bitplane_same_packed_b4_array(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                            const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                            const uint8_t *b_packed, size_t n);
void pqc_kernel_gf16_mul_bitplane_same_packed_b4_ilp_array(uint8_t *s0, uint8_t *s1, uint8_t *s2, uint8_t *s3,
                                                 const uint8_t *a0, const uint8_t *a1, const uint8_t *a2, const uint8_t *a3,
                                                 const uint8_t *b_packed, size_t n);

/* GF(16) 求逆（正文 3.9 节）。
 * 三条路线：16 项逆元表、复用 log/exp 表、费马幂 a^14 = a^2 * a^4 * a^8。
 * 标量版按元素查表；SIMD 版一次 shuffle 处理整向量（表在 16 项索引空间内）。
 * 约定 inv(0) = 0，因此不需要额外的非零掩码。
 */
extern const uint8_t pqc_kernel_gf16_inv16_table[16];
extern pqc_kernel_gf16_vec pqc_kernel_gf16_inv_vec;
uint8_t pqc_kernel_gf16_inv_table_scalar(uint8_t a);
uint8_t pqc_kernel_gf16_inv_logexp_scalar(uint8_t a);
uint8_t pqc_kernel_gf16_inv_fermat_scalar(uint8_t a);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_inv_table_shuffle(pqc_kernel_gf16_vec a);
pqc_kernel_gf16_vec pqc_kernel_gf16_vec_inv_logexp(pqc_kernel_gf16_vec a);
void pqc_kernel_gf16_inv_table_shuffle_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf16_inv_logexp_vec_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf16_inv_fermat_array(uint8_t *out, const uint8_t *a, size_t n);
/* CLMUL/PMULL packed GF(16) backend.
 * Two logical values occupy the low two bytes of one 64-bit word; the upper
 * six bytes are zero. Cross terms are discarded before GF(16) reduction.
 * Pack/unpack are explicit so conversion cost can be timed separately. */
size_t pqc_kernel_gf16_clmul_packed_bytes(size_t n);
void pqc_kernel_gf16_clmul_pack_array(uint8_t *packed, const uint8_t *raw, size_t n);
void pqc_kernel_gf16_clmul_unpack_array(uint8_t *raw, const uint8_t *packed, size_t n);
void pqc_kernel_gf16_clmul_product_array(uint8_t *product, const uint8_t *a_packed,
                              const uint8_t *b_packed, size_t n);
void pqc_kernel_gf16_clmul_reduce_array(uint8_t *out_packed, const uint8_t *product,
                             size_t n);
void pqc_kernel_gf16_clmul_mul_packed_array(uint8_t *out_packed, const uint8_t *a_packed,
                                 const uint8_t *b_packed, size_t n);
void pqc_kernel_gf16_clmul_ps_packed(uint8_t *out_packed, const uint64_t *coeff_words,
                          const uint8_t *sig_packed, size_t rows, size_t cols,
                          size_t lanes);
void pqc_kernel_gf16_clmul_sts_packed(uint8_t *out_packed, const uint8_t *left_packed,
                           const uint8_t *right_packed, size_t rows,
                           size_t cols, size_t lanes);

/* SNOVA 风格 expanded 表示。
 * 输入：raw byte GF16 或 uint16_t 扩展态。
 * 输出：expand / compress / raw×expanded muladd。
 * 场景：模拟verify 的延迟压缩模型。
 */
uint16_t pqc_kernel_gf16_official_expand(uint8_t a);
uint8_t pqc_kernel_gf16_official_compress(uint16_t a);
void pqc_kernel_gf16_official_expand_array(uint16_t *out, const uint8_t *in, size_t n);
void pqc_kernel_gf16_official_compress_array(uint8_t *out, const uint16_t *in, size_t n);
void pqc_kernel_gf16_official_raw_expanded_muladd(uint16_t *acc, const uint8_t *raw, const uint16_t *expanded, size_t n);
void pqc_kernel_gf16_official_const_expanded_muladd(uint16_t *acc, uint8_t k, const uint16_t *expanded, size_t n);

/* SNOVA-like workload wrappers。
 * 这些不是新的数学算法，而是把上面的 GF16 原语拼成更接近 SNOVA 的局部任务。
 * P*S：常数系数矩阵/向量收缩。
 * S^T*T：变量/变量收缩。
 * packed-a：packed 左输入的 split/cmp 探测。
 */
void pqc_kernel_gf16_snova_ps_table(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                         size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_modexp_scalar(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                                 size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_shuffle(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                           size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_modexp_vec(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                              size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_bitplane(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                            size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_bitplane_ilp(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                                size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_mtk2_pair(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                              size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_ps_official(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                            size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_packeda_table(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                              size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_packeda_split(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                              size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_packeda_split_ilp(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                                  size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_packeda_cmp(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                            size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_packeda_cmp_ilp(uint8_t *out_lo, uint8_t *out_hi, const uint8_t *coeff, const uint8_t *sig_packed,
                                size_t rows, size_t cols, size_t lane_count);

/* packed / raw 交叉场景。
 * 输入：left/right 为原始矩阵式工作负载。
 * 输出：out 为收缩结果。
 * 场景：把常数/变量乘法原语放进 SNOVA 风格任务里做横向比较。
 */
void pqc_kernel_gf16_snova_sts_table(uint8_t *out, const uint8_t *left, const uint8_t *right,
                          size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_modexp_scalar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                                  size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_bitplane(uint8_t *out, const uint8_t *left, const uint8_t *right,
                             size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_bitplane_ilp(uint8_t *out, const uint8_t *left, const uint8_t *right,
                                 size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_same_b4(uint8_t *out, const uint8_t *left, const uint8_t *right,
                            size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_modexp_vec(uint8_t *out, const uint8_t *left, const uint8_t *right,
                                size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_official(uint8_t *out, const uint8_t *left, const uint8_t *right,
                              size_t rows, size_t cols, size_t lane_count);
/* 全 packed S^T*T：left/right 都是 packed，输出也是 packed，全程不拆包。 */
void pqc_kernel_gf16_snova_sts_packed_ab_packed(uint8_t *out_packed, const uint8_t *left_packed,
                                     const uint8_t *right_packed,
                                     size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_sts_packed_ab_dual(uint8_t *out_straight_packed, uint8_t *out_cross_packed,
                                   const uint8_t *left_packed, const uint8_t *right_packed,
                                   size_t rows, size_t cols, size_t lane_count);

/* 两路 packed 端到端收缩。
 * coeff_pair 的每个系数项为 [low_stream_coeff, high_stream_coeff]；
 * in_packed/out_packed 的一个 byte 低/高 nibble 分别承载两路 GF(16) 数据。
 * 该接口不做 pack/unpack，只负责保持 packed 状态完成一次收缩。
 */
void pqc_kernel_gf16_snova_packed_pair_contract(uint8_t *out_packed, const uint8_t *coeff_pair,
                                     const uint8_t *in_packed, size_t rows, size_t cols,
                                     size_t lane_count);

/* 两级 packed 链：中间结果始终保持 packed。 */
void pqc_kernel_gf16_snova_packed_pair_chain(uint8_t *out_packed,
                                   const uint8_t *coeff_stage1,
                                   const uint8_t *coeff_stage2,
                                   const uint8_t *in_packed,
                                   uint8_t *scratch_packed,
                                   size_t rows, size_t cols, size_t lane_count);

/* 两路 expanded 端到端收缩。
 * 输入是两路普通 GF(16) byte，内部扩展后保持 uint16_t expanded，
 * 每个收缩边界按收缩流程 compress，再进入下一级；函数最终输出 expanded，
 * 由调用方决定是否把最终结果 compress 回普通 GF(16)。
 */
void pqc_kernel_gf16_snova_expanded_pair_chain(uint16_t *out_lo, uint16_t *out_hi,
                                    const uint8_t *coeff_stage1,
                                    const uint8_t *coeff_stage2,
                                    const uint8_t *in_lo, const uint8_t *in_hi,
                                    uint16_t *expanded_in_lo, uint16_t *expanded_in_hi,
                                    uint16_t *scratch_lo, uint16_t *scratch_hi,
                                    uint8_t *raw_mid_lo, uint8_t *raw_mid_hi,
                                    size_t rows, size_t cols, size_t lane_count);
void pqc_kernel_gf16_snova_expanded_pair_chain_preexpanded(uint16_t *out_lo, uint16_t *out_hi,
                                                const uint8_t *coeff_stage1,
                                                const uint8_t *coeff_stage2,
                                                const uint16_t *in_lo, const uint16_t *in_hi,
                                                uint16_t *scratch_lo, uint16_t *scratch_hi,
                                                uint8_t *raw_mid_lo, uint8_t *raw_mid_hi,
                                                size_t rows, size_t cols, size_t lane_count);

#endif
