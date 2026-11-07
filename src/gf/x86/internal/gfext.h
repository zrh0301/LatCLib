#ifndef GFEXT_H
#define GFEXT_H

#include "gf16.h"

#include <stddef.h>
#include <stdint.h>

/*
 * The x86 extension layer keeps the existing GF16 API intact and adds the
 * other field domains used by the code-based and mixed-field experiments.
 * GF256 uses the Rijndael polynomial x^8+x^4+x^3+x+1 (0x11b).
 */
typedef pqc_kernel_gf16_vec pqc_kernel_gf256_vec;

/* GF(2) word operations. */
uint8_t pqc_kernel_gf2_parity_u64(uint64_t x);
uint8_t pqc_kernel_gf2_dot_u64(uint64_t a, uint64_t b);
void pqc_kernel_gf2_xor_words(uint64_t *dst, const uint64_t *src, size_t words);
uint8_t pqc_kernel_gf2_dot_words(const uint64_t *a, const uint64_t *b, size_t words);
void pqc_kernel_gf2_mat_vec_u64(uint64_t *out, const uint64_t *matrix_rows,
                     const uint64_t *vector, size_t rows, size_t words_per_row);

/* GF(256), polynomial basis with the 0x11b reduction polynomial. */
uint8_t pqc_kernel_gf256_mul_ref(uint8_t a, uint8_t b);
/* MQOM's default reference backend: four byte lanes packed in one uint32_t
 * (SWAR), using the same Rijndael reduction polynomial. */
uint8_t pqc_kernel_gf256_mul_swar_ref(uint8_t a, uint8_t b);
void pqc_kernel_gf256_mul_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_mul_swar_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_mul_const_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf256_mulacc_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);

int pqc_kernel_gf256_gfni_available(void);
pqc_kernel_gf256_vec pqc_kernel_gf256_vec_mul_avx2(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b);
/* Direct vector GFNI calls require GFNI hardware; array wrappers probe once. */
pqc_kernel_gf256_vec pqc_kernel_gf256_vec_mul_gfni(pqc_kernel_gf256_vec a, pqc_kernel_gf256_vec b);
void pqc_kernel_gf256_mul_gfni_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_mulacc_gfni_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_mul_const_gfni_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
uint8_t pqc_kernel_gf256_dot_gfni_array(const uint8_t *a, const uint8_t *b, size_t n);

/* GF(16) coefficient embedded into the GFNI/AES GF(256) byte field. */
uint8_t pqc_kernel_gf16_embed_gf256(uint8_t a);
void pqc_kernel_gf16_to_gf256_array(uint8_t *out, const uint8_t *in, size_t n);
void pqc_kernel_gf16_gf256_mul_array(uint8_t *out, const uint8_t *gf16,
                          const uint8_t *gf256, size_t n);
void pqc_kernel_gf16_gf256_mulacc_array(uint8_t *acc, const uint8_t *gf16,
                             const uint8_t *gf256, size_t n);
void pqc_kernel_gf16_gf256_mulacc_gfni_array(uint8_t *acc, const uint8_t *gf16,
                                  const uint8_t *gf256, size_t n);

/* MQOM v2.1.1-compatible packed GF(16)-GF(256) path.  The input contains
 * two GF(16) elements per byte (low nibble first); n is the logical element
 * count and the GF(256) vector/output contain n bytes.  The MQOM embedding is
 * kept separate from pqc_kernel_gf16_embed_gf256() because both maps are valid but are
 * different field-coordinate conventions. */
uint8_t pqc_kernel_mqom_gf16_embed_gf256(uint8_t a);
void pqc_kernel_mqom_gf16_to_gf256_packed_array(uint8_t *out, const uint8_t *packed,
                                     size_t n);
void pqc_kernel_mqom_gf16_gf256_const_mul_array(uint8_t *out, uint8_t pqc_kernel_gf16_const,
                                     const uint8_t *gf256, size_t n);
void pqc_kernel_mqom_gf16_gf256_mul_packed_array(uint8_t *out, const uint8_t *pqc_kernel_gf16_packed,
                                      const uint8_t *gf256, size_t n);
void pqc_kernel_mqom_gf16_gf256_mulacc_packed_array(uint8_t *acc, const uint8_t *pqc_kernel_gf16_packed,
                                         const uint8_t *gf256, size_t n);
uint8_t pqc_kernel_mqom_gf16_gf256_dot_packed(const uint8_t *pqc_kernel_gf16_packed,
                                   const uint8_t *gf256, size_t n);

/* ---------------------------------------------------------------- GF(256) 基础补充 */
/* 4 路 SWAR 乘法的字级原语：一个 uint32_t 内四个独立字节通道。 */
uint32_t pqc_kernel_gf256_mul4_swar_word(uint32_t x, uint32_t y);
void pqc_kernel_gf256_mul_const_swar_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf256_mulacc_swar_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);
uint8_t pqc_kernel_gf256_dot_swar_array(const uint8_t *a, const uint8_t *b, size_t n);

/* 标量/逐字节参考：求逆、平方、幂、异或。 */
uint8_t pqc_kernel_gf256_inv_logexp_ref(uint8_t a);
uint8_t pqc_kernel_gf256_pow_ref(uint8_t a, uint32_t e);
void pqc_kernel_gf256_inv_logexp_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_pow_scalar_array(uint8_t *out, const uint8_t *a, uint32_t e, size_t n);
void pqc_kernel_gf256_square_ref_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_xor_scalar_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_xor_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);

/* ---------------------------------------------------------------- GF(256) 塔分解层 */
/* 表示层：*_array 用塔坐标，*_11b_array 用 Rijndael 0x11B 字节（边界含 phi/psi）。 */
void pqc_kernel_gf256_tower_init(void);
int pqc_kernel_gf256_tower_selftest(void);
uint8_t pqc_kernel_gf256_tower_phi(uint8_t a);
uint8_t pqc_kernel_gf256_tower_psi(uint8_t a);
void pqc_kernel_gf256_tower_phi_array(uint8_t *dst, const uint8_t *src, size_t n);
void pqc_kernel_gf256_tower_psi_array(uint8_t *dst, const uint8_t *src, size_t n);

uint8_t pqc_kernel_gf256_tower_mul(uint8_t a, uint8_t b);
uint8_t pqc_kernel_gf256_tower_mul_11b(uint8_t a, uint8_t b);
uint8_t pqc_kernel_gf256_tower_square(uint8_t a);
uint8_t pqc_kernel_gf256_tower_square_11b(uint8_t a);
uint8_t pqc_kernel_gf256_tower_frobenius(uint8_t a);
uint8_t pqc_kernel_gf256_tower_frobenius_11b(uint8_t a);
uint8_t pqc_kernel_gf256_tower_inv(uint8_t a);
uint8_t pqc_kernel_gf256_tower_inv_11b(uint8_t a);
uint8_t pqc_kernel_gf256_tower_pow(uint8_t a, uint32_t e);
uint8_t pqc_kernel_gf256_tower_pow_11b(uint8_t a, uint32_t e);

/* 标量参考数组。 */
void pqc_kernel_gf256_tower_mul_ref_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_mul_11b_ref_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_mulacc_ref_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_mulacc_11b_ref_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);
uint8_t pqc_kernel_gf256_tower_dot_ref_array(const uint8_t *a, const uint8_t *b, size_t n);
uint8_t pqc_kernel_gf256_tower_dot_11b_ref_array(const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_const_mul_ref_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf256_tower_const_mul_11b_ref_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf256_tower_square_ref_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_square_11b_ref_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_frobenius_ref_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_frobenius_11b_ref_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_inv_ref_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_inv_11b_ref_array(uint8_t *out, const uint8_t *a, size_t n);

/* SIMD 数组（塔坐标）。 */
void pqc_kernel_gf256_tower_mul_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_mulacc_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);
uint8_t pqc_kernel_gf256_tower_dot_array(const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_const_mul_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf256_tower_square_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_frobenius_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_inv_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_pow_array(uint8_t *out, const uint8_t *a, uint32_t e, size_t n);

/* SIMD 数组（0x11B 表示，边界含 phi/psi）。 */
void pqc_kernel_gf256_tower_mul_11b_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_mulacc_11b_array(uint8_t *acc, const uint8_t *a, const uint8_t *b, size_t n);
uint8_t pqc_kernel_gf256_tower_dot_11b_array(const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_const_mul_11b_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
void pqc_kernel_gf256_tower_square_11b_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_frobenius_11b_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_inv_11b_array(uint8_t *out, const uint8_t *a, size_t n);
void pqc_kernel_gf256_tower_pow_11b_array(uint8_t *out, const uint8_t *a, uint32_t e, size_t n);

/* 向量级塔内核：输入输出都是塔坐标，边界转换由调用方负责。
 * 供第三方代码（例如 MQOM 字段层）把塔直接编进内循环。 */
pqc_kernel_gf16_vec pqc_kernel_gf256_tower_phi_block(pqc_kernel_gf16_vec v);
pqc_kernel_gf16_vec pqc_kernel_gf256_tower_psi_block(pqc_kernel_gf16_vec v);
pqc_kernel_gf16_vec pqc_kernel_gf256_tower_mul_block(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b);
void pqc_kernel_gf256_tower_mulacc_block(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b, pqc_kernel_gf16_vec *lo, pqc_kernel_gf16_vec *hi);
pqc_kernel_gf16_vec pqc_kernel_gf256_tower_const_mul_block(pqc_kernel_gf16_vec a, uint8_t k);

/* 仿射基变换版（GF2P8AFFINEQB 代替 shuffle 查表做 φ/ψ，0x11B 进/出）。
 * 无 GFNI 时这些符号回退到对应的 shuffle 版。 */
void pqc_kernel_gf256_tower_mul_11b_affine_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);
void pqc_kernel_gf256_tower_const_mul_11b_affine_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n);
uint8_t pqc_kernel_gf256_tower_dot_11b_affine_array(const uint8_t *a, const uint8_t *b, size_t n);

/* ---------------------------------------------------------------- GF(256) 矩阵级 */
/* 平台 SIMD 字节乘别名：矩阵层用它保持 x86/ARMv8 同一份内核代码。 */
#define GF256_SIMD_MUL pqc_kernel_gf256_vec_mul_avx2
#define GF256_MATRIX_HAS_GFNI 1

/* M1 常数系数收缩：out[r][lane] = XOR_c coeff[r][c] * sig[c][lane]。
 * M2 变量矩阵收缩：out[r][lane] = XOR_c left[(r*cols+c)][lane] * right[c][lane]。
 * 两者都按 3.0 节的收缩口径归一化，独立收缩测试用 l=1、R=C=L=64。 */
void pqc_kernel_gf256_ps_scalar(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                     size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_ps_swar(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                   size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_ps_simd(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                   size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_ps_gfni(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                   size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_sts_scalar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                      size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_sts_swar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                    size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_sts_simd(uint8_t *out, const uint8_t *left, const uint8_t *right,
                    size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_sts_gfni(uint8_t *out, const uint8_t *left, const uint8_t *right,
                    size_t rows, size_t cols, size_t lanes);

/* 塔坐标内核（转换在外面做，scratch 需要 max(rows*cols, cols*lanes, rows*lanes)
 * 或 rows*cols*lanes + cols*lanes + rows*lanes 字节）。 */
void pqc_kernel_gf256_ps_tower(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                    size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_ps_tower_11b(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                        size_t rows, size_t cols, size_t lanes, uint8_t *scratch);
void pqc_kernel_gf256_sts_tower(uint8_t *out, const uint8_t *left, const uint8_t *right,
                     size_t rows, size_t cols, size_t lanes);
void pqc_kernel_gf256_sts_tower_11b(uint8_t *out, const uint8_t *left, const uint8_t *right,
                         size_t rows, size_t cols, size_t lanes, uint8_t *scratch);

#endif
