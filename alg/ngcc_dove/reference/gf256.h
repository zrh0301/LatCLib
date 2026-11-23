#ifndef GF256_H
#define GF256_H

#include <stdint.h>
#include <stddef.h>

// 有限域 GF(2^8)，不可约多项式: x^8 + x^4 + x^3 + x + 1 (0x11b)
//
// 本实现基于 UOV 参考实现 (DOVE/2/ref) 中的基本运算：
//   - 单元素 gf256_mul / gf256_inv（来自 uov/gf16.h，static inline）
//   - 向量 gf256v_add / gf256v_madd / gf256v_mul_scalar /
//     gf256v_set_zero / gf256v_conditional_add（由本目录 gf256.c 实现，
//     内部使用 uov/blas_u32.h 中的 32 位向量原语）
// 与原始 DOVE 实现的不可约多项式一致（均为 0x11b），
// 因此所有计算过程与产生的结果均与原实现保持一致。

// 复用库的 GF(256) 内核，模多项式同为 0x11b
// 提供：gf256_mul / gf256_inv
#include "gf/gf256.h"

// ----------------------- 单元素运算 -----------------------

// 域加法（等价于异或）
uint8_t gf256_add(uint8_t a, uint8_t b);

// 域乘法：GF(2^8) 内的乘法（由 uov/gf16.h 提供 static inline 版本）
// （不重复声明，避免与 uov/gf16.h 的 static inline 冲突）

// 域求逆（由 uov/gf16.h 提供 static inline 版本）

// ----------------------- 向量点积 -----------------------

// 向量点积：sum(a[i] * b[i]) over GF(256)
uint8_t gf256_dot(const uint8_t *a, const uint8_t *b, size_t len);

// ----------------------- 向量基本运算 -----------------------
// 下列函数与 UOV 参考实现 (DOVE/2/ref) 中的 gf256v_* 完全等价。
// （接口命名沿用 UOV 风格：b ^= a；c ^= b * a；等）

// b ^= a ，长度为 _num_byte 字节
void gf256v_add(uint8_t *b, const uint8_t *a, unsigned _num_byte);

// accu_c ^= b * a  （对每个字节逐元素做 accu_c[k] ^= b * a[k]），长度为 _num_byte 字节
void gf256v_madd(uint8_t *accu_c, const uint8_t *a, uint8_t gf256_b, unsigned _num_byte);

// a *= b  （对每个字节逐元素做 a[k] *= b），长度为 _num_byte 字节
void gf256v_mul_scalar(uint8_t *a, uint8_t b, unsigned _num_byte);

// a = 0 ，长度为 _num_byte 字节
void gf256v_set_zero(uint8_t *a, unsigned _num_byte);

// 若 condition 非零则 accu_b ^= a ，否则不变（与 UOV 实现一致，常数时间）
void gf256v_conditional_add(uint8_t *accu_b, uint8_t condition, const uint8_t *a, unsigned _num_byte);

#endif
