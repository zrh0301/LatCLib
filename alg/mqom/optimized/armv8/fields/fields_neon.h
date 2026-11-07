#ifndef __FIELDS_NEON_H__
#define __FIELDS_NEON_H__

/* ARMv8-A NEON field backend for MQOM's GF(16)->GF(256) parameter. */
#include "fields_common.h"
#include "fields_ref.h"

#include <arm_neon.h>
#include <string.h>

#include "mqom_gf256_neon.h"

/* MOQM_NEON_FAST_CONST = 1（默认，V2.0）：常数路径走 nibble 拆分寄存器查表。
 * MOQM_NEON_FAST_CONST = 0：退回 V1.0 的通用 8 轮位串行实现，用于在 ARM 上 A/B 复现。
 * 编译期加 -DMQOM_NEON_FAST_CONST=0 切换，其余源码不动；两条路径逐字节等价。 */
#ifndef MOQM_NEON_FAST_CONST
#define MOQM_NEON_FAST_CONST 1
#endif

typedef uint8x16_t mqom_neon_vec;

static inline mqom_neon_vec mqom_neon_mul(mqom_neon_vec a, mqom_neon_vec b) {
	const mqom_neon_vec one = vdupq_n_u8(1);
	const mqom_neon_vec reduction = vdupq_n_u8(0x1b);
	mqom_neon_vec result = vdupq_n_u8(0);
	for (unsigned i = 0; i < 8; ++i) {
		const mqom_neon_vec take = vceqq_u8(vandq_u8(b, one), one);
		const mqom_neon_vec high = vshrq_n_u8(a, 7);
		const mqom_neon_vec reduce = vceqq_u8(high, one);
		result = veorq_u8(result, vandq_u8(a, take));
		a = veorq_u8(vshlq_n_u8(a, 1), vandq_u8(reduction, reduce));
		b = vshrq_n_u8(b, 1);
	}
	return result;
}

static inline mqom_neon_vec mqom_neon_load_tail(const uint8_t *p, uint32_t len) {
	uint8_t tmp[16] = { 0 };
	if (len != 0) memcpy(tmp, p, len);
	return vld1q_u8(tmp);
}

static inline uint8_t mqom_neon_reduce(mqom_neon_vec v) {
	uint8_t tmp[16];
	uint8_t result = 0;
	vst1q_u8(tmp, v);
	for (unsigned i = 0; i < 16; ++i) result ^= tmp[i];
	return result;
}

static inline uint8_t gf256_mult_neon(uint8_t x, uint8_t y) {
	return gf256_mult_ref(x, y);
}

static inline void gf256_constant_vect_mult_neon(uint8_t b, const uint8_t *a,
		uint8_t *c, uint32_t len) {
#if MOQM_NEON_FAST_CONST
	/* 常数乘 = 低/高 nibble 两次 16 项寄存器查表 + 异或，共 2 次 vqtbl1q_u8。 */
	mqom_gf256_const_vect_mul_neon(b, a, c, len);
#else
	const mqom_neon_vec vb = vdupq_n_u8(b);
	uint32_t i = 0;
	for (; i + 16u <= len; i += 16u)
		vst1q_u8(c + i, mqom_neon_mul(vld1q_u8(a + i), vb));
	if (i < len) {
		uint8_t tmp[16];
		vst1q_u8(tmp, mqom_neon_mul(mqom_neon_load_tail(a + i, len - i), vb));
		memcpy(c + i, tmp, len - i);
	}
#endif
}

static inline uint8_t gf256_vect_mult_neon(const uint8_t *a, const uint8_t *b,
		uint32_t len) {
	mqom_neon_vec acc0 = vdupq_n_u8(0);
	mqom_neon_vec acc1 = vdupq_n_u8(0);
	mqom_neon_vec acc2 = vdupq_n_u8(0);
	mqom_neon_vec acc3 = vdupq_n_u8(0);
	uint32_t i = 0;
	if (len < 64u) {
		for (; i + 16u <= len; i += 16u)
			acc0 = veorq_u8(acc0, mqom_neon_mul(vld1q_u8(a + i), vld1q_u8(b + i)));
		if (i < len) {
			const uint32_t tail = len - i;
			acc0 = veorq_u8(acc0, mqom_neon_mul(mqom_neon_load_tail(a + i, tail),
				mqom_neon_load_tail(b + i, tail)));
		}
		return mqom_neon_reduce(acc0);
	}
	for (; i + 64u <= len; i += 64u) {
		acc0 = veorq_u8(acc0, mqom_neon_mul(vld1q_u8(a + i), vld1q_u8(b + i)));
		acc1 = veorq_u8(acc1, mqom_neon_mul(vld1q_u8(a + i + 16), vld1q_u8(b + i + 16)));
		acc2 = veorq_u8(acc2, mqom_neon_mul(vld1q_u8(a + i + 32), vld1q_u8(b + i + 32)));
		acc3 = veorq_u8(acc3, mqom_neon_mul(vld1q_u8(a + i + 48), vld1q_u8(b + i + 48)));
	}
	for (; i + 16u <= len; i += 16u)
		acc0 = veorq_u8(acc0, mqom_neon_mul(vld1q_u8(a + i), vld1q_u8(b + i)));
	if (i < len) {
		const uint32_t tail = len - i;
		acc0 = veorq_u8(acc0, mqom_neon_mul(mqom_neon_load_tail(a + i, tail),
			mqom_neon_load_tail(b + i, tail)));
	}
	return mqom_neon_reduce(veorq_u8(veorq_u8(acc0, acc1), veorq_u8(acc2, acc3)));
}

static inline void gf256_vect_mult_multiple_public_neon(uint8_t *const *c,
		const uint8_t *a, const uint8_t *const *b, uint32_t len, uint32_t n) {
	uint32_t i = 0;
	for (; i + 4u <= n; i += 4u) {
		mqom_neon_vec acc0 = vdupq_n_u8(0);
		mqom_neon_vec acc1 = vdupq_n_u8(0);
		mqom_neon_vec acc2 = vdupq_n_u8(0);
		mqom_neon_vec acc3 = vdupq_n_u8(0);
		for (uint32_t j = 0; j + 16u <= len; j += 16u) {
			const mqom_neon_vec va = vld1q_u8(a + j);
			acc0 = veorq_u8(acc0, mqom_neon_mul(va, vld1q_u8(b[i] + j)));
			acc1 = veorq_u8(acc1, mqom_neon_mul(va, vld1q_u8(b[i + 1] + j)));
			acc2 = veorq_u8(acc2, mqom_neon_mul(va, vld1q_u8(b[i + 2] + j)));
			acc3 = veorq_u8(acc3, mqom_neon_mul(va, vld1q_u8(b[i + 3] + j)));
		}
		if (len % 16u != 0) {
			const uint32_t j = len & ~15u;
			const uint32_t tail = len - j;
			const mqom_neon_vec va = mqom_neon_load_tail(a + j, tail);
			acc0 = veorq_u8(acc0, mqom_neon_mul(va, mqom_neon_load_tail(b[i] + j, tail)));
			acc1 = veorq_u8(acc1, mqom_neon_mul(va, mqom_neon_load_tail(b[i + 1] + j, tail)));
			acc2 = veorq_u8(acc2, mqom_neon_mul(va, mqom_neon_load_tail(b[i + 2] + j, tail)));
			acc3 = veorq_u8(acc3, mqom_neon_mul(va, mqom_neon_load_tail(b[i + 3] + j, tail)));
		}
		*c[i] = mqom_neon_reduce(acc0);
		*c[i + 1] = mqom_neon_reduce(acc1);
		*c[i + 2] = mqom_neon_reduce(acc2);
		*c[i + 3] = mqom_neon_reduce(acc3);
	}
	for (; i < n; ++i) *c[i] = gf256_vect_mult_neon(a, b[i], len);
}

static inline void gf256_mat_mult_neon(const uint8_t *A, const uint8_t *X,
		uint8_t *Y, uint32_t n, matrix_type mtype) {
	if (mtype == TRI_INF) {
		uint32_t i = 0;
		for (; i + 4u <= n; i += 4u) {
			const uint8_t *rows[4] = { A + n * i, A + n * (i + 1),
				A + n * (i + 2), A + n * (i + 3) };
			uint8_t *outs[4] = { Y + i, Y + i + 1, Y + i + 2, Y + i + 3 };
			gf256_vect_mult_multiple_public_neon(outs, X, rows, i + 4u, 4u);
		}
		for (; i < n; ++i) Y[i] = gf256_vect_mult_neon(A + n * i, X, i + 1u);
		return;
	}
	GF256_MAT_MULT(A, X, Y, n, mtype, gf256_vect_mult_neon);
}

static inline uint8_t gf16_gf256_mult_neon(uint8_t a_gf16, uint8_t b_gf256) {
	return gf256_mult_neon((uint8_t)((a_gf16 & 1u) ^
			(-((a_gf16 >> 1) & 1u) & 0xe0u) ^
			(-((a_gf16 >> 2) & 1u) & 0x5du) ^
			(-((a_gf16 >> 3) & 1u) & 0xb0u)), b_gf256);
}

/* Lift packed GF(16) elements to their GF(256) embedding.  The input stores
 * two GF(16) elements per byte, low nibble first, while the output stores one
 * GF(256) element per byte.  Keep the conversion branch-free with respect to
 * the field value because x is derived from the secret key seed. */
static inline void gf256_vect_lift_from_gf16_neon(const uint8_t *a_gf16,
		uint8_t *b_gf256, uint32_t n) {
#if MOQM_NEON_FAST_CONST
	/* 一次 vqtbl1q_u8 处理 16 个元素（原实现是逐元素标量查表）。 */
	mqom_gf16_lift_vect_neon(a_gf16, b_gf256, n);
#else
	for (uint32_t i = 0; i < n; ++i) {
		const uint8_t x = (uint8_t)((a_gf16[i / 2u] >> (4u * (i & 1u))) & 0x0fu);
		b_gf256[i] = gf16_gf256_mult_neon(x, 1u);
	}
#endif
}

static inline mqom_neon_vec mqom_gf16_lift_block(const uint8_t *packed, uint32_t count) {
#if MOQM_NEON_FAST_CONST
	/* 原实现每次调用都要跑 16 次标量乘重建表 + memcpy + 逐字节查表，这里只要 1 次 shuffle。
	 * 注意：count 之后的 lane 必须清零（原实现靠 `= {0}` 初始化保证）。 */
	return mqom_gf16_lift_block_neon(packed, count);
#else
	uint8_t src[8] = { 0 };
	uint8_t lifted[16] = { 0 };
	uint8_t table[16];
	for (unsigned x = 0; x < 16; ++x) table[x] = gf16_gf256_mult_neon((uint8_t)x, 1u);
	memcpy(src, packed, (count + 1u) / 2u);
	for (uint32_t i = 0; i < count; ++i)
		lifted[i] = table[(src[i / 2u] >> (4u * (i & 1u))) & 0x0fu];
	return vld1q_u8(lifted);
#endif
}

static inline uint8_t gf16_gf256_vect_mult_neon(const uint8_t *a_gf16,
		const uint8_t *b_gf256, uint32_t len) {
	mqom_neon_vec acc0 = vdupq_n_u8(0);
	mqom_neon_vec acc1 = vdupq_n_u8(0);
	mqom_neon_vec acc2 = vdupq_n_u8(0);
	mqom_neon_vec acc3 = vdupq_n_u8(0);
	uint32_t i = 0;
	for (; i + 64u <= len; i += 64u) {
		acc0 = veorq_u8(acc0, mqom_neon_mul(mqom_gf16_lift_block(a_gf16 + i / 2u, 16u), vld1q_u8(b_gf256 + i)));
		acc1 = veorq_u8(acc1, mqom_neon_mul(mqom_gf16_lift_block(a_gf16 + i / 2u + 8u, 16u), vld1q_u8(b_gf256 + i + 16)));
		acc2 = veorq_u8(acc2, mqom_neon_mul(mqom_gf16_lift_block(a_gf16 + i / 2u + 16u, 16u), vld1q_u8(b_gf256 + i + 32)));
		acc3 = veorq_u8(acc3, mqom_neon_mul(mqom_gf16_lift_block(a_gf16 + i / 2u + 24u, 16u), vld1q_u8(b_gf256 + i + 48)));
	}
	for (; i + 16u <= len; i += 16u)
		acc0 = veorq_u8(acc0, mqom_neon_mul(mqom_gf16_lift_block(a_gf16 + i / 2u, 16u), vld1q_u8(b_gf256 + i)));
	if (i < len) {
		const uint32_t tail = len - i;
		acc0 = veorq_u8(acc0, mqom_neon_mul(mqom_gf16_lift_block(a_gf16 + i / 2u, tail), mqom_neon_load_tail(b_gf256 + i, tail)));
	}
	return mqom_neon_reduce(veorq_u8(veorq_u8(acc0, acc1), veorq_u8(acc2, acc3)));
}

static inline void gf16_gf256_constant_vect_mult_neon(uint8_t a_gf16,
		const uint8_t *b, uint8_t *c, uint32_t n) {
#if MOQM_NEON_FAST_CONST
	/* GF(16) 常数先嵌入 GF(256)，再走 nibble 拆分常数乘。 */
	mqom_gf16_gf256_const_vect_mul_neon(a_gf16, b, c, n);
#else
	gf256_constant_vect_mult_neon(gf16_gf256_mult_neon(a_gf16, 1u), b, c, n);
#endif
}

static inline void gf256_gf16_constant_vect_mult_neon(uint8_t a_gf256,
		const uint8_t *b, uint8_t *c, uint32_t n) {
#if MOQM_NEON_FAST_CONST
	/* 常数乘烘进提升表 ⇒ 整个函数只剩 1 次 vqtbl1q_u8 / 16 个输出。 */
	mqom_gf256_gf16_const_vect_mul_neon(a_gf256, b, c, n);
#else
	for (uint32_t i = 0; i < n; ++i)
		c[i] = gf16_gf256_mult_neon((uint8_t)((b[i / 2u] >> (4u * (i & 1u))) & 0x0fu), a_gf256);
#endif
}

static inline uint8_t gf256_gf16_vect_mult_neon(const uint8_t *a,
		const uint8_t *b, uint32_t n) {
	return gf16_gf256_vect_mult_neon(b, a, n);
}

static inline void gf256_gf16_mat_mult_neon(const uint8_t *A, const uint8_t *X,
		uint8_t *Y, uint32_t n, matrix_type mtype) {
	GF256_GF16_MAT_MULT(A, X, Y, n, mtype, gf256_gf16_vect_mult_neon);
}

static inline void gf256_mat_transpose_neon(const uint8_t *A, uint8_t *B,
		uint32_t n, matrix_type mtype) {
	GF256_MAT_TRANSPOSE(A, B, n, mtype);
}

#endif /* __FIELDS_NEON_H__ */
