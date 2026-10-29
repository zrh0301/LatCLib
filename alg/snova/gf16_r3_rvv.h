/* SPDX-License-Identifier: MIT */
#ifndef SNOVA_R3_GF16_RVV_H
#define SNOVA_R3_GF16_RVV_H

#include <stddef.h>
#include <stdint.h>

#if defined(__riscv_vector)
#include <riscv_vector.h>
#endif

static inline uint16_t snova_r3_gf16_compress_scalar(uint16_t value) {
	uint16_t folded = (uint16_t)((value & 0xf) ^ ((value & 0xf0) >> 3) ^
		((value & 0xf00) >> 6) ^ ((value & 0xf000) >> 9));
	return (uint16_t)((folded ^ ((folded & 0xf0) >> 3) ^ (folded >> 4)) & 0xf);
}

/* Set VL once when a caller is about to issue many AXPY operations of the
 * same length. */
static inline size_t snova_r3_gf16_setvl(size_t count) {
#if defined(__riscv_vector)
	return __riscv_vsetvl_e16m2(count);
#else
	return count;
#endif
}

/* AXPY using a caller-provided VL. Common kernel lengths are exact multiples
 * of VL, so repeated calls avoid a vsetvl instruction. */
static inline void snova_r3_gf16_axpy_vl(uint16_t *dst, const uint16_t *src,
	uint16_t scalar, size_t count, size_t vl) {
#if defined(__riscv_vector)
	while (count != 0) {
		size_t chunk = count < vl ? count : vl;
		if (chunk != vl)
			chunk = __riscv_vsetvl_e16m2(chunk);
		vuint16m2_t vd = __riscv_vle16_v_u16m2(dst, chunk);
		vuint16m2_t vs = __riscv_vle16_v_u16m2(src, chunk);
		vd = __riscv_vxor_vv_u16m2(vd,
			__riscv_vmul_vx_u16m2(vs, scalar, chunk), chunk);
		__riscv_vse16_v_u16m2(dst, vd, chunk);
		dst += chunk;
		src += chunk;
		count -= chunk;
	}
#else
	(void)vl;
	for (size_t i = 0; i < count; ++i)
		dst[i] ^= (uint16_t)(scalar * src[i]);
#endif
}

/* Expanded-domain GF(16) AXPY: dst[i] ^= scalar * src[i]. */
static inline void snova_r3_gf16_axpy(uint16_t *dst, const uint16_t *src,
		uint16_t scalar, size_t count) {
#if defined(__riscv_vector)
	while (count != 0) {
		size_t vl = snova_r3_gf16_setvl(count);
		snova_r3_gf16_axpy_vl(dst, src, scalar, vl, vl);
		dst += vl;
		src += vl;
		count -= vl;
	}
#else
	for (size_t i = 0; i < count; ++i)
		dst[i] ^= (uint16_t)(scalar * src[i]);
#endif
}

/* Accumulate many expanded-domain GF(16) products into one vector.  Keeping
 * the destination in a vector register avoids loading and storing it for
 * every short AXPY term.  Each source points to `count` contiguous elements. */
static inline void snova_r3_gf16_accumulate_terms(uint16_t *dst,
	const uint16_t *const *srcs, const uint16_t *scalars,
	size_t terms, size_t count) {
#if defined(__riscv_vector)
	size_t offset = 0;
	while (count != 0) {
		size_t vl = __riscv_vsetvl_e16m2(count);
		vuint16m2_t acc = __riscv_vle16_v_u16m2(dst + offset, vl);
		for (size_t term = 0; term < terms; ++term) {
			vuint16m2_t src = __riscv_vle16_v_u16m2(srcs[term] + offset, vl);
			acc = __riscv_vxor_vv_u16m2(acc,
				__riscv_vmul_vx_u16m2(src, scalars[term], vl), vl);
		}
		__riscv_vse16_v_u16m2(dst + offset, acc, vl);
		offset += vl;
		count -= vl;
	}
#else
	for (size_t term = 0; term < terms; ++term)
		for (size_t i = 0; i < count; ++i)
			dst[i] ^= (uint16_t)(scalars[term] * srcs[term][i]);
#endif
}

static inline void snova_r3_gf16_compress(uint16_t *row, size_t count) {
#if defined(__riscv_vector)
	while (count != 0) {
		size_t vl = __riscv_vsetvl_e16m2(count);
		vuint16m2_t value = __riscv_vle16_v_u16m2(row, vl);
		vuint16m2_t folded = __riscv_vxor_vv_u16m2(
			__riscv_vand_vx_u16m2(value, 0xf, vl),
			__riscv_vsrl_vx_u16m2(__riscv_vand_vx_u16m2(value, 0xf0, vl), 3, vl), vl);
		folded = __riscv_vxor_vv_u16m2(folded,
			__riscv_vsrl_vx_u16m2(__riscv_vand_vx_u16m2(value, 0xf00, vl), 6, vl), vl);
		folded = __riscv_vxor_vv_u16m2(folded,
			__riscv_vsrl_vx_u16m2(__riscv_vand_vx_u16m2(value, 0xf000, vl), 9, vl), vl);
		vuint16m2_t compressed = __riscv_vxor_vv_u16m2(folded,
			__riscv_vxor_vv_u16m2(
				__riscv_vsrl_vx_u16m2(__riscv_vand_vx_u16m2(folded, 0xf0, vl), 3, vl),
				__riscv_vsrl_vx_u16m2(folded, 4, vl), vl), vl);
		compressed = __riscv_vand_vx_u16m2(compressed, 0xf, vl);
		__riscv_vse16_v_u16m2(row, compressed, vl);
		row += vl;
		count -= vl;
	}
#else
	for (size_t i = 0; i < count; ++i)
		row[i] = snova_r3_gf16_compress_scalar(row[i]);
#endif
}

static inline void snova_r3_gf16_mul_compress(uint16_t *row, uint16_t scalar,
		size_t count) {
#if defined(__riscv_vector)
	while (count != 0) {
		size_t vl = __riscv_vsetvl_e16m2(count);
		vuint16m2_t value = __riscv_vle16_v_u16m2(row, vl);
		value = __riscv_vmul_vx_u16m2(value, scalar, vl);
		__riscv_vse16_v_u16m2(row, value, vl);
		snova_r3_gf16_compress(row, vl);
		row += vl;
		count -= vl;
	}
#else
	for (size_t i = 0; i < count; ++i)
		row[i] = snova_r3_gf16_compress_scalar((uint16_t)(row[i] * scalar));
#endif
}

/* Constant-time row merge used by pivot selection: dst ^= mask & src. */
static inline void snova_r3_xor_mask(uint16_t *dst, const uint16_t *src,
		uint16_t mask, size_t count) {
#if defined(__riscv_vector)
	while (count != 0) {
		size_t vl = __riscv_vsetvl_e16m2(count);
		vuint16m2_t vd = __riscv_vle16_v_u16m2(dst, vl);
		vuint16m2_t vs = __riscv_vle16_v_u16m2(src, vl);
		vs = __riscv_vand_vx_u16m2(vs, mask, vl);
		vd = __riscv_vxor_vv_u16m2(vd, vs, vl);
		__riscv_vse16_v_u16m2(dst, vd, vl);
		dst += vl;
		src += vl;
		count -= vl;
	}
#else
	for (size_t i = 0; i < count; ++i)
		dst[i] ^= (uint16_t)(mask & src[i]);
#endif
}

#endif
