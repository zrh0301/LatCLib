#include <string.h>

#if defined(VALGRIND)
#include <valgrind/memcheck.h>
#endif

#include "api.h"
#include "rng.h"
#include "snova_internal.h"

#ifndef SNOVA_VERIFY_VTL
#define SNOVA_VERIFY_VTL 0
#endif

#ifndef SNOVA_VERIFY_EXPANDED_NEON
#define SNOVA_VERIFY_EXPANDED_NEON 1
#endif

#ifndef SNOVA_VERIFY_EXPANDED_ILP4
#define SNOVA_VERIFY_EXPANDED_ILP4 0
#endif

#ifndef SNOVA_VERIFY_DIRECT_MTK2
#define SNOVA_VERIFY_DIRECT_MTK2 1
#endif

#if SNOVA_VERIFY_DIRECT_MTK2 && SNOVA_VERIFY_VTL
#error "SNOVA_VERIFY_DIRECT_MTK2 and SNOVA_VERIFY_VTL are mutually exclusive"
#endif

#define mvl_SNOVA32 ((m_SNOVA * v_SNOVA * l_SNOVA + VEC_LENGTH - 1) / VEC_LENGTH)
#define mvl_SNOVA8 (mvl_SNOVA32 * VEC_LENGTH)
#define mol_SNOVA32 ((m_SNOVA * o_SNOVA * l_SNOVA + VEC_LENGTH - 1) / VEC_LENGTH)
#define mol_SNOVA8 (mol_SNOVA32 * VEC_LENGTH)
#define vl_SNOVA (v_SNOVA * l_SNOVA)
#define ol_SNOVA (o_SNOVA * l_SNOVA)
#define ol_SNOVA_floor2 ((ol_SNOVA) / 2 * 2)
#define ol_SNOVA_pair ((ol_SNOVA) / 2)
#ifndef GEN_F_PAIRS
#define GEN_F_PAIRS 2
#endif
#ifndef GEN_F_CHUNK
#define GEN_F_CHUNK 6
#endif
#if ((ol_SNOVA_pair % GEN_F_PAIRS) == 0) && ((mvl_SNOVA32 % GEN_F_CHUNK) == 0)
#define GEN_F_TILED 1
#else
#define GEN_F_TILED 0
#endif
#if SNOVA_NEON_RANK4
static inline uint8x16_t gather_rank_rows(const uint8_t* p, int row_stride) {
	uint32_t w0, w1, w2, w3;
	memcpy(&w0, p + 0 * row_stride, 4);
	memcpy(&w1, p + 1 * row_stride, 4);
	memcpy(&w2, p + 2 * row_stride, 4);
	memcpy(&w3, p + 3 * row_stride, 4);
	uint32x4_t v = vdupq_n_u32(w0);
	v = vsetq_lane_u32(w1, v, 1);
	v = vsetq_lane_u32(w2, v, 2);
	v = vsetq_lane_u32(w3, v, 3);
	return vreinterpretq_u8_u32(v);
}
static inline void xor_block_from_rows(uint8_t* block, const uint8_t* p, int row_stride) {
	uint8x16_t acc = vld1q_u8(block);
	vst1q_u8(block, veorq_u8(acc, gather_rank_rows(p, row_stride)));
}
static inline void xor_block_from_rows_tr(uint8_t* block, const uint8_t* p, int row_stride) {
	const uint8x16_t idx = vld1q_u8(snova_tr_block_idx);
	uint8x16_t acc = vld1q_u8(block);
	vst1q_u8(block, veorq_u8(acc, vqtbl1q_u8(gather_rank_rows(p, row_stride), idx)));
}
#elif SNOVA_AVX2_RANK4
static inline __m128i gather_rank_rows_avx2(const uint8_t* p, int row_stride) {
	uint32_t w0, w1, w2, w3;
	memcpy(&w0, p + 0 * row_stride, sizeof(w0));
	memcpy(&w1, p + 1 * row_stride, sizeof(w1));
	memcpy(&w2, p + 2 * row_stride, sizeof(w2));
	memcpy(&w3, p + 3 * row_stride, sizeof(w3));
	return _mm_set_epi32((int)w3, (int)w2, (int)w1, (int)w0);
}
static inline void xor_block_from_rows(uint8_t* block, const uint8_t* p, int row_stride) {
	__m128i acc = _mm_loadu_si128((const __m128i*)block);
	_mm_storeu_si128((__m128i*)block, _mm_xor_si128(acc, gather_rank_rows_avx2(p, row_stride)));
}
static inline void xor_block_from_rows_tr(uint8_t* block, const uint8_t* p, int row_stride) {
	const __m128i idx = _mm_loadu_si128((const __m128i*)snova_tr_block_idx);
	__m128i acc = _mm_loadu_si128((const __m128i*)block);
	__m128i rows = _mm_shuffle_epi8(gather_rank_rows_avx2(p, row_stride), idx);
	_mm_storeu_si128((__m128i*)block, _mm_xor_si128(acc, rows));
}
#endif
void snova_gen_seeds_and_T12(T12_t T12, const uint8_t* seed) {
	SNOVA_PROFILE_BEGIN(genkeys_t12_t);
	uint8_t prng_output_private[bytes_prng_private];
	gf16_t GF16_prng_output_private[GF16s_prng_private];
	gf16_t* pt_array;
	shake256(prng_output_private, bytes_prng_private, seed, seed_length_private);
	snova_bytes_to_gf16(prng_output_private, GF16_prng_output_private, GF16s_prng_private);
	pt_array = GF16_prng_output_private;
	for (int j = 0; j < v_SNOVA; ++j) {
		for (int k = 0; k < o_SNOVA; ++k) {
			snova_gen_a_FqS_ct(pt_array, T12[j][k]);
			pt_array += rank;
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_T12, genkeys_t12_t);
}
#if FIXED_ABQ
static uint8_t fixed_abq[4 * m_SNOVA * alpha_SNOVA * lsq_SNOVA] = {0};
void snova_gen_fixed_ABQ(const char* abq_seed) {
	uint8_t rng_out[m_SNOVA * alpha_SNOVA * (lsq_SNOVA + l_SNOVA)];
	uint8_t q12[2 * m_SNOVA * alpha_SNOVA * l_SNOVA];
	shake256(rng_out, sizeof(rng_out), (const uint8_t*)abq_seed, strlen(abq_seed));
	snova_bytes_to_gf16(rng_out, fixed_abq, 2 * m_SNOVA * alpha_SNOVA * lsq_SNOVA);
	snova_bytes_to_gf16(&rng_out[m_SNOVA * alpha_SNOVA * lsq_SNOVA], q12, 2 * m_SNOVA * alpha_SNOVA * l_SNOVA);
	for (int pi = 0; pi < m_SNOVA; ++pi) {
		for (int index = 0; index < alpha_SNOVA; ++index) {
			snova_be_invertible_by_add_aS(&fixed_abq[(pi * alpha_SNOVA + index) * lsq_SNOVA]);
		}
		for (int index = 0; index < alpha_SNOVA; ++index) {
			snova_be_invertible_by_add_aS(&fixed_abq[((m_SNOVA + pi) * alpha_SNOVA + index) * lsq_SNOVA]);
		}
		for (int index = 0; index < alpha_SNOVA; ++index) {
			snova_gen_a_FqS(&q12[(pi * alpha_SNOVA + index) * l_SNOVA],
			                 &fixed_abq[((2 * m_SNOVA + pi) * alpha_SNOVA + index) * lsq_SNOVA]);
		}
		for (int index = 0; index < alpha_SNOVA; ++index) {
			snova_gen_a_FqS(&q12[((m_SNOVA + pi) * alpha_SNOVA + index) * l_SNOVA],
			                 &fixed_abq[((3 * m_SNOVA + pi) * alpha_SNOVA + index) * lsq_SNOVA]);
		}
	}
}
#endif
typedef char snova_public_key_expand_fits_expanded_PK[(sizeof(public_key_expand) <= sizeof(expanded_PK)) ? 1 : -1];

#define DIRECT_MTK2_P_BYTES (SNOVA_m1 * SNOVA_n * SNOVA_n * SNOVA_l2)
#define DIRECT_MTK2_AM_BYTES (SNOVA_o * SNOVA_alpha * SNOVA_r2)
#define DIRECT_MTK2_BM_BYTES (SNOVA_o * SNOVA_alpha * SNOVA_lr)
#define DIRECT_MTK2_Q1_BYTES (SNOVA_o * SNOVA_alpha * SNOVA_l)
#define DIRECT_MTK2_Q2_BYTES (SNOVA_o * SNOVA_alpha * SNOVA_l)
#define DIRECT_MTK2_P_OFFSET 0
#define DIRECT_MTK2_AM_OFFSET (DIRECT_MTK2_P_OFFSET + DIRECT_MTK2_P_BYTES)
#define DIRECT_MTK2_BM_OFFSET (DIRECT_MTK2_AM_OFFSET + DIRECT_MTK2_AM_BYTES)
#define DIRECT_MTK2_Q1_OFFSET (DIRECT_MTK2_BM_OFFSET + DIRECT_MTK2_BM_BYTES)
#define DIRECT_MTK2_Q2_OFFSET (DIRECT_MTK2_Q1_OFFSET + DIRECT_MTK2_Q1_BYTES)
#define DIRECT_MTK2_PK_SEED_OFFSET (DIRECT_MTK2_Q2_OFFSET + DIRECT_MTK2_Q2_BYTES)
#define DIRECT_MTK2_TOTAL_BYTES (DIRECT_MTK2_PK_SEED_OFFSET + SEED_LENGTH_PUBLIC)

typedef char snova_direct_mtk2_pk_fits_expanded_PK[
	(DIRECT_MTK2_TOTAL_BYTES <= sizeof(expanded_PK)) ? 1 : -1];

static inline uint16_t official_gf16_expand(gf16_t a) {
	uint16_t val = (uint16_t)(a | (a << 3) | (a << 6) | (a << 9));
	return val & 0x1111;
}

static inline uint16_t official_gf16_compress(uint16_t a) {
	uint16_t val = (uint16_t)((a & 0x000f) ^ ((a & 0x00f0) >> 3) ^ ((a & 0x0f00) >> 6) ^ ((a & 0xf000) >> 9));
	return (uint16_t)((val ^ ((val & 0x00f0) >> 3) ^ (val >> 4)) & 0x000f);
}

static inline uint16_t official_gf16_cleanup(uint16_t val) {
	return official_gf16_expand((gf16_t)official_gf16_compress(val));
}

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
/*
 * The representation stores one GF(16) value in the four bit
 * positions 0, 4, 8 and 12 of a uint16_t.  The verify path can therefore
 * keep eight independent values in one uint16x8_t lane group.  Multiplication
 * by a raw GF(16) coefficient is ordinary uint16 multiplication followed by
 * XOR accumulation; compression is delayed until the same boundaries as the
 * scalar routine.
 */
static inline uint16x8_t official_expand_u8x8_neon(const uint8_t* src) {
	const uint8x8_t raw = vld1_u8(src);
	const uint16x8_t x = vmovl_u8(raw);
	const uint16x8_t expanded = vorrq_u16(
		vorrq_u16(x, vshlq_n_u16(x, 3)),
		vorrq_u16(vshlq_n_u16(x, 6), vshlq_n_u16(x, 9)));
	return vandq_u16(expanded, vdupq_n_u16(0x1111));
}

static inline uint16x8_t official_compress_u16x8_neon(uint16x8_t x) {
	uint16x8_t v = veorq_u16(
		veorq_u16(vandq_u16(x, vdupq_n_u16(0x000f)),
		          vshrq_n_u16(vandq_u16(x, vdupq_n_u16(0x00f0)), 3)),
		veorq_u16(vshrq_n_u16(vandq_u16(x, vdupq_n_u16(0x0f00)), 6),
		          vshrq_n_u16(vandq_u16(x, vdupq_n_u16(0xf000)), 9)));
	v = veorq_u16(v, vshrq_n_u16(vandq_u16(v, vdupq_n_u16(0x00f0)), 3));
	v = veorq_u16(v, vshrq_n_u16(v, 4));
	return vandq_u16(v, vdupq_n_u16(0x000f));
}

static inline void official_expand_array_neon(uint16_t* dst, const uint8_t* src, size_t count) {
	size_t i = 0;
	for (; i + 8 <= count; i += 8) {
		vst1q_u16(dst + i, official_expand_u8x8_neon(src + i));
	}
	for (; i < count; ++i) {
		dst[i] = official_gf16_expand(src[i]);
	}
}

static inline void official_compress_array_neon(uint16_t* dst, size_t count) {
	size_t i = 0;
	for (; i + 8 <= count; i += 8) {
		vst1q_u16(dst + i, official_compress_u16x8_neon(vld1q_u16(dst + i)));
	}
	for (; i < count; ++i) {
		dst[i] = official_gf16_compress(dst[i]);
	}
}

static inline void official_raw_expanded_xor_array_neon(uint16_t* acc, const uint8_t* raw,
	                                                       const uint16_t* expanded, size_t count) {
	size_t i = 0;
	for (; i + 8 <= count; i += 8) {
		const uint16x8_t av = vld1q_u16(acc + i);
		const uint16x8_t bv = vld1q_u16(expanded + i);
		const uint16x8_t coeff = vmovl_u8(vld1_u8(raw + i));
		vst1q_u16(acc + i, veorq_u16(av, vmulq_u16(coeff, bv)));
	}
	for (; i < count; ++i) {
		acc[i] ^= (uint16_t)raw[i] * expanded[i];
	}
}

static inline void official_const_expanded_xor_array_neon(uint16_t* acc, uint16_t coeff,
	                                                         const uint16_t* expanded, size_t count) {
	size_t i = 0;
	/* The hot path always processes 16 lanes. Keep two independent vectors
	 * visible to the compiler so the multiply and load streams can overlap. */
	for (; i + 16 <= count; i += 16) {
		const uint16x8_t a0 = vld1q_u16(acc + i);
		const uint16x8_t a1 = vld1q_u16(acc + i + 8);
		const uint16x8_t b0 = vld1q_u16(expanded + i);
		const uint16x8_t b1 = vld1q_u16(expanded + i + 8);
		vst1q_u16(acc + i, veorq_u16(a0, vmulq_n_u16(b0, coeff)));
		vst1q_u16(acc + i + 8, veorq_u16(a1, vmulq_n_u16(b1, coeff)));
	}
	for (; i + 8 <= count; i += 8) {
		vst1q_u16(acc + i, veorq_u16(vld1q_u16(acc + i),
		                              vmulq_n_u16(vld1q_u16(expanded + i), coeff)));
	}
	for (; i < count; ++i) {
		acc[i] ^= coeff * expanded[i];
	}
}

static inline void official_const4_expanded_xor_array_neon(
	uint16_t* restrict acc0, uint16_t* restrict acc1,
	uint16_t* restrict acc2, uint16_t* restrict acc3,
	uint16_t coeff0, uint16_t coeff1, uint16_t coeff2, uint16_t coeff3,
	const uint16_t* restrict expanded, size_t count) {
	size_t i = 0;
	for (; i + 8 <= count; i += 8) {
		const uint16x8_t b = vld1q_u16(expanded + i);
		const uint16x8_t a0 = vld1q_u16(acc0 + i);
		const uint16x8_t a1 = vld1q_u16(acc1 + i);
		const uint16x8_t a2 = vld1q_u16(acc2 + i);
		const uint16x8_t a3 = vld1q_u16(acc3 + i);
		vst1q_u16(acc0 + i, veorq_u16(a0, vmulq_n_u16(b, coeff0)));
		vst1q_u16(acc1 + i, veorq_u16(a1, vmulq_n_u16(b, coeff1)));
		vst1q_u16(acc2 + i, veorq_u16(a2, vmulq_n_u16(b, coeff2)));
		vst1q_u16(acc3 + i, veorq_u16(a3, vmulq_n_u16(b, coeff3)));
	}
	for (; i < count; ++i) {
		const uint16_t b = expanded[i];
		acc0[i] ^= coeff0 * b;
		acc1[i] ^= coeff1 * b;
		acc2[i] ^= coeff2 * b;
		acc3[i] ^= coeff3 * b;
	}
}
#endif

#if !FIXED_ABQ
static inline gf16_t snova_adjust_q_coeff(const gf16_t* coeff, int idx) {
	if (idx != l_SNOVA - 1 || coeff[idx] != 0) {
		return coeff[idx];
	}
	return (gf16_t)(SNOVA_q - (coeff[0] + (coeff[0] == 0)));
}

static inline void snova_store_q_coeffs(gf16_t* restrict dst, const gf16_t* restrict src) {
	for (int idx = 0; idx < l_SNOVA; ++idx) {
		dst[idx] = snova_adjust_q_coeff(src, idx);
	}
}
#endif

static void snova_gen_ABQ_ex(map_group1* map, const uint8_t* pt_public_key_seed, Qalpha_coeff_t q1_out,
                             Qalpha_coeff_t q2_out) {
	uint8_t prng_output_public[bytes_prng_public];
	uint8_t Q_temp[(sizeof(Qalpha1_t) + sizeof(Qalpha2_t)) / l_SNOVA];
	(void)q1_out;
	(void)q2_out;
	snova_pk_expand_seed(pt_public_key_seed, prng_output_public);
#if FIXED_ABQ
	snova_bytes_to_gf16(prng_output_public, (uint8_t*)map, GF16s_prng_public - sizeof(Q_temp));
	memcpy(map->Aalpha, fixed_abq, 4 * m_SNOVA * alpha_SNOVA * lsq_SNOVA);
#else
	snova_bytes_to_gf16(prng_output_public, (uint8_t*)map, GF16s_prng_public - sizeof(Q_temp));
	snova_bytes_to_gf16(prng_output_public + sizeof(prng_output_public) - ((sizeof(Q_temp) + 1) >> 1), Q_temp,
	                     sizeof(Q_temp));
	for (int pi = 0; pi < m_SNOVA; ++pi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			snova_be_invertible_by_add_aS(map->Aalpha[pi][alpha]);
		}
	}
	for (int pi = 0; pi < m_SNOVA; ++pi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			snova_be_invertible_by_add_aS(map->Balpha[pi][alpha]);
		}
	}
	gf16_t* pt_array = Q_temp;
	for (int pi = 0; pi < m_SNOVA; ++pi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			if (q1_out != NULL) {
				snova_store_q_coeffs(q1_out[pi][alpha], pt_array);
			}
			snova_gen_a_FqS(pt_array, map->Qalpha1[pi][alpha]);
			pt_array += l_SNOVA;
		}
	}
	for (int pi = 0; pi < m_SNOVA; ++pi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			if (q2_out != NULL) {
				snova_store_q_coeffs(q2_out[pi][alpha], pt_array);
			}
			snova_gen_a_FqS(pt_array, map->Qalpha2[pi][alpha]);
			pt_array += l_SNOVA;
		}
	}
#endif
}

void snova_gen_ABQ(map_group1* map, const uint8_t* pt_public_key_seed) {
	snova_gen_ABQ_ex(map, pt_public_key_seed, NULL, NULL);
}

static void snova_gen_public_P_only(map_group1* map, const uint8_t* public_seed) {
	uint8_t prng_output_public[bytes_prng_public];
	const size_t p_elements = sizeof(P11_t) + sizeof(P12_t) + sizeof(P21_t);
	snova_pk_expand_seed(public_seed, prng_output_public);
	snova_bytes_to_gf16(prng_output_public, (uint8_t*)map, (int)p_elements);
}

static void gen_F_contract(VECTOR* restrict res256, const VECTOR* restrict in256,
                           const uint8_t* restrict pair_idx) {
#if GEN_F_TILED
	for (int p_base = 0; p_base < ol_SNOVA_pair; p_base += GEN_F_PAIRS) {
		for (int c_base = 0; c_base < mvl_SNOVA32; c_base += GEN_F_CHUNK) {
			VECTOR acc[GEN_F_PAIRS][GEN_F_CHUNK];
			for (int p = 0; p < GEN_F_PAIRS; ++p) {
				for (int c = 0; c < GEN_F_CHUNK; ++c) {
					acc[p][c] = zero256;
				}
			}
			for (int dk_k1 = 0; dk_k1 < vl_SNOVA; ++dk_k1) {
				VECTOR ktab[GEN_F_PAIRS];
				const VECTOR* in_row = in256 + (size_t)dk_k1 * mvl_SNOVA32 + c_base;
				for (int p = 0; p < GEN_F_PAIRS; ++p) {
					ktab[p] = mtk2_16[pair_idx[dk_k1 * ol_SNOVA_pair + p_base + p]];
				}
				for (int c = 0; c < GEN_F_CHUNK; ++c) {
					VECTOR in_vec = in_row[c];
					for (int p = 0; p < GEN_F_PAIRS; ++p) {
						acc[p][c] ^= VEC_SHUFFLE(ktab[p], in_vec);
					}
				}
			}
			for (int p = 0; p < GEN_F_PAIRS; ++p) {
				const int dj_j1 = 2 * (p_base + p);
				VECTOR* lo = res256 + (size_t)dj_j1 * mvl_SNOVA32 + c_base;
				VECTOR* hi = res256 + (size_t)(dj_j1 + 1) * mvl_SNOVA32 + c_base;
				for (int c = 0; c < GEN_F_CHUNK; ++c) {
					lo[c] = acc[p][c] & l_mask;
					hi[c] = (acc[p][c] >> 4) & l_mask;
				}
			}
		}
	}
#else
	for (int p = 0; p < ol_SNOVA_pair; ++p) {
		const int dj_j1 = 2 * p;
		for (int c = 0; c < mvl_SNOVA32; ++c) {
			res256[dj_j1 * mvl_SNOVA32 + c] = zero256;
			res256[(dj_j1 + 1) * mvl_SNOVA32 + c] = zero256;
		}
		for (int dk_k1 = 0; dk_k1 < vl_SNOVA; ++dk_k1) {
			VECTOR ktab = mtk2_16[pair_idx[dk_k1 * ol_SNOVA_pair + p]];
			for (int c = 0; c < mvl_SNOVA32; ++c) {
				VECTOR prod = VEC_SHUFFLE(ktab, in256[dk_k1 * mvl_SNOVA32 + c]);
				res256[dj_j1 * mvl_SNOVA32 + c] ^= prod & l_mask;
				res256[(dj_j1 + 1) * mvl_SNOVA32 + c] ^= (prod >> 4) & l_mask;
			}
		}
	}
#endif
}
#if ol_SNOVA % 2
static void gen_F_contract_tail(VECTOR* restrict res256, const VECTOR* restrict in256,
                                const uint8_t* restrict single_idx) {
	VECTOR* out = res256 + (size_t)ol_SNOVA_floor2 * mvl_SNOVA32;
	for (int c = 0; c < mvl_SNOVA32; ++c) {
		out[c] = zero256;
	}
	for (int dk_k1 = 0; dk_k1 < vl_SNOVA; ++dk_k1) {
		VECTOR ktab = SNOVA_MULTTAB(single_idx[dk_k1]);
		for (int c = 0; c < mvl_SNOVA32; ++c) {
			out[c] ^= VEC_SHUFFLE(ktab, in256[dk_k1 * mvl_SNOVA32 + c]);
		}
	}
}
#endif
static void gen_F_pack_pair_idx(uint8_t* restrict pair_idx, const uint8_t* restrict t12_8) {
	for (int dk_k1 = 0; dk_k1 < vl_SNOVA; ++dk_k1) {
		for (int p = 0; p < ol_SNOVA_pair; ++p) {
			pair_idx[dk_k1 * ol_SNOVA_pair + p] =
			    (uint8_t)(t12_8[dk_k1 * ol_SNOVA + 2 * p] | (t12_8[dk_k1 * ol_SNOVA + 2 * p + 1] << 4));
		}
	}
}
void snova_gen_F_vtl(map_group2* map2, map_group1* map1, T12_t T12) {
	SNOVA_PROFILE_BEGIN(gen_F_total_t);
	VECTOR p11_256[mvl_SNOVA32 * vl_SNOVA] = {0};
	VECTOR res256[mvl_SNOVA32 * ol_SNOVA] = {0};
	alignas(VEC_LENGTH) uint8_t t12_8[vl_SNOVA * ol_SNOVA] = {0};
	alignas(VEC_LENGTH) uint8_t pair_idx[vl_SNOVA * ol_SNOVA_pair] = {0};
	uint8_t* p11_8 = (uint8_t*)p11_256;
	uint8_t* res8 = (uint8_t*)res256;
	memcpy(map2->F11, map1->P11, sizeof(F11_t));
	memcpy(map2->F12, map1->P12, sizeof(F12_t));
	memcpy(map2->F21, map1->P21, sizeof(F21_t));
	SNOVA_PROFILE_BEGIN(gen_F_f12_pack_t);
	for (int dk = 0; dk < v_SNOVA; ++dk) {
		snova_pack_block_run(t12_8 + (dk * rank) * ol_SNOVA, ol_SNOVA, (const uint8_t*)T12[dk], o_SNOVA);
	}
	gen_F_pack_pair_idx(pair_idx, t12_8);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int dk = 0; dk < v_SNOVA; ++dk) {
			snova_pack_block_run_tr_ex(p11_8 + (dk * rank) * mvl_SNOVA8 + mi * vl_SNOVA, mvl_SNOVA8,
			                            (const uint8_t*)map1->P11[mi][0][dk], v_SNOVA * sq_rank, v_SNOVA);
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_F_F12_PACK, gen_F_f12_pack_t);
	SNOVA_PROFILE_BEGIN(gen_F_f12_mul_t);
	gen_F_contract(res256, p11_256, pair_idx);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_F_F12_MUL, gen_F_f12_mul_t);
#if ol_SNOVA % 2
	{
		alignas(VEC_LENGTH) uint8_t single_idx[vl_SNOVA];
		for (int dk_k1 = 0; dk_k1 < vl_SNOVA; ++dk_k1) {
			single_idx[dk_k1] = t12_8[dk_k1 * ol_SNOVA + ol_SNOVA_floor2];
		}
		gen_F_contract_tail(res256, p11_256, single_idx);
	}
#endif
	SNOVA_PROFILE_BEGIN(gen_F_f12_scatter_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int dj = 0; dj < o_SNOVA; ++dj) {
			for (int di = 0; di < v_SNOVA; ++di) {
				const uint8_t* src = res8 + (dj * rank) * mvl_SNOVA8 + mi * vl_SNOVA + di * rank;
#if SNOVA_RANK4_FAST
				xor_block_from_rows_tr(map2->F12[mi][di][dj], src, mvl_SNOVA8);
#else
				for (int i1 = 0; i1 < rank; ++i1) {
					for (int j1 = 0; j1 < rank; ++j1) {
						map2->F12[mi][di][dj][i1 * rank + j1] ^= src[j1 * mvl_SNOVA8 + i1];
					}
				}
#endif
			}
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_F_F12_SCATTER, gen_F_f12_scatter_t);
	SNOVA_PROFILE_BEGIN(gen_F_f21_pack_t);
	for (int dj = 0; dj < v_SNOVA; ++dj) {
		snova_pack_block_run_tr(t12_8 + (dj * rank) * ol_SNOVA, ol_SNOVA, (const uint8_t*)T12[dj], o_SNOVA);
	}
	gen_F_pack_pair_idx(pair_idx, t12_8);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int di = 0; di < v_SNOVA; ++di) {
			snova_pack_block_run(p11_8 + (di * rank) * mvl_SNOVA8 + mi * vl_SNOVA, mvl_SNOVA8,
			                      (const uint8_t*)map1->P11[mi][di], v_SNOVA);
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_F_F21_PACK, gen_F_f21_pack_t);
	SNOVA_PROFILE_BEGIN(gen_F_f21_mul_t);
	gen_F_contract(res256, p11_256, pair_idx);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_F_F21_MUL, gen_F_f21_mul_t);
#if ol_SNOVA % 2
	{
		alignas(VEC_LENGTH) uint8_t single_idx[vl_SNOVA];
		for (int dk_k1 = 0; dk_k1 < vl_SNOVA; ++dk_k1) {
			single_idx[dk_k1] = t12_8[dk_k1 * ol_SNOVA + ol_SNOVA_floor2];
		}
		gen_F_contract_tail(res256, p11_256, single_idx);
	}
#endif
	SNOVA_PROFILE_BEGIN(gen_F_f21_scatter_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int dj = 0; dj < o_SNOVA; ++dj) {
			for (int di = 0; di < v_SNOVA; ++di) {
				const uint8_t* src = res8 + (dj * rank) * mvl_SNOVA8 + mi * vl_SNOVA + di * rank;
#if SNOVA_RANK4_FAST
				xor_block_from_rows(map2->F21[mi][dj][di], src, mvl_SNOVA8);
#else
				for (int i1 = 0; i1 < rank; ++i1) {
					for (int j1 = 0; j1 < rank; ++j1) {
						map2->F21[mi][dj][di][i1 * rank + j1] ^= src[i1 * mvl_SNOVA8 + j1];
					}
				}
#endif
			}
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_F_F21_SCATTER, gen_F_f21_scatter_t);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_GEN_F, gen_F_total_t);
}
static void gen_P22_contract(uint8_t* restrict temp8, const VECTOR* restrict in256, T12_t T12, int transpose_t12) {
	VECTOR* temp256 = (VECTOR*)temp8;
	for (int dout = 0; dout < o_SNOVA; ++dout) {
		for (int r0 = 0; r0 < rank_floor2; r0 += 2) {
			VECTOR acc[mol_SNOVA32];
			for (int c = 0; c < mol_SNOVA32; ++c) {
				acc[c] = zero256;
			}
			for (int di = 0; di < v_SNOVA; ++di) {
				const gf16_t* t = T12[di][dout];
				for (int k = 0; k < rank; ++k) {
					const uint8_t lo = transpose_t12 ? t[k * rank + r0] : t[r0 * rank + k];
					const uint8_t hi = transpose_t12 ? t[k * rank + r0 + 1] : t[(r0 + 1) * rank + k];
					VECTOR ktab = mtk2_16[lo | (hi << 4)];
					const VECTOR* in_row = in256 + (size_t)(di * rank + k) * mol_SNOVA32;
					for (int c = 0; c < mol_SNOVA32; ++c) {
						acc[c] ^= VEC_SHUFFLE(ktab, in_row[c]);
					}
				}
			}
			VECTOR* out_lo = temp256 + (size_t)(dout * rank + r0) * mol_SNOVA32;
			VECTOR* out_hi = temp256 + (size_t)(dout * rank + r0 + 1) * mol_SNOVA32;
			for (int c = 0; c < mol_SNOVA32; ++c) {
				out_lo[c] = acc[c] & l_mask;
				out_hi[c] = (acc[c] >> 4) & l_mask;
			}
		}
#if rank % 2
		{
			VECTOR acc[mol_SNOVA32];
			for (int c = 0; c < mol_SNOVA32; ++c) {
				acc[c] = zero256;
			}
			for (int di = 0; di < v_SNOVA; ++di) {
				const gf16_t* t = T12[di][dout];
				for (int k = 0; k < rank; ++k) {
					const uint8_t lo = transpose_t12 ? t[k * rank + rank_floor2] : t[rank_floor2 * rank + k];
					VECTOR ktab = SNOVA_MULTTAB(lo);
					const VECTOR* in_row = in256 + (size_t)(di * rank + k) * mol_SNOVA32;
					for (int c = 0; c < mol_SNOVA32; ++c) {
						acc[c] ^= VEC_SHUFFLE(ktab, in_row[c]);
					}
				}
			}
			VECTOR* out = temp256 + (size_t)(dout * rank + rank_floor2) * mol_SNOVA32;
			for (int c = 0; c < mol_SNOVA32; ++c) {
				out[c] = acc[c];
			}
		}
#endif
	}
}
void snova_gen_P22_vtl(P22_byte_t outP22, T12_t T12, P21_t P21, F12_t F12) {
	SNOVA_PROFILE_BEGIN(gen_P22_total_t);
	P22_t P22 = {0};
	alignas(VEC_LENGTH) uint8_t temp_8[mol_SNOVA8 * o_SNOVA * l_SNOVA] = {0};
	alignas(VEC_LENGTH) uint8_t in_8[mol_SNOVA8 * v_SNOVA * l_SNOVA] = {0};
	VECTOR* in_256 = (VECTOR*)in_8;
	SNOVA_PROFILE_BEGIN(gen_P22_f12_pack_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int di = 0; di < v_SNOVA; ++di) {
			snova_pack_block_run(in_8 + (di * rank) * mol_SNOVA8 + mi * (o_SNOVA * rank), mol_SNOVA8,
			                      (const uint8_t*)F12[mi][di], o_SNOVA);
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_F12_PACK, gen_P22_f12_pack_t);
	SNOVA_PROFILE_BEGIN(gen_P22_f12_mul_t);
	gen_P22_contract(temp_8, in_256, T12, 0);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_F12_MUL, gen_P22_f12_mul_t);
	SNOVA_PROFILE_BEGIN(gen_P22_f12_fold_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int dj = 0; dj < o_SNOVA; ++dj) {
			for (int dk = 0; dk < o_SNOVA; ++dk) {
				const uint8_t* src = temp_8 + (dj * rank) * mol_SNOVA8 + mi * (o_SNOVA * rank) + dk * rank;
#if SNOVA_RANK4_FAST
				xor_block_from_rows(P22[mi][dj][dk], src, mol_SNOVA8);
#else
				for (int i1 = 0; i1 < rank; ++i1) {
					for (int j1 = 0; j1 < rank; ++j1) {
						P22[mi][dj][dk][i1 * rank + j1] ^= src[i1 * mol_SNOVA8 + j1];
					}
				}
#endif
			}
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_F12_FOLD, gen_P22_f12_fold_t);
	SNOVA_PROFILE_BEGIN(gen_P22_p21_pack_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int di = 0; di < v_SNOVA; ++di) {
			snova_pack_block_run_tr_ex(in_8 + (di * rank) * mol_SNOVA8 + mi * (o_SNOVA * rank), mol_SNOVA8,
			                            (const uint8_t*)P21[mi][0][di], v_SNOVA * sq_rank, o_SNOVA);
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_P21_PACK, gen_P22_p21_pack_t);
	SNOVA_PROFILE_BEGIN(gen_P22_p21_mul_t);
	gen_P22_contract(temp_8, in_256, T12, 1);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_P21_MUL, gen_P22_p21_mul_t);
	SNOVA_PROFILE_BEGIN(gen_P22_p21_fold_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int dj = 0; dj < o_SNOVA; ++dj) {
			for (int dk = 0; dk < o_SNOVA; ++dk) {
				const uint8_t* src = temp_8 + (dk * rank) * mol_SNOVA8 + mi * (o_SNOVA * rank) + dj * rank;
#if SNOVA_RANK4_FAST
				xor_block_from_rows_tr(P22[mi][dj][dk], src, mol_SNOVA8);
#else
				for (int i1 = 0; i1 < rank; ++i1) {
					for (int j1 = 0; j1 < rank; ++j1) {
						P22[mi][dj][dk][i1 * rank + j1] ^= src[j1 * mol_SNOVA8 + i1];
					}
				}
#endif
			}
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_P21_FOLD, gen_P22_p21_fold_t);
	SNOVA_PROFILE_BEGIN(gen_P22_compress_t);
	snova_gf16_to_bytes(outP22, (const uint8_t*)P22, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GEN_P22_COMPRESS, gen_P22_compress_t);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_GEN_P22, gen_P22_total_t);
}

#define GAUSS_ROW (m_SNOVA * lsq_SNOVA + 1)
#define GAUSS_ROW32 ((GAUSS_ROW + VEC_LENGTH - 1) / VEC_LENGTH)
#define GAUSS_ROW_mult32 (GAUSS_ROW32 * VEC_LENGTH)
#define VROW_BYTES (vtl_v_len * VEC_LENGTH)
typedef uint8_t lr_pair_t[m_SNOVA][alpha_SNOVA][rank_next32][VROW_BYTES];
static void calc_LR_J_vtl(lr_pair_t L_J_pair, lr_pair_t R_tr_J_pair, const sk_gf16* sk, gf16m_t* X_in_GF16Matrix) {
	alignas(VEC_LENGTH) VECTOR X_J[rank][vtl_v_len] = {0};
	alignas(VEC_LENGTH) VECTOR X_tr_J[rank][vtl_v_len] = {0};
	snova_pack_blocks((VECTOR*)X_J, (const uint8_t*)X_in_GF16Matrix, 1, v_SNOVA);
	snova_pack_blocks_tr((VECTOR*)X_tr_J, (const uint8_t*)X_in_GF16Matrix, 1, v_SNOVA);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		alignas(VEC_LENGTH) VECTOR AxS_tr_256[alpha_SNOVA_next2][rank][vtl_v_len] = {0};
		alignas(VEC_LENGTH) VECTOR Q2xS_tr_256[alpha_SNOVA_next2][rank][vtl_v_len] = {0};
		gf16m_t A_pair[alpha_SNOVA_next2] = {0};
		gf16m_t Q2_pair[alpha_SNOVA_next2] = {0};
		memcpy((uint8_t*)A_pair, (const uint8_t*)sk->Aalpha[mi], sizeof(sk->Aalpha[mi]));
		memcpy((uint8_t*)Q2_pair, (const uint8_t*)sk->Qalpha2[mi], sizeof(sk->Qalpha2[mi]));
		for (int alpha = 0; alpha < alpha_SNOVA_next2; alpha += 2) {
			alignas(VEC_LENGTH) VECTOR AxS_256[rank][vtl_v_len] = {0};
			alignas(VEC_LENGTH) VECTOR Q2xS_256[rank][vtl_v_len] = {0};
			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					VECTOR k1_lh = mtk2_16[get_gf16m(A_pair[alpha], ni, nj) | (get_gf16m(A_pair[alpha + 1], ni, nj) << 4)];
					VECTOR k2_lh = mtk2_16[get_gf16m(Q2_pair[alpha], ni, nj) | (get_gf16m(Q2_pair[alpha + 1], ni, nj) << 4)];
					for (int nk = 0; nk < vtl_v_len; ++nk) {
						AxS_256[ni][nk] ^= VEC_SHUFFLE(k1_lh, X_tr_J[nj][nk]);
						Q2xS_256[ni][nk] ^= VEC_SHUFFLE(k2_lh, X_J[nj][nk]);
					}
				}
			}
			snova_vtl_transpose((VECTOR*)AxS_tr_256[alpha], (const VECTOR*)AxS_256, 1, v_SNOVA);
			snova_vtl_transpose((VECTOR*)Q2xS_tr_256[alpha], (const VECTOR*)Q2xS_256, 1, v_SNOVA);
			for (int ni = 0; ni < rank; ++ni) {
				for (int nk = 0; nk < vtl_v_len; ++nk) {
					AxS_tr_256[alpha + 1][ni][nk] = (AxS_tr_256[alpha][ni][nk] >> 4) & l_mask;
					AxS_tr_256[alpha][ni][nk] &= l_mask;
					Q2xS_tr_256[alpha + 1][ni][nk] = (Q2xS_tr_256[alpha][ni][nk] >> 4) & l_mask;
					Q2xS_tr_256[alpha][ni][nk] &= l_mask;
				}
			}
		}
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mj = i_prime(mi, alpha);
			alignas(VEC_LENGTH) VECTOR L_tr_J_256[rank][vtl_v_len] = {0};
			VECTOR* R_out = (VECTOR*)R_tr_J_pair[mj][alpha];
			for (int nj = 0; nj < rank; ++nj) {
				VECTOR ktab[rank_next32];
				for (int p = 0; p < rank_floor2 / 2; ++p) {
					ktab[p] = mtk2_16[get_gf16m(sk->Qalpha1[mi][alpha], 2 * p, nj) |
					                  (get_gf16m(sk->Qalpha1[mi][alpha], 2 * p + 1, nj) << 4)];
				}
#if rank % 2
				ktab[rank_next32 - 1] = SNOVA_MULTTAB(get_gf16m(sk->Qalpha1[mi][alpha], rank_floor2, nj));
#endif
				for (int nk = 0; nk < vtl_v_len; ++nk) {
					VECTOR src = AxS_tr_256[alpha][nj][nk];
					for (int p = 0; p < rank_next32; ++p) {
						L_tr_J_256[2 * p][nk] ^= VEC_SHUFFLE(ktab[p], src);
					}
				}
			}
			for (int nk = 0; nk < vtl_v_len; ++nk) {
				for (int p = 0; p < rank_next32; ++p) {
					R_out[p * vtl_v_len + nk] = zero256;
				}
			}
			for (int nj = 0; nj < rank; ++nj) {
				VECTOR ktab[rank_next32];
				for (int p = 0; p < rank_floor2 / 2; ++p) {
					ktab[p] = mtk2_16[get_gf16m(sk->Balpha[mi][alpha], nj, 2 * p) |
					                  (get_gf16m(sk->Balpha[mi][alpha], nj, 2 * p + 1) << 4)];
				}
#if rank % 2
				ktab[rank_next32 - 1] = SNOVA_MULTTAB(get_gf16m(sk->Balpha[mi][alpha], nj, rank_floor2));
#endif
				for (int nk = 0; nk < vtl_v_len; ++nk) {
					VECTOR src = Q2xS_tr_256[alpha][nj][nk];
					for (int p = 0; p < rank_next32; ++p) {
						R_out[p * vtl_v_len + nk] ^= VEC_SHUFFLE(ktab[p], src);
					}
				}
			}
			snova_pack_rankpair_tr_rows((uint8_t*)L_J_pair[mj][alpha], (const VECTOR*)L_tr_J_256, v_SNOVA);
		}
	}
}
static void sign_fvv_accum_vtl(gf16m_t Fvv_in_GF16Matrix[m_SNOVA], lr_pair_t L_J_pair, lr_pair_t R_tr_J_pair,
                               const sk_gf16* sk) {
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		gf16m_set_zero(Fvv_in_GF16Matrix[mi]);
	}
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		alignas(VEC_LENGTH) VECTOR F11_J_256[v_SNOVA * rank][vtl_v_len] = {0};
		snova_pack_blocks((VECTOR*)F11_J_256, (const uint8_t*)sk->F11[mi], v_SNOVA, v_SNOVA);
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime = i_prime_inv(mi, alpha);
			alignas(VEC_LENGTH) VECTOR LJxF11J_256[rank][vtl_v_len];
			VECTOR acc[rank_next32][vtl_v_len];
			for (int p = 0; p < rank_next32; ++p) {
				for (int vk = 0; vk < vtl_v_len; ++vk) {
					acc[p][vk] = zero256;
				}
			}
			for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
				VECTOR ktab[rank_next32];
				for (int p = 0; p < rank_next32; ++p) {
					ktab[p] = mtk2_16[L_J_pair[mi][alpha][p][vj]];
				}
				for (int vk = 0; vk < vtl_v_len; ++vk) {
					VECTOR f = F11_J_256[vj][vk];
					for (int p = 0; p < rank_next32; ++p) {
						acc[p][vk] ^= VEC_SHUFFLE(ktab[p], f);
					}
				}
			}
			for (int vk = 0; vk < vtl_v_len; ++vk) {
				for (int p = 0; p < rank_floor2 / 2; ++p) {
					LJxF11J_256[2 * p][vk] = acc[p][vk] & l_mask;
					LJxF11J_256[2 * p + 1][vk] = (acc[p][vk] >> 4) & l_mask;
				}
#if rank % 2
				LJxF11J_256[rank_floor2][vk] = acc[rank_next32 - 1][vk];
#endif
			}
			uint8_t LJxF11JxRJ[rank][rank];
			const VECTOR* R_packed = (const VECTOR*)R_tr_J_pair[mi][alpha];
			for (int vj = 0; vj < rank; ++vj) {
				VECTOR sum[rank];
				for (int vi = 0; vi < rank; ++vi) {
					sum[vi] = zero256;
				}
				for (int vk = 0; vk < vtl_v_len; ++vk) {
					VECTOR packed = R_packed[(vj / 2) * vtl_v_len + vk];
					VECTOR b = (vj & 1) ? ((packed >> 4) & l_mask) : (packed & l_mask);
#if rank == 4
					gf16_32_mul_32_add_4_same_b_values(LJxF11J_256[0][vk], LJxF11J_256[1][vk], LJxF11J_256[2][vk],
					                                   LJxF11J_256[3][vk], &b, &sum[0], &sum[1], &sum[2], &sum[3]);
#else
					for (int vi = 0; vi < rank; ++vi) {
						gf16_32_mul_32_add(&LJxF11J_256[vi][vk], &b, &sum[vi]);
					}
#endif
				}
				for (int vi = 0; vi < rank; ++vi) {
					LJxF11JxRJ[vi][vj] = horizontal_xor_256(sum[vi]);
				}
			}
			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					gf16_t t = get_gf16m(Fvv_in_GF16Matrix[mi_prime], ni, nj);
					set_gf16m(Fvv_in_GF16Matrix[mi_prime], ni, nj, t ^ LJxF11JxRJ[ni][nj]);
				}
			}
		}
	}
}
static void sign_linear_accum_vtl(uint8_t* restrict Temp1, uint8_t* restrict Temp2, lr_pair_t L_J_pair,
                                  lr_pair_t R_tr_J_pair, const sk_gf16* sk) {
	const int o_row32 = vtl_mainRow_x_rank32(o_SNOVA);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		gf16m_t F21_vo[v_SNOVA][o_SNOVA];
		alignas(VEC_LENGTH) VECTOR F12_J_256[v_SNOVA * rank][vtl_o_len] = {0};
		alignas(VEC_LENGTH) VECTOR F21_vo_tr_J_256[v_SNOVA * rank][vtl_o_len] = {0};
		for (int oi = 0; oi < o_SNOVA; ++oi) {
			for (int vi = 0; vi < v_SNOVA; ++vi) {
				gf16m_clone(F21_vo[vi][oi], sk->F21[mi][oi][vi]);
			}
		}
		snova_pack_blocks((VECTOR*)F12_J_256, (const uint8_t*)sk->F12[mi], v_SNOVA, o_SNOVA);
		snova_pack_blocks_tr((VECTOR*)F21_vo_tr_J_256, (const uint8_t*)F21_vo, v_SNOVA, o_SNOVA);
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime_inv = i_prime_inv(mi, alpha);
			alignas(VEC_LENGTH) VECTOR LJxF12J_256[rank][vtl_o_len];
			alignas(VEC_LENGTH) VECTOR LJxF12J_tr_256[rank][vtl_o_len] = {0};
			alignas(VEC_LENGTH) VECTOR LJxF12JxQJ_tr_256[rank][vtl_o_len] = {0};
			VECTOR acc[rank_next32][vtl_o_len];
			for (int p = 0; p < rank_next32; ++p) {
				for (int oi = 0; oi < vtl_o_len; ++oi) {
					acc[p][oi] = zero256;
				}
			}
			for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
				VECTOR ktab[rank_next32];
				for (int p = 0; p < rank_next32; ++p) {
					ktab[p] = mtk2_16[L_J_pair[mi][alpha][p][vj]];
				}
				for (int oi = 0; oi < vtl_o_len; ++oi) {
					VECTOR f = F12_J_256[vj][oi];
					for (int p = 0; p < rank_next32; ++p) {
						acc[p][oi] ^= VEC_SHUFFLE(ktab[p], f);
					}
				}
			}
			for (int oi = 0; oi < vtl_o_len; ++oi) {
				for (int p = 0; p < rank_floor2 / 2; ++p) {
					LJxF12J_256[2 * p][oi] = acc[p][oi] & l_mask;
					LJxF12J_256[2 * p + 1][oi] = (acc[p][oi] >> 4) & l_mask;
				}
#if rank % 2
				LJxF12J_256[rank_floor2][oi] = acc[rank_next32 - 1][oi];
#endif
			}
			snova_vtl_transpose((VECTOR*)LJxF12J_tr_256, (const VECTOR*)LJxF12J_256, 1, o_SNOVA);
			for (int ri = 0; ri < rank; ++ri) {
				for (int rj = 0; rj < rank; ++rj) {
					VECTOR k = SNOVA_MULTTAB(get_gf16m(sk->Qalpha2[mi_prime_inv][alpha], rj, ri));
					for (int oi = 0; oi < vtl_o_len; ++oi) {
						LJxF12JxQJ_tr_256[ri][oi] ^= VEC_SHUFFLE(k, LJxF12J_tr_256[rj][oi]);
					}
				}
			}
			for (int ti2 = 0; ti2 < l_SNOVA; ++ti2) {
				for (int tj2 = 0; tj2 < l_SNOVA; ++tj2) {
					VECTOR k = SNOVA_MULTTAB(sk->Balpha[mi_prime_inv][alpha][tj2 * rank + ti2]);
					VECTOR* Temp1_256 = (VECTOR*)Temp1;
					const VECTOR* AJ_256 = (const VECTOR*)LJxF12JxQJ_tr_256;
					for (int tj1 = 0; tj1 < l_SNOVA; ++tj1) {
						for (int toi = 0; toi < o_row32; ++toi) {
							Temp1_256[((mi_prime_inv * lsq_SNOVA + ti2 * rank + tj2) * l_SNOVA + tj1) * o_row32 + toi] ^=
							    VEC_SHUFFLE(k, AJ_256[(tj1 * o_row32) + toi]);
						}
					}
				}
			}
			alignas(VEC_LENGTH) VECTOR RTRJ_x_F21TRJ_256[rank][vtl_o_len];
			alignas(VEC_LENGTH) VECTOR F21JxRJ_256[rank][vtl_o_len] = {0};
			alignas(VEC_LENGTH) VECTOR Q1xF21JxRJ_256[rank][vtl_o_len] = {0};
			for (int p = 0; p < rank_next32; ++p) {
				for (int oi = 0; oi < vtl_o_len; ++oi) {
					acc[p][oi] = zero256;
				}
			}
			for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
				VECTOR ktab[rank_next32];
				for (int p = 0; p < rank_next32; ++p) {
					ktab[p] = mtk2_16[R_tr_J_pair[mi][alpha][p][vj]];
				}
				for (int oi = 0; oi < vtl_o_len; ++oi) {
					VECTOR f = F21_vo_tr_J_256[vj][oi];
					for (int p = 0; p < rank_next32; ++p) {
						acc[p][oi] ^= VEC_SHUFFLE(ktab[p], f);
					}
				}
			}
			for (int oi = 0; oi < vtl_o_len; ++oi) {
				for (int p = 0; p < rank_floor2 / 2; ++p) {
					RTRJ_x_F21TRJ_256[2 * p][oi] = acc[p][oi] & l_mask;
					RTRJ_x_F21TRJ_256[2 * p + 1][oi] = (acc[p][oi] >> 4) & l_mask;
				}
#if rank % 2
				RTRJ_x_F21TRJ_256[rank_floor2][oi] = acc[rank_next32 - 1][oi];
#endif
			}
			snova_vtl_transpose((VECTOR*)F21JxRJ_256, (const VECTOR*)RTRJ_x_F21TRJ_256, 1, o_SNOVA);
			for (int ri = 0; ri < rank; ++ri) {
				for (int rj = 0; rj < rank; ++rj) {
					VECTOR k = SNOVA_MULTTAB(get_gf16m(sk->Qalpha1[mi_prime_inv][alpha], ri, rj));
					for (int oi = 0; oi < vtl_o_len; ++oi) {
						Q1xF21JxRJ_256[ri][oi] ^= VEC_SHUFFLE(k, F21JxRJ_256[rj][oi]);
					}
				}
			}
			for (int ti1 = 0; ti1 < l_SNOVA; ++ti1) {
				for (int tj2 = 0; tj2 < l_SNOVA; ++tj2) {
					VECTOR k = SNOVA_MULTTAB(sk->Aalpha[mi_prime_inv][alpha][ti1 * rank + tj2]);
					VECTOR* Temp2_256 = (VECTOR*)Temp2;
					const VECTOR* AJ_256 = (const VECTOR*)Q1xF21JxRJ_256;
					for (int tj1 = 0; tj1 < l_SNOVA; ++tj1) {
						for (int toi = 0; toi < o_row32; ++toi) {
							Temp2_256[((mi_prime_inv * lsq_SNOVA + ti1 * l_SNOVA + tj2) * l_SNOVA + tj1) * o_row32 + toi] ^=
							    VEC_SHUFFLE(k, AJ_256[(tj1 * o_row32) + toi]);
						}
					}
				}
			}
		}
	}
}
static int sign_gauss_eliminate(uint8_t Gauss[m_SNOVA * lsq_SNOVA][GAUSS_ROW_mult32]) {
	int flag_redo = 0;
	for (int mi2 = 0; mi2 < m_SNOVA * lsq_SNOVA; ++mi2) {
		const int kstart = (mi2 / VEC_LENGTH) * VEC_LENGTH;
		int swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
		for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
			int8x16_t swap256 = vld1q_dup_s8((int8_t*)&swap);
#else
			VECTOR swap256 = _mm256_set1_epi32(swap);
#endif
			VECTOR* gdest = (VECTOR*)&Gauss[mi2][kstart];
			VECTOR* gsource = (VECTOR*)&Gauss[j2][kstart];
			for (int k2 = 0; k2 < GAUSS_ROW32 - kstart / VEC_LENGTH; ++k2) {
				gdest[k2] ^= gsource[k2] & swap256;
			}
			swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
		}
		flag_redo |= swap;
		VECTOR res256[1];
		VECTOR g256[1];
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
		g256[0] = vld1q_dup_u8(&Gauss[mi2][mi2]);
#else
		g256[0] = _mm256_set1_epi8(Gauss[mi2][mi2]);
#endif
		res256[0] = VEC_SHUFFLE(vector_inv_table, g256[0]);
		uint8_t t_GF16 = VEC_BYTE(res256, 0);
		VECTOR pivot_tab = SNOVA_MULTTAB(t_GF16);
		for (int k = kstart; k < GAUSS_ROW_mult32; k += VEC_LENGTH) {
			gf16_32_mul_tab((VECTOR*)(Gauss[mi2] + k), pivot_tab, (VECTOR*)(Gauss[mi2] + k));
		}
		for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2) {
			VECTOR row_tab = SNOVA_MULTTAB(Gauss[j2][mi2]);
			for (int k2 = kstart; k2 < GAUSS_ROW_mult32; k2 += VEC_LENGTH) {
				gf16_32_mul_tab_add((VECTOR*)(Gauss[mi2] + k2), row_tab, (VECTOR*)(Gauss[j2] + k2));
			}
		}
	}
	return flag_redo;
}
int snova_sign_core(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, const uint8_t* array_salt,
                     const sk_gf16* sk) {
	alignas(VEC_LENGTH) uint8_t Gauss[m_SNOVA * lsq_SNOVA][GAUSS_ROW_mult32];
	gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
	gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
	gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
	gf16m_t signature_in_GF16Matrix[n_SNOVA];
	uint8_t signed_hash[bytes_hash];
	uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];
	int flag_redo = 1;
	uint8_t num_sign = 0;
	SNOVA_PROFILE_BEGIN(sign_prehash_t);
	snova_signed_hash(digest, bytes_digest, sk->pt_public_key_seed, array_salt, signed_hash);
	SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_PREHASH, sign_prehash_t);
	SNOVA_PROFILE_BEGIN(sign_hash_to_gf_t);
	snova_bytes_to_gf16(signed_hash, hash_in_GF16, GF16s_hash);
	SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_HASH_TO_GF, sign_hash_to_gf_t);
	do {
		SNOVA_PROFILE_BEGIN(sign_gauss_init_t);
		memset(Gauss, 0, sizeof(Gauss));
		num_sign++;
		if (num_sign == 255) {
			memset(pt_signature, 0, bytes_sig_with_salt);
			SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_GAUSS_INIT, sign_gauss_init_t);
			return -1;
		}
		flag_redo = 0;
		for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
			Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
		}
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_GAUSS_INIT, sign_gauss_init_t);
		SNOVA_PROFILE_BEGIN(sign_vinegar_t);
		Keccak_HashInstance hashInstance;
		Keccak_HashInitialize_SHAKE256(&hashInstance);
		Keccak_HashUpdate(&hashInstance, sk->pt_private_key_seed, 8 * seed_length_private);
		Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
		Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
		Keccak_HashUpdate(&hashInstance, &num_sign, 8);
		Keccak_HashFinal(&hashInstance, NULL);
		Keccak_HashSqueeze(&hashInstance, vinegar_in_byte, 8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));
		snova_bytes_to_gf16(vinegar_in_byte, (uint8_t*)X_in_GF16Matrix, v_SNOVA * lsq_SNOVA);
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_VINEGAR, sign_vinegar_t);
		alignas(VEC_LENGTH) lr_pair_t L_J_pair;
		alignas(VEC_LENGTH) lr_pair_t R_tr_J_pair;
		SNOVA_PROFILE_BEGIN(sign_lr_t);
		calc_LR_J_vtl(L_J_pair, R_tr_J_pair, sk, X_in_GF16Matrix);
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_LR, sign_lr_t);
		SNOVA_PROFILE_BEGIN(sign_fvv_t);
		sign_fvv_accum_vtl(Fvv_in_GF16Matrix, L_J_pair, R_tr_J_pair, sk);
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_FVV, sign_fvv_t);
		SNOVA_PROFILE_BEGIN(sign_fvv_rhs_t);
		for (int i = 0; i < m_SNOVA; ++i) {
			for (int j = 0; j < rank; ++j) {
				for (int k = 0; k < rank; ++k) {
					int index1 = i * lsq_SNOVA + j * rank + k;
					int index2 = m_SNOVA * lsq_SNOVA;
					Gauss[index1][index2] = gf16_get_add(Gauss[index1][index2], get_gf16m(Fvv_in_GF16Matrix[i], j, k));
				}
			}
		}
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_FVV_RHS, sign_fvv_rhs_t);
		alignas(VEC_LENGTH) uint8_t Temp1[o_SNOVA * l_SNOVA * lsq_SNOVA * vtl_mainRow_x_rank(o_SNOVA)] = {0};
		alignas(VEC_LENGTH) uint8_t Temp2[o_SNOVA * l_SNOVA * lsq_SNOVA * vtl_mainRow_x_rank(o_SNOVA)] = {0};
		SNOVA_PROFILE_BEGIN(sign_linear_accum_t);
		sign_linear_accum_vtl(Temp1, Temp2, L_J_pair, R_tr_J_pair, sk);
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_LINEAR_ACCUM, sign_linear_accum_t);
		SNOVA_PROFILE_BEGIN(sign_linear_fold1_t);
		for (int mi_prime_inv = 0; mi_prime_inv < o_SNOVA; ++mi_prime_inv)
			for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
				for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
					for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
						for (int oi = 0; oi < o_SNOVA; ++oi)
							for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
							Gauss[mi_prime_inv * lsq_SNOVA + ti1 * rank + ti2][oi * lsq_SNOVA + tj1 * rank + tj2] ^=
							    Temp1[((mi_prime_inv * lsq_SNOVA + ti2 * rank + tj2) * l_SNOVA + tj1) *
							          vtl_mainRow_x_rank(o_SNOVA) +
							          (oi * rank + ti1)];
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_LINEAR_FOLD1, sign_linear_fold1_t);
		SNOVA_PROFILE_BEGIN(sign_linear_fold2_t);
		for (int mi_prime_inv = 0; mi_prime_inv < o_SNOVA; ++mi_prime_inv)
			for (int oi = 0; oi < o_SNOVA; ++oi)
				for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
					for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
						for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
							for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
							Gauss[mi_prime_inv * lsq_SNOVA + ti1 * rank + ti2][oi * lsq_SNOVA + tj1 * rank + tj2] ^=
							    Temp2[((mi_prime_inv * lsq_SNOVA + ti1 * l_SNOVA + tj2) * l_SNOVA + tj1) *
							          vtl_mainRow_x_rank(o_SNOVA) +
							          (oi * rank + ti2)];
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_LINEAR_FOLD2, sign_linear_fold2_t);
		SNOVA_PROFILE_BEGIN(sign_gauss_t);
		flag_redo |= sign_gauss_eliminate(Gauss);
		SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_GAUSS, sign_gauss_t);
	} while (flag_redo);
	alignas(VEC_LENGTH) uint8_t solution[GAUSS_ROW_mult32] = {0};
	uint32_t xSolution[m_SNOVA * lsq_SNOVA] = {0};
	uint8_t t_GF16 = 0;
	uint8_t Gauss_last_col;
	SNOVA_PROFILE_BEGIN(sign_backsub_t);
	for (int mil2 = m_SNOVA * lsq_SNOVA - 1; mil2 >= 0; --mil2) {
		Gauss_last_col = Gauss[mil2][m_SNOVA * lsq_SNOVA];
		VECTOR t_GF16_256 = zero256;
		Gauss[mil2][m_SNOVA * lsq_SNOVA] = 0;
		int kstart = ((mil2 + 1) / VEC_LENGTH) * VEC_LENGTH;
		for (int k2 = kstart; k2 < GAUSS_ROW_mult32; k2 += VEC_LENGTH) {
			gf16_32_mul_32_add((VECTOR*)(&Gauss[mil2][k2]), (VECTOR*)(&solution[k2]), &t_GF16_256);
		}
		uint8_t* t_GF16_256_8_ptr = (uint8_t*)(&t_GF16_256);
		for (int k2 = 0; k2 < VEC_LENGTH; k2++) {
			t_GF16 ^= t_GF16_256_8_ptr[k2];
		}
		solution[mil2] = Gauss_last_col ^ t_GF16;
		t_GF16 = 0;
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_BACKSUB, sign_backsub_t);
	SNOVA_PROFILE_BEGIN(sign_build_sig_t);
	for (int mi2 = m_SNOVA * lsq_SNOVA - 1; mi2 >= 0; --mi2) {
		xSolution[mi2] = gf16_from_nibble(solution[mi2]);
	}
	uint32_t xT12[v_SNOVA][o_SNOVA][lsq_SNOVA] = {0};
	for (int dj = 0; dj < v_SNOVA; ++dj)
		for (int dk = 0; dk < o_SNOVA; ++dk)
			for (int idx = 0; idx < lsq_SNOVA; ++idx) {
				xT12[dj][dk][idx] = gf16_from_nibble(sk->T12[dj][dk][idx]);
			}
	uint32_t xSig[lsq_SNOVA] = {0};
	for (int dj = 0; dj < v_SNOVA; ++dj) {
		for (int dk = 0; dk < o_SNOVA; ++dk)
			for (int i1 = 0; i1 < l_SNOVA; ++i1)
				for (int j1 = 0; j1 < l_SNOVA; ++j1)
					for (int k1 = 0; k1 < l_SNOVA; ++k1) {
						xSig[i1 * l_SNOVA + j1] ^=
						    xT12[dj][dk][i1 * l_SNOVA + k1] * xSolution[dk * lsq_SNOVA + k1 * l_SNOVA + j1];
					}
		for (int idx = 0; idx < lsq_SNOVA; ++idx) {
			signature_in_GF16Matrix[dj][idx] = X_in_GF16Matrix[dj][idx] ^ gf16_to_nibble(xSig[idx]);
		}
		memset(xSig, 0, sizeof(xSig));
	}
	for (int index = 0; index < o_SNOVA; ++index)
		for (int idx = 0; idx < lsq_SNOVA; ++idx) {
			signature_in_GF16Matrix[v_SNOVA + index][idx] = solution[index * lsq_SNOVA + idx];
		}
	SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_BUILD_SIG, sign_build_sig_t);
	SNOVA_PROFILE_BEGIN(sign_pack_t);
	snova_gf16_to_bytes(pt_signature, (const gf16_t*)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
	for (int i1 = 0; i1 < bytes_salt; ++i1) {
		pt_signature[bytes_signature + i1] = array_salt[i1];
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_SIGN_PACK, sign_pack_t);
	return 0;
}

#define NROW_BYTES (vtl_len * VEC_LENGTH)
typedef uint8_t verify_pair_t[m_SNOVA][alpha_SNOVA][rank_next32][NROW_BYTES];
#ifndef SNOVA_VERIFY_PACKED_SECONDARY
#define SNOVA_VERIFY_PACKED_SECONDARY 1
#endif
#ifndef SNOVA_VERIFY_PACKED_MASK_SCHEDULE
#define SNOVA_VERIFY_PACKED_MASK_SCHEDULE 1
#endif
#if SNOVA_VERIFY_PACKED_SECONDARY && rank == 4
#define VERIFY_H_COLS rank_next32
#else
#define VERIFY_H_COLS rank
#endif
static void verify_eval_lr(const public_key_expand* pkx, const VECTOR S_J_256[rank][vtl_len],
                           const VECTOR S_tr_J_256[rank][vtl_len], verify_pair_t L_J_pair, verify_pair_t R_tr_J_pair) {
	/* The caller owns the coarse eval_lr timer; this helper reports only micro stages. */
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		VECTOR AxS_tr_256[alpha_SNOVA_next2][rank][vtl_len] = {0};
		VECTOR Q2xS_tr_256[alpha_SNOVA_next2][rank][vtl_len] = {0};
		gf16m_t Aalpha[alpha_SNOVA_next2] = {0};
		gf16m_t Qalpha2[alpha_SNOVA_next2] = {0};
		memcpy((uint8_t*)Aalpha, (const uint8_t*)(pkx->map1.Aalpha[mi]), sizeof(pkx->map1.Aalpha[mi]));
		memcpy((uint8_t*)Qalpha2, (const uint8_t*)(pkx->map1.Qalpha2[mi]), sizeof(pkx->map1.Qalpha2[mi]));
		for (int alpha = 0; alpha < alpha_SNOVA_next2; alpha += 2) {
			VECTOR AxS_256[rank][vtl_len] = {0};
			VECTOR Q2xS_256[rank][vtl_len] = {0};
			VECTOR A_tab[rank][rank];
			VECTOR Q2_tab[rank][rank];
			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					A_tab[ni][nj] = mtk2_16[get_gf16m(Aalpha[alpha], ni, nj) |
					                         (get_gf16m(Aalpha[alpha + 1], ni, nj) << 4)];
					Q2_tab[ni][nj] = mtk2_16[get_gf16m(Qalpha2[alpha], ni, nj) |
					                          (get_gf16m(Qalpha2[alpha + 1], ni, nj) << 4)];
				}
			}
			SNOVA_MICRO_PROFILE_BEGIN(verify_lr_contract_t);
			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					for (int nk = 0; nk < vtl_len; ++nk) {
						AxS_256[ni][nk] ^= VEC_SHUFFLE(A_tab[ni][nj], S_tr_J_256[nj][nk]);
						Q2xS_256[ni][nk] ^= VEC_SHUFFLE(Q2_tab[ni][nj], S_J_256[nj][nk]);
					}
				}
			}
			SNOVA_MICRO_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_LR_CONTRACT, verify_lr_contract_t);
			SNOVA_MICRO_PROFILE_BEGIN(verify_lr_pack_t);
			snova_vtl_transpose((VECTOR*)AxS_tr_256[alpha], (const VECTOR*)AxS_256, 1, n_SNOVA);
			snova_vtl_transpose((VECTOR*)Q2xS_tr_256[alpha], (const VECTOR*)Q2xS_256, 1, n_SNOVA);
			for (int ni = 0; ni < rank; ++ni) {
				for (int nk = 0; nk < vtl_len; ++nk) {
					AxS_tr_256[alpha + 1][ni][nk] = (AxS_tr_256[alpha][ni][nk] >> 4) & l_mask;
					AxS_tr_256[alpha][ni][nk] &= l_mask;
					Q2xS_tr_256[alpha + 1][ni][nk] = (Q2xS_tr_256[alpha][ni][nk] >> 4) & l_mask;
					Q2xS_tr_256[alpha][ni][nk] &= l_mask;
				}
			}
			SNOVA_MICRO_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_LR_PACK, verify_lr_pack_t);
		}
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mj = i_prime(mi, alpha);
			VECTOR L_tr_J_256[rank][vtl_len] = {0};
			VECTOR* R_out = (VECTOR*)R_tr_J_pair[mj][alpha];
			VECTOR Q1_tab[rank][rank_next32];
			VECTOR B_tab[rank][rank_next32];
			for (int nj = 0; nj < rank; ++nj) {
				for (int p = 0; p < rank_floor2 / 2; ++p) {
					Q1_tab[nj][p] = mtk2_16[get_gf16m(pkx->map1.Qalpha1[mi][alpha], 2 * p, nj) |
					                        (get_gf16m(pkx->map1.Qalpha1[mi][alpha], 2 * p + 1, nj) << 4)];
					B_tab[nj][p] = mtk2_16[get_gf16m(pkx->map1.Balpha[mi][alpha], nj, 2 * p) |
					                       (get_gf16m(pkx->map1.Balpha[mi][alpha], nj, 2 * p + 1) << 4)];
				}
#if rank % 2
				Q1_tab[nj][rank_next32 - 1] =
				    SNOVA_MULTTAB(get_gf16m(pkx->map1.Qalpha1[mi][alpha], rank_floor2, nj));
				B_tab[nj][rank_next32 - 1] =
				    SNOVA_MULTTAB(get_gf16m(pkx->map1.Balpha[mi][alpha], nj, rank_floor2));
#endif
			}
			for (int nj = 0; nj < rank; ++nj) {
				for (int nk = 0; nk < vtl_len; ++nk) {
					VECTOR src = AxS_tr_256[alpha][nj][nk];
					for (int p = 0; p < rank_next32; ++p) {
						L_tr_J_256[2 * p][nk] ^= VEC_SHUFFLE(Q1_tab[nj][p], src);
					}
				}
			}
			for (int nk = 0; nk < vtl_len; ++nk) {
				for (int p = 0; p < rank_next32; ++p) {
					R_out[p * vtl_len + nk] = zero256;
				}
			}
			for (int nj = 0; nj < rank; ++nj) {
				for (int nk = 0; nk < vtl_len; ++nk) {
					VECTOR src = Q2xS_tr_256[alpha][nj][nk];
					for (int p = 0; p < rank_next32; ++p) {
						R_out[p * vtl_len + nk] ^= VEC_SHUFFLE(B_tab[nj][p], src);
					}
				}
			}
			snova_pack_rankpair_tr_rows((uint8_t*)L_J_pair[mj][alpha], (const VECTOR*)L_tr_J_256, n_SNOVA);
		}
	}
}
static inline void verify_eval_ljpj(const uint8_t* restrict LJ_pair, const VECTOR (*restrict P_J_256)[vtl_len],
                                    VECTOR LJxPJ_nibble[rank_next32][vtl_len]) {
	/* The caller owns the coarse eval_ljpj timer; this helper reports its mul slice. */
#if rank == 4
	const uint8_t* restrict LJ0 = LJ_pair;
	const uint8_t* restrict LJ1 = LJ_pair + NROW_BYTES;
	VECTOR* restrict out0 = LJxPJ_nibble[0];
	VECTOR* restrict out1 = LJxPJ_nibble[1];
	SNOVA_MICRO_PROFILE_BEGIN(verify_ljpj_mul_t);
	for (int nj = 0; nj < vtl_len; ++nj) {
		out0[nj] = zero256;
		out1[nj] = zero256;
	}
	for (int nk = 0; nk < n_SNOVA * rank; ++nk) {
		VECTOR k0 = mtk2_16[LJ0[nk]];
		VECTOR k1 = mtk2_16[LJ1[nk]];
		const VECTOR* restrict P_row = P_J_256[nk];
		for (int nj = 0; nj < vtl_len; ++nj) {
			VECTOR p = P_row[nj];
			out0[nj] ^= VEC_SHUFFLE(k0, p);
			out1[nj] ^= VEC_SHUFFLE(k1, p);
		}
	}
	SNOVA_MICRO_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_LJPJ_MUL, verify_ljpj_mul_t);
#else
	for (int p = 0; p < rank_next32; ++p) {
		for (int nj = 0; nj < vtl_len; ++nj) {
			LJxPJ_nibble[p][nj] = zero256;
		}
	}
	for (int nk = 0; nk < n_SNOVA * rank; ++nk) {
		VECTOR ktab[rank_next32];
		for (int p = 0; p < rank_next32; ++p) {
			ktab[p] = mtk2_16[LJ_pair[p * NROW_BYTES + nk]];
		}
		const VECTOR* P_row = P_J_256[nk];
		for (int nj = 0; nj < vtl_len; ++nj) {
			VECTOR pv = P_row[nj];
			for (int p = 0; p < rank_next32; ++p) {
				LJxPJ_nibble[p][nj] ^= VEC_SHUFFLE(ktab[p], pv);
			}
		}
	}
#endif
}
static inline void verify_eval_secondary(const VECTOR LJxPJ_nibble[rank_next32][vtl_len],
                                         const uint8_t* restrict R_pair, VECTOR h[rank][VERIFY_H_COLS]) {
	/* The caller owns the coarse eval_secondary timer; only packed-a substeps are sampled here. */
	const VECTOR* R_vec = (const VECTOR*)R_pair;
#if SNOVA_VERIFY_PACKED_SECONDARY && rank == 4
	for (int nk = 0; nk < vtl_len; ++nk) {
		const VECTOR a01 = LJxPJ_nibble[0][nk];
		const VECTOR a23 = LJxPJ_nibble[1][nk];
#if SNOVA_VERIFY_PACKED_MASK_SCHEDULE == 1
		SNOVA_MICRO_PROFILE_BEGIN(verify_secondary_packa_t);
		VECTOR m01_1, m01_2, m01_4, m01_8;
		VECTOR m01_1_hi, m01_2_hi, m01_4_hi, m01_8_hi;
		VECTOR m23_1, m23_2, m23_4, m23_8;
		VECTOR m23_1_hi, m23_2_hi, m23_4_hi, m23_8_hi;
		gf16_32_prepare_packed_a_masks(a01,
		                              &m01_1, &m01_2, &m01_4, &m01_8,
		                              &m01_1_hi, &m01_2_hi, &m01_4_hi, &m01_8_hi);
		gf16_32_prepare_packed_a_masks(a23,
		                              &m23_1, &m23_2, &m23_4, &m23_8,
		                              &m23_1_hi, &m23_2_hi, &m23_4_hi, &m23_8_hi);
		SNOVA_MICRO_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_SECONDARY_PACKA, verify_secondary_packa_t);
		for (int p = 0; p < rank_next32; ++p) {
			VECTOR b_packed = R_vec[p * vtl_len + nk];
			SNOVA_MICRO_PROFILE_BEGIN(verify_secondary_mul_t);
			gf16_32_mul_32_add_4_same_packed_a_values_masked(
			    &b_packed,
			    m01_1, m01_2, m01_4, m01_8,
			    m01_1_hi, m01_2_hi, m01_4_hi, m01_8_hi,
			    m23_1, m23_2, m23_4, m23_8,
			    m23_1_hi, m23_2_hi, m23_4_hi, m23_8_hi,
			    &h[0][p], &h[1][p], &h[2][p], &h[3][p]);
			SNOVA_MICRO_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_SECONDARY_MUL, verify_secondary_mul_t);
		}
#elif SNOVA_VERIFY_PACKED_MASK_SCHEDULE == 2
		for (int p = 0; p < rank_next32; ++p) {
			VECTOR b_packed = R_vec[p * vtl_len + nk];
			SNOVA_MICRO_PROFILE_BEGIN(verify_secondary_mul_t);
			gf16_32_mul_32_add_4_same_packed_a_values(
			    a01, a23, &b_packed,
			    &h[0][p], &h[1][p], &h[2][p], &h[3][p]);
			SNOVA_MICRO_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_SECONDARY_MUL, verify_secondary_mul_t);
		}
#else
#error "SNOVA_VERIFY_PACKED_MASK_SCHEDULE must be 1 (precompute) or 2 (streaming)"
#endif
	}
#else
	for (int nk = 0; nk < vtl_len; ++nk) {
		VECTOR a[rank];
		for (int p = 0; p < rank_floor2 / 2; ++p) {
			a[2 * p] = LJxPJ_nibble[p][nk] & l_mask;
			a[2 * p + 1] = (LJxPJ_nibble[p][nk] >> 4) & l_mask;
		}
#if rank % 2
		a[rank_floor2] = LJxPJ_nibble[rank_next32 - 1][nk];
#endif
		for (int nj = 0; nj < rank; ++nj) {
			VECTOR packed = R_vec[(nj / 2) * vtl_len + nk];
			VECTOR b = (nj & 1) ? ((packed >> 4) & l_mask) : (packed & l_mask);
#if rank == 4
			gf16_32_mul_32_add_4_same_b_values(a[0], a[1], a[2], a[3], &b, &h[0][nj], &h[1][nj], &h[2][nj], &h[3][nj]);
#else
			for (int ni = 0; ni < rank; ++ni) {
				gf16_32_mul_32_add(&a[ni], &b, &h[ni][nj]);
			}
#endif
		}
	}
#endif
}

#ifndef SNOVA_VERIFY_DIRECT_OFFICIAL
#define SNOVA_VERIFY_DIRECT_OFFICIAL 0
#endif

static inline gf16_t verify_public_p_entry(const public_key_expand* restrict pkx, int mi, int ni, int nj, int row,
                                           int col) {
	if (ni < v_SNOVA) {
		if (nj < v_SNOVA) {
			return get_gf16m(pkx->map1.P11[mi][ni][nj], row, col);
		}
		return get_gf16m(pkx->map1.P12[mi][ni][nj - v_SNOVA], row, col);
	}
	if (nj < v_SNOVA) {
		return get_gf16m(pkx->map1.P21[mi][ni - v_SNOVA][nj], row, col);
	}
	return get_gf16m(pkx->P22[mi][ni - v_SNOVA][nj - v_SNOVA], row, col);
}

static inline void verify_copy_public_blocks(expanded_PK* restrict dst,
                                              const public_key_expand* restrict src) {
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < n_SNOVA; ++ni) {
			for (int nj = 0; nj < n_SNOVA; ++nj) {
				for (int row = 0; row < l_SNOVA; ++row) {
					for (int col = 0; col < l_SNOVA; ++col) {
						dst->P[((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2 +
						       row * l_SNOVA + col] =
						    verify_public_p_entry(src, mi, ni, nj, row, col);
					}
				}
			}
		}
	}
}

static inline void verify_copy_public_blocks_mtk2(uint8_t* restrict dst,
                                                   const public_key_expand* restrict src) {
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < n_SNOVA; ++ni) {
			for (int nj = 0; nj < n_SNOVA; ++nj) {
				for (int row = 0; row < l_SNOVA; ++row) {
					for (int col = 0; col < l_SNOVA; ++col) {
					dst[((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2 +
					    row * l_SNOVA + col] =
					    verify_public_p_entry(src, mi, ni, nj, row, col);
					}
				}
			}
		}
	}
}

#if SNOVA_VERIFY_DIRECT_OFFICIAL && !FIXED_ABQ
static void evaluation_vector_direct_official(gf16m_t* restrict hash_in_GF16Matrix,
                                              const public_key_expand* restrict pkx,
                                              gf16m_t* restrict signature_in_GF16Matrix) {
	uint8_t whipped_sig[n_SNOVA * l_SNOVA * SNOVA_lr] = {0};
	uint8_t sum_t1[m_SNOVA * SNOVA_l2 * SNOVA_r2] = {0};
	uint8_t sum_t1s[m_SNOVA * SNOVA_l2 * SNOVA_r2] = {0};
	uint8_t hash_in_GF[m_SNOVA * SNOVA_lr] = {0};

	for (int ab = 0; ab < l_SNOVA; ++ab) {
		for (int idx = 0; idx < n_SNOVA; ++idx) {
			for (int i1 = 0; i1 < l_SNOVA; ++i1) {
				for (int j1 = 0; j1 < rank; ++j1) {
					uint8_t acc = 0;
					for (int k1 = 0; k1 < l_SNOVA; ++k1) {
						acc ^= mt(get_gf16m(snova_S[ab], i1, k1), get_gf16m(signature_in_GF16Matrix[idx], k1, j1));
					}
					whipped_sig[idx * l_SNOVA * SNOVA_lr + i1 * SNOVA_lr + ab * rank + j1] = acc;
				}
			}
		}
	}

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < n_SNOVA; ++ni) {
			uint8_t sum_t0[l_SNOVA * SNOVA_lr] = {0};

			for (int nj = 0; nj < n_SNOVA; ++nj) {
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					for (int i1 = 0; i1 < l_SNOVA; ++i1) {
						const gf16_t p = verify_public_p_entry(pkx, mi, ni, nj, i1, k1);
						if (p == 0) {
							continue;
						}
						for (int b1 = 0; b1 < SNOVA_lr; ++b1) {
							sum_t0[i1 * SNOVA_lr + b1] ^=
							    mt(p, whipped_sig[nj * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + b1]);
						}
					}
				}
			}

			for (int a1 = 0; a1 < l_SNOVA; ++a1) {
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					for (int i1 = 0; i1 < rank; ++i1) {
						const gf16_t s = whipped_sig[ni * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + a1 * rank + i1];
						if (s == 0) {
							continue;
						}
						for (int b1 = 0; b1 < SNOVA_lr; ++b1) {
							sum_t1[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr + b1] ^=
							    mt(s, sum_t0[k1 * SNOVA_lr + b1]);
						}
					}
				}
			}
		}
	}

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int a1 = 0; a1 < l_SNOVA; ++a1) {
			for (int i1 = 0; i1 < rank; ++i1) {
				for (int b1 = 0; b1 < l_SNOVA; ++b1) {
					for (int j1 = 0; j1 < rank; ++j1) {
						sum_t1s[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + b1 * SNOVA_r2 + i1 * rank + j1] =
						    sum_t1[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr + b1 * rank + j1];
					}
				}
			}
		}
	}

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime = i_prime(mi, alpha);
			uint8_t gfm_temp1[SNOVA_r2] = {0};
			uint8_t gfm_temp2[SNOVA_lr] = {0};

			for (int a1 = 0; a1 < l_SNOVA; ++a1) {
				uint8_t gfm_temp0[SNOVA_r2] = {0};
				const gf16_t q1 = pkx->q1coeff[mi][alpha][a1];

				for (int i1 = 0; i1 < rank; ++i1) {
					for (int b1 = 0; b1 < l_SNOVA; ++b1) {
						const gf16_t q2 = pkx->q2coeff[mi][alpha][b1];
						if (q2 == 0) {
							continue;
						}
						for (int j1 = 0; j1 < rank; ++j1) {
							gfm_temp0[i1 * rank + j1] ^=
							    mt(q2, sum_t1s[(mi_prime * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 +
							                    b1 * SNOVA_r2 + i1 * rank + j1]);
						}
					}
				}

				if (q1 != 0) {
					for (int idx = 0; idx < SNOVA_r2; ++idx) {
						gfm_temp1[idx] ^= mt(gfm_temp0[idx], q1);
					}
				}
			}

			for (int i1 = 0; i1 < rank; ++i1) {
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					for (int k1 = 0; k1 < rank; ++k1) {
						gfm_temp2[i1 * l_SNOVA + j1] ^=
						    mt(gfm_temp1[i1 * rank + k1], get_gf16m(pkx->map1.Balpha[mi][alpha], k1, j1));
					}
				}
			}

			for (int i1 = 0; i1 < rank; ++i1) {
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					for (int k1 = 0; k1 < rank; ++k1) {
						hash_in_GF[mi * SNOVA_lr + i1 * l_SNOVA + j1] ^=
						    mt(get_gf16m(pkx->map1.Aalpha[mi][alpha], i1, k1), gfm_temp2[k1 * l_SNOVA + j1]);
					}
				}
			}
		}
	}

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int idx = 0; idx < SNOVA_lr; ++idx) {
			hash_in_GF16Matrix[mi][idx] = hash_in_GF[mi * SNOVA_lr + idx];
		}
	}
}
#endif

#if !SNOVA_VERIFY_VTL && !SNOVA_VERIFY_DIRECT_MTK2
static int snova_verify_official_after_pk_expand_scalar(const expanded_PK* pkx, const uint8_t* sig,
                                                        const uint8_t* digest, const size_t len_digest) {
	gf16_t signature_in_GF[n_SNOVA * SNOVA_lr];
	SNOVA_PROFILE_BEGIN(verify_unpack_t);
	if (snova_bytes_to_gf16(sig, signature_in_GF, NUMGF_SIGNATURE)) {
		SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);
		return -1;
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);

	const uint8_t* salt = sig + bytes_signature;
	uint8_t signed_hash[bytes_hash];
	uint8_t signed_gf[GF16_HASH] = {0};
	SNOVA_PROFILE_BEGIN(verify_hash_t);
	snova_signed_hash(digest, len_digest, pkx->pk_seed, salt, signed_hash);
#if (o_SNOVA * l_SNOVA) & 0x1 == 1
	signed_hash[bytes_hash - 1] &= 0x0f;
#endif
	snova_bytes_to_gf16(signed_hash, signed_gf, GF16_HASH);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_HASH, verify_hash_t);

	uint16_t whipped_sig[l_SNOVA * n_SNOVA * SNOVA_lr] = {0};
	for (int ab = 0; ab < l_SNOVA; ++ab)
		for (int idx = 0; idx < n_SNOVA; ++idx)
			for (int i1 = 0; i1 < l_SNOVA; ++i1)
				for (int j1 = 0; j1 < rank; ++j1)
					for (int k1 = 0; k1 < l_SNOVA; ++k1)
						whipped_sig[idx * l_SNOVA * SNOVA_lr + i1 * SNOVA_lr + ab * rank + j1] ^=
						    official_gf16_expand(get_gf16m(snova_S[ab], i1, k1)) *
						    signature_in_GF[idx * SNOVA_lr + k1 * rank + j1];

	for (int idx = 0; idx < l_SNOVA * n_SNOVA * SNOVA_lr; ++idx) {
		whipped_sig[idx] = official_gf16_cleanup(whipped_sig[idx]);
	}

	uint16_t hash_in_GF[o_SNOVA * SNOVA_lr] = {0};
	uint16_t sum_t1[m_SNOVA * SNOVA_l2 * SNOVA_r2] = {0};
	uint16_t sum_t1s[m_SNOVA * SNOVA_l2 * SNOVA_r2] = {0};

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < n_SNOVA; ++ni) {
			uint16_t sum_t0[l_SNOVA * SNOVA_lr] = {0};

			for (int nj = 0; nj < n_SNOVA; ++nj)
				for (int k1 = 0; k1 < l_SNOVA; ++k1)
					for (int i1 = 0; i1 < l_SNOVA; ++i1)
						for (int b1 = 0; b1 < SNOVA_lr; ++b1)
							sum_t0[i1 * SNOVA_lr + b1] ^=
							    pkx->P[((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2 + i1 * l_SNOVA + k1] *
							    whipped_sig[nj * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + b1];

			for (int idx = 0; idx < l_SNOVA * SNOVA_lr; ++idx) {
				sum_t0[idx] = official_gf16_compress(sum_t0[idx]);
			}

			for (int a1 = 0; a1 < l_SNOVA; ++a1)
				for (int k1 = 0; k1 < l_SNOVA; ++k1)
					for (int i1 = 0; i1 < rank; ++i1)
						for (int b1 = 0; b1 < SNOVA_lr; ++b1)
							sum_t1[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr + b1] ^=
							    whipped_sig[ni * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + a1 * rank + i1] *
							    sum_t0[k1 * SNOVA_lr + b1];
		}
	}

	for (int idx = 0; idx < m_SNOVA * SNOVA_l2 * SNOVA_r2; ++idx) {
		sum_t1[idx] = official_gf16_compress(sum_t1[idx]);
	}

	for (int mi = 0; mi < m_SNOVA; ++mi)
		for (int a1 = 0; a1 < l_SNOVA; ++a1)
			for (int i1 = 0; i1 < rank; ++i1)
				for (int b1 = 0; b1 < l_SNOVA; ++b1)
					for (int j1 = 0; j1 < rank; ++j1)
						sum_t1s[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + b1 * SNOVA_r2 + i1 * rank + j1] =
						    sum_t1[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr + b1 * rank + j1];

	uint16_t Amx[o_SNOVA * alpha_SNOVA * SNOVA_r2];
	uint16_t Bmx[o_SNOVA * alpha_SNOVA * SNOVA_lr];
	uint16_t q1x[o_SNOVA * alpha_SNOVA * l_SNOVA];
	uint16_t q2x[o_SNOVA * alpha_SNOVA * l_SNOVA];

	for (int idx = 0; idx < o_SNOVA * alpha_SNOVA * SNOVA_r2; ++idx) {
		Amx[idx] = official_gf16_expand(pkx->Am[idx]);
	}
	for (int idx = 0; idx < o_SNOVA * alpha_SNOVA * SNOVA_lr; ++idx) {
		Bmx[idx] = official_gf16_expand(pkx->Bm[idx]);
	}
	for (int idx = 0; idx < o_SNOVA * alpha_SNOVA * l_SNOVA; ++idx) {
		q1x[idx] = official_gf16_expand(pkx->q1[idx]);
		q2x[idx] = official_gf16_expand(pkx->q2[idx]);
	}

	for (int mi = 0; mi < o_SNOVA; ++mi)
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime = i_prime(mi, alpha);
			uint16_t gfm_temp1[SNOVA_r2] = {0};
			uint16_t gfm_temp2[SNOVA_lr] = {0};

			for (int a1 = 0; a1 < l_SNOVA; ++a1) {
				uint16_t gfm_temp0[SNOVA_r2] = {0};

				for (int i1 = 0; i1 < rank; ++i1)
					for (int b1 = 0; b1 < l_SNOVA; ++b1)
						for (int j1 = 0; j1 < rank; ++j1)
							gfm_temp0[i1 * rank + j1] ^=
							    q2x[(mi * alpha_SNOVA + alpha) * l_SNOVA + b1] *
							    sum_t1s[(mi_prime * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + b1 * SNOVA_r2 + i1 * rank + j1];

				for (int idx = 0; idx < SNOVA_r2; ++idx) {
					gfm_temp0[idx] = official_gf16_compress(gfm_temp0[idx]);
				}

				for (int idx = 0; idx < SNOVA_r2; ++idx) {
					gfm_temp1[idx] ^= gfm_temp0[idx] * q1x[(mi * alpha_SNOVA + alpha) * l_SNOVA + a1];
				}
			}

			for (int idx = 0; idx < SNOVA_r2; ++idx) {
				gfm_temp1[idx] = official_gf16_compress(gfm_temp1[idx]);
			}

			for (int i1 = 0; i1 < rank; ++i1)
				for (int j1 = 0; j1 < l_SNOVA; ++j1)
					for (int k1 = 0; k1 < rank; ++k1)
						gfm_temp2[i1 * l_SNOVA + j1] ^=
						    gfm_temp1[i1 * rank + k1] *
						    Bmx[(mi * alpha_SNOVA + alpha) * SNOVA_lr + k1 * l_SNOVA + j1];

			for (int idx = 0; idx < SNOVA_lr; ++idx) {
				gfm_temp2[idx] = official_gf16_compress(gfm_temp2[idx]);
			}

			for (int i1 = 0; i1 < rank; ++i1)
				for (int j1 = 0; j1 < l_SNOVA; ++j1)
					for (int k1 = 0; k1 < rank; ++k1)
						hash_in_GF[mi * SNOVA_lr + i1 * l_SNOVA + j1] ^=
						    Amx[(mi * alpha_SNOVA + alpha) * SNOVA_r2 + i1 * rank + k1] *
						    gfm_temp2[k1 * l_SNOVA + j1];
		}

	for (int idx = 0; idx < o_SNOVA * SNOVA_lr; ++idx) {
		hash_in_GF[idx] = official_gf16_compress(hash_in_GF[idx]);
	}

	int result = 0;
	for (int idx = 0; idx < GF16_HASH; ++idx) {
		if (hash_in_GF[idx] != signed_gf[idx]) {
			result = -1;
			break;
		}
	}
	return result;
}

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#if SNOVA_VERIFY_EXPANDED_NEON
static int snova_verify_official_after_pk_expand_neon(const expanded_PK* pkx, const uint8_t* sig,
                                                      const uint8_t* digest, const size_t len_digest) {
	gf16_t signature_in_GF[n_SNOVA * SNOVA_lr];
	SNOVA_PROFILE_BEGIN(verify_unpack_t);
	if (snova_bytes_to_gf16(sig, signature_in_GF, NUMGF_SIGNATURE)) {
		SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);
		return -1;
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);

	const uint8_t* salt = sig + bytes_signature;
	uint8_t signed_hash[bytes_hash];
	uint8_t signed_gf[GF16_HASH] = {0};
	SNOVA_PROFILE_BEGIN(verify_hash_t);
	snova_signed_hash(digest, len_digest, pkx->pk_seed, salt, signed_hash);
#if (o_SNOVA * l_SNOVA) & 0x1 == 1
	signed_hash[bytes_hash - 1] &= 0x0f;
#endif
	snova_bytes_to_gf16(signed_hash, signed_gf, GF16_HASH);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_HASH, verify_hash_t);

	/* The whip stage keeps the uint16 representation. */
	uint16_t whipped_sig[l_SNOVA * n_SNOVA * SNOVA_lr] = {0};
	SNOVA_PROFILE_BEGIN(verify_expanded_whip_t);
	for (int ab = 0; ab < l_SNOVA; ++ab) {
		for (int idx = 0; idx < n_SNOVA; ++idx) {
			for (int i1 = 0; i1 < l_SNOVA; ++i1) {
				for (int j1 = 0; j1 < rank; ++j1) {
					uint16_t acc = 0;
					for (int k1 = 0; k1 < l_SNOVA; ++k1) {
						acc ^= official_gf16_expand(get_gf16m(snova_S[ab], i1, k1)) *
						        signature_in_GF[idx * SNOVA_lr + k1 * rank + j1];
					}
					whipped_sig[idx * l_SNOVA * SNOVA_lr + i1 * SNOVA_lr + ab * rank + j1] = acc;
				}
			}
		}
	}
	for (size_t idx = 0; idx < sizeof(whipped_sig) / sizeof(whipped_sig[0]); ++idx) {
		whipped_sig[idx] = official_gf16_cleanup(whipped_sig[idx]);
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EXPANDED_WHIP, verify_expanded_whip_t);

	uint16_t hash_in_GF[o_SNOVA * SNOVA_lr] = {0};
	uint16_t sum_t1[m_SNOVA * SNOVA_l2 * SNOVA_r2] = {0};
	uint16_t sum_t1s[m_SNOVA * SNOVA_l2 * SNOVA_r2] = {0};

	SNOVA_PROFILE_BEGIN(verify_expanded_contract_t);
	/* Right contraction: P * whipped_sig, vectorized across the lr lanes. */
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < n_SNOVA; ++ni) {
			uint16_t sum_t0[l_SNOVA * SNOVA_lr] = {0};
			for (int nj = 0; nj < n_SNOVA; ++nj) {
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					const uint16_t* src = whipped_sig + nj * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr;
#if SNOVA_VERIFY_EXPANDED_ILP4 && l_SNOVA == 4
					const uint16_t* p = pkx->P + ((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2 + k1;
					official_const4_expanded_xor_array_neon(
						sum_t0, sum_t0 + SNOVA_lr, sum_t0 + 2 * SNOVA_lr, sum_t0 + 3 * SNOVA_lr,
						p[0 * l_SNOVA], p[1 * l_SNOVA], p[2 * l_SNOVA], p[3 * l_SNOVA],
						src, SNOVA_lr);
#else
					for (int i1 = 0; i1 < l_SNOVA; ++i1) {
						const uint16_t p = pkx->P[((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2 +
						                         i1 * l_SNOVA + k1];
						official_const_expanded_xor_array_neon(sum_t0 + i1 * SNOVA_lr, p, src, SNOVA_lr);
					}
#endif
				}
			}
			official_compress_array_neon(sum_t0, l_SNOVA * SNOVA_lr);

			/* Left contraction: whipped_sig^T * sum_t0, also vectorized over lr. */
#if SNOVA_VERIFY_EXPANDED_ILP4 && l_SNOVA == 4
			for (int k1 = 0; k1 < l_SNOVA; ++k1) {
				const uint16_t* src = sum_t0 + k1 * SNOVA_lr;
				for (int i1 = 0; i1 < rank; ++i1) {
					const uint16_t* s = whipped_sig + ni * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + i1;
					official_const4_expanded_xor_array_neon(
						sum_t1 + (mi * l_SNOVA + 0) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr,
						sum_t1 + (mi * l_SNOVA + 1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr,
						sum_t1 + (mi * l_SNOVA + 2) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr,
						sum_t1 + (mi * l_SNOVA + 3) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr,
						s[0 * rank], s[1 * rank], s[2 * rank], s[3 * rank], src, SNOVA_lr);
				}
			}
#else
			for (int a1 = 0; a1 < l_SNOVA; ++a1) {
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					for (int i1 = 0; i1 < rank; ++i1) {
						const uint16_t s = whipped_sig[ni * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + a1 * rank + i1];
						uint16_t* dst = sum_t1 + (mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr;
						official_const_expanded_xor_array_neon(dst, s, sum_t0 + k1 * SNOVA_lr, SNOVA_lr);
					}
				}
			}
#endif
		}
	}
	official_compress_array_neon(sum_t1, m_SNOVA * SNOVA_l2 * SNOVA_r2);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EXPANDED_CONTRACT, verify_expanded_contract_t);

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int a1 = 0; a1 < l_SNOVA; ++a1) {
			for (int i1 = 0; i1 < rank; ++i1) {
				for (int b1 = 0; b1 < l_SNOVA; ++b1) {
					for (int j1 = 0; j1 < rank; ++j1) {
						sum_t1s[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + b1 * SNOVA_r2 + i1 * rank + j1] =
							sum_t1[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 + i1 * SNOVA_lr + b1 * rank + j1];
					}
				}
			}
		}
	}

	uint16_t Amx[o_SNOVA * alpha_SNOVA * SNOVA_r2];
	uint16_t Bmx[o_SNOVA * alpha_SNOVA * SNOVA_lr];
	uint16_t q1x[o_SNOVA * alpha_SNOVA * l_SNOVA];
	uint16_t q2x[o_SNOVA * alpha_SNOVA * l_SNOVA];
	SNOVA_PROFILE_BEGIN(verify_expanded_secondary_t);
	official_expand_array_neon(Amx, pkx->Am, sizeof(Amx) / sizeof(Amx[0]));
	official_expand_array_neon(Bmx, pkx->Bm, sizeof(Bmx) / sizeof(Bmx[0]));
	official_expand_array_neon(q1x, pkx->q1, sizeof(q1x) / sizeof(q1x[0]));
	official_expand_array_neon(q2x, pkx->q2, sizeof(q2x) / sizeof(q2x[0]));

	for (int mi = 0; mi < o_SNOVA; ++mi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime = i_prime(mi, alpha);
			uint16_t gfm_temp1[SNOVA_r2] = {0};
			uint16_t gfm_temp2[SNOVA_lr] = {0};
			for (int a1 = 0; a1 < l_SNOVA; ++a1) {
				uint16_t gfm_temp0[SNOVA_r2] = {0};
				for (int i1 = 0; i1 < rank; ++i1) {
					for (int b1 = 0; b1 < l_SNOVA; ++b1) {
						const uint16_t q2 = q2x[(mi * alpha_SNOVA + alpha) * l_SNOVA + b1];
						for (int j1 = 0; j1 < rank; ++j1) {
							gfm_temp0[i1 * rank + j1] ^=
								q2 * sum_t1s[(mi_prime * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 +
								             b1 * SNOVA_r2 + i1 * rank + j1];
						}
					}
				}
				official_compress_array_neon(gfm_temp0, SNOVA_r2);
				const uint16_t q1 = q1x[(mi * alpha_SNOVA + alpha) * l_SNOVA + a1];
				for (int idx = 0; idx < SNOVA_r2; ++idx) {
					gfm_temp1[idx] ^= gfm_temp0[idx] * q1;
				}
			}
			official_compress_array_neon(gfm_temp1, SNOVA_r2);

			for (int i1 = 0; i1 < rank; ++i1) {
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					for (int k1 = 0; k1 < rank; ++k1) {
						gfm_temp2[i1 * l_SNOVA + j1] ^=
							gfm_temp1[i1 * rank + k1] *
							Bmx[(mi * alpha_SNOVA + alpha) * SNOVA_lr + k1 * l_SNOVA + j1];
					}
				}
			}
			official_compress_array_neon(gfm_temp2, SNOVA_lr);
			for (int i1 = 0; i1 < rank; ++i1) {
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					for (int k1 = 0; k1 < rank; ++k1) {
						hash_in_GF[mi * SNOVA_lr + i1 * l_SNOVA + j1] ^=
							Amx[(mi * alpha_SNOVA + alpha) * SNOVA_r2 + i1 * rank + k1] *
							gfm_temp2[k1 * l_SNOVA + j1];
					}
				}
			}
		}
	}
	official_compress_array_neon(hash_in_GF, o_SNOVA * SNOVA_lr);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EXPANDED_SECONDARY, verify_expanded_secondary_t);

	SNOVA_PROFILE_BEGIN(verify_compare_t);
	for (int idx = 0; idx < GF16_HASH; ++idx) {
		if ((uint8_t)hash_in_GF[idx] != signed_gf[idx]) {
			SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_COMPARE, verify_compare_t);
			return -1;
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_COMPARE, verify_compare_t);
	return 0;
}
#endif
#endif

static int snova_verify_official_after_pk_expand(const expanded_PK* pkx, const uint8_t* sig,
                                                   const uint8_t* digest, const size_t len_digest) {
#if SNOVA_VERIFY_EXPANDED_NEON && (defined(__ARM_NEON) || defined(__ARM_NEON__))
	return snova_verify_official_after_pk_expand_neon(pkx, sig, digest, len_digest);
#else
	return snova_verify_official_after_pk_expand_scalar(pkx, sig, digest, len_digest);
#endif
}
#endif

#if SNOVA_VERIFY_DIRECT_MTK2
#if rank != 4 || l_SNOVA != 4 || SNOVA_lr != 16
#error "SNOVA_VERIFY_DIRECT_MTK2 currently targets rank=4, l=4"
#endif
static inline VECTOR direct_mtk2_load16(const uint8_t* src) {
#if defined(__AVX2__)
	return _mm256_inserti128_si256(_mm256_setzero_si256(), _mm_loadu_si128((const __m128i*)src), 0);
#else
	return vld1q_u8(src);
#endif
}

#if defined(__AVX2__)
static inline VECTOR direct_mtk2_load16_pair(const uint8_t* lo, const uint8_t* hi) {
	return _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i*)lo)),
	                               _mm_loadu_si128((const __m128i*)hi), 1);
}

static inline VECTOR direct_mtk2_table_pair(uint8_t lo, uint8_t hi) {
	const __m128i tab_lo = _mm_load_si128((const __m128i*)(mtk2_16 + lo));
	const __m128i tab_hi = _mm_load_si128((const __m128i*)(mtk2_16 + hi));
	return _mm256_inserti128_si256(_mm256_castsi128_si256(tab_lo), tab_hi, 1);
}

static inline VECTOR direct_mtk2_fold_halves(VECTOR value) {
	const __m128i folded = _mm_xor_si128(_mm256_castsi256_si128(value),
	                                      _mm256_extracti128_si256(value, 1));
	return _mm256_inserti128_si256(_mm256_setzero_si256(), folded, 0);
}
#endif

static inline void direct_mtk2_store16(uint8_t* dst, VECTOR value) {
#if defined(__AVX2__)
	_mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(value));
#else
	vst1q_u8(dst, value);
#endif
}

static inline void direct_mtk2_unpack_gf16_run(uint8_t* restrict dst, const uint8_t* restrict src,
                                               size_t count) {
	size_t idx = 0;
#if defined(__AVX2__)
	const __m128i low_mask = _mm_set1_epi8(0x0f);
	for (; idx + 32 <= count; idx += 32) {
		const __m128i packed = _mm_loadu_si128((const __m128i*)(src + idx / 2));
		const __m128i low = _mm_and_si128(packed, low_mask);
		const __m128i high = _mm_and_si128(_mm_srli_epi16(packed, 4), low_mask);
		_mm_storeu_si128((__m128i*)(dst + idx), _mm_unpacklo_epi8(low, high));
		_mm_storeu_si128((__m128i*)(dst + idx + 16), _mm_unpackhi_epi8(low, high));
	}
#endif
	for (; idx < count; ++idx) {
		const uint8_t packed = src[idx / 2];
		dst[idx] = (idx & 1) ? (packed >> 4) : (packed & 0x0f);
	}
}

static inline void direct_mtk2_expand_array(uint16_t* dst, const uint8_t* src, size_t count) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	official_expand_array_neon(dst, src, count);
#else
	for (size_t i = 0; i < count; ++i) {
		dst[i] = official_gf16_expand(src[i]);
	}
#endif
}

static inline void direct_mtk2_compress_array(uint16_t* dst, size_t count) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	official_compress_array_neon(dst, count);
#else
	for (size_t i = 0; i < count; ++i) {
		dst[i] = official_gf16_compress(dst[i]);
	}
#endif
}

static int snova_verify_direct_mtk2_after_pk_expand(const uint8_t* pkx, const uint8_t* sig,
                                                     const uint8_t* digest, const size_t len_digest) {
	const uint8_t* const pk_P = pkx + DIRECT_MTK2_P_OFFSET;
	const uint8_t* const pk_Am = pkx + DIRECT_MTK2_AM_OFFSET;
	const uint8_t* const pk_Bm = pkx + DIRECT_MTK2_BM_OFFSET;
	const uint8_t* const pk_q1 = pkx + DIRECT_MTK2_Q1_OFFSET;
	const uint8_t* const pk_q2 = pkx + DIRECT_MTK2_Q2_OFFSET;
	const uint8_t* const pk_seed = pkx + DIRECT_MTK2_PK_SEED_OFFSET;
	gf16_t signature_in_GF[n_SNOVA * SNOVA_lr];
	SNOVA_PROFILE_BEGIN(verify_unpack_t);
	if (snova_bytes_to_gf16(sig, signature_in_GF, NUMGF_SIGNATURE)) {
		SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);
		return -1;
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);

	const uint8_t* salt = sig + bytes_signature;
	uint8_t signed_hash[bytes_hash];
	uint8_t signed_gf[GF16_HASH] = {0};
	SNOVA_PROFILE_BEGIN(verify_hash_t);
	snova_signed_hash(digest, len_digest, pk_seed, salt, signed_hash);
#if (o_SNOVA * l_SNOVA) & 0x1 == 1
	signed_hash[bytes_hash - 1] &= 0x0f;
#endif
	snova_bytes_to_gf16(signed_hash, signed_gf, GF16_HASH);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_HASH, verify_hash_t);

	uint8_t whipped_sig[n_SNOVA * l_SNOVA * SNOVA_lr] = {0};
	SNOVA_PROFILE_BEGIN(verify_direct_mtk2_whip_t);
	for (int ab = 0; ab < l_SNOVA; ++ab) {
		for (int idx = 0; idx < n_SNOVA; ++idx) {
			for (int i1 = 0; i1 < l_SNOVA; ++i1) {
				for (int j1 = 0; j1 < rank; ++j1) {
					uint8_t acc = 0;
					for (int k1 = 0; k1 < l_SNOVA; ++k1) {
						acc ^= mt(get_gf16m(snova_S[ab], i1, k1),
						          signature_in_GF[idx * SNOVA_lr + k1 * rank + j1]);
					}
					whipped_sig[idx * l_SNOVA * SNOVA_lr + i1 * SNOVA_lr + ab * rank + j1] = acc;
				}
			}
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_DIRECT_MTK2_WHIP, verify_direct_mtk2_whip_t);

	VECTOR sum_t1_pair[m_SNOVA][rank_next32][rank];
	memset(sum_t1_pair, 0, sizeof(sum_t1_pair));
	SNOVA_PROFILE_BEGIN(verify_direct_mtk2_contract_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < n_SNOVA; ++ni) {
			VECTOR sum_t0_01 = zero256;
			VECTOR sum_t0_23 = zero256;
#if defined(__AVX2__)
			int nj = 0;
			for (; nj + 1 < n_SNOVA; nj += 2) {
				const uint8_t* const p0 = pk_P + ((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2;
				const uint8_t* const p1 = p0 + SNOVA_l2;
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					const VECTOR src = direct_mtk2_load16_pair(
						whipped_sig + nj * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr,
						whipped_sig + (nj + 1) * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr);
					const VECTOR tab01 = direct_mtk2_table_pair(
						p0[0 * l_SNOVA + k1] | (p0[1 * l_SNOVA + k1] << 4),
						p1[0 * l_SNOVA + k1] | (p1[1 * l_SNOVA + k1] << 4));
					const VECTOR tab23 = direct_mtk2_table_pair(
						p0[2 * l_SNOVA + k1] | (p0[3 * l_SNOVA + k1] << 4),
						p1[2 * l_SNOVA + k1] | (p1[3 * l_SNOVA + k1] << 4));
					sum_t0_01 ^= VEC_SHUFFLE(tab01, src);
					sum_t0_23 ^= VEC_SHUFFLE(tab23, src);
				}
			}
			if (nj < n_SNOVA) {
				const uint8_t* const p = pk_P + ((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2;
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					const VECTOR src = direct_mtk2_load16(
						whipped_sig + nj * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr);
					sum_t0_01 ^= VEC_SHUFFLE(mtk2_16[p[0 * l_SNOVA + k1] |
					                                           (p[1 * l_SNOVA + k1] << 4)], src);
					sum_t0_23 ^= VEC_SHUFFLE(mtk2_16[p[2 * l_SNOVA + k1] |
					                                           (p[3 * l_SNOVA + k1] << 4)], src);
				}
			}
			sum_t0_01 = direct_mtk2_fold_halves(sum_t0_01);
			sum_t0_23 = direct_mtk2_fold_halves(sum_t0_23);
#else
			for (int nj = 0; nj < n_SNOVA; ++nj) {
				const uint8_t* p = pk_P + ((mi * n_SNOVA + ni) * n_SNOVA + nj) * SNOVA_l2;
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					const VECTOR src = direct_mtk2_load16(whipped_sig + nj * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr);
					const VECTOR tab01 = mtk2_16[p[0 * l_SNOVA + k1] | (p[1 * l_SNOVA + k1] << 4)];
					const VECTOR tab23 = mtk2_16[p[2 * l_SNOVA + k1] | (p[3 * l_SNOVA + k1] << 4)];
					sum_t0_01 ^= VEC_SHUFFLE(tab01, src);
					sum_t0_23 ^= VEC_SHUFFLE(tab23, src);
				}
			}
#endif

			for (int k1 = 0; k1 < l_SNOVA; ++k1) {
				const VECTOR pair = k1 < 2 ? sum_t0_01 : sum_t0_23;
				const VECTOR src = (k1 & 1) ? ((pair >> 4) & l_mask) : (pair & l_mask);
				for (int i1 = 0; i1 < rank; ++i1) {
					const uint8_t* s = whipped_sig + ni * l_SNOVA * SNOVA_lr + k1 * SNOVA_lr + i1;
					const VECTOR tab01 = mtk2_16[s[0 * rank] | (s[1 * rank] << 4)];
					const VECTOR tab23 = mtk2_16[s[2 * rank] | (s[3 * rank] << 4)];
					sum_t1_pair[mi][0][i1] ^= VEC_SHUFFLE(tab01, src);
					sum_t1_pair[mi][1][i1] ^= VEC_SHUFFLE(tab23, src);
				}
			}
		}
	}

	uint16_t sum_t1s[m_SNOVA * SNOVA_l2 * SNOVA_r2];
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int pair = 0; pair < rank_next32; ++pair) {
			for (int i1 = 0; i1 < rank; ++i1) {
				uint8_t lanes[SNOVA_lr];
				direct_mtk2_store16(lanes, sum_t1_pair[mi][pair][i1]);
				for (int half = 0; half < 2; ++half) {
					const int a1 = 2 * pair + half;
					for (int b1 = 0; b1 < l_SNOVA; ++b1) {
						for (int j1 = 0; j1 < rank; ++j1) {
							const uint8_t packed = lanes[b1 * rank + j1];
							const uint8_t raw = half ? (packed >> 4) : (packed & 0x0f);
							sum_t1s[(mi * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 +
							          b1 * SNOVA_r2 + i1 * rank + j1] = raw;
						}
					}
				}
			}
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_DIRECT_MTK2_CONTRACT, verify_direct_mtk2_contract_t);

	uint16_t hash_in_GF[o_SNOVA * SNOVA_lr] = {0};
	uint16_t Amx[o_SNOVA * alpha_SNOVA * SNOVA_r2];
	uint16_t Bmx[o_SNOVA * alpha_SNOVA * SNOVA_lr];
	uint16_t q1x[o_SNOVA * alpha_SNOVA * l_SNOVA];
	uint16_t q2x[o_SNOVA * alpha_SNOVA * l_SNOVA];
	SNOVA_PROFILE_BEGIN(verify_direct_mtk2_secondary_t);
	direct_mtk2_expand_array(Amx, pk_Am, sizeof(Amx) / sizeof(Amx[0]));
	direct_mtk2_expand_array(Bmx, pk_Bm, sizeof(Bmx) / sizeof(Bmx[0]));
	direct_mtk2_expand_array(q1x, pk_q1, sizeof(q1x) / sizeof(q1x[0]));
	direct_mtk2_expand_array(q2x, pk_q2, sizeof(q2x) / sizeof(q2x[0]));
	for (int mi = 0; mi < o_SNOVA; ++mi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime = i_prime(mi, alpha);
			uint16_t gfm_temp1[SNOVA_r2] = {0};
			uint16_t gfm_temp2[SNOVA_lr] = {0};
			for (int a1 = 0; a1 < l_SNOVA; ++a1) {
				uint16_t gfm_temp0[SNOVA_r2] = {0};
				for (int i1 = 0; i1 < rank; ++i1) {
					for (int b1 = 0; b1 < l_SNOVA; ++b1) {
						const uint16_t q2 = q2x[(mi * alpha_SNOVA + alpha) * l_SNOVA + b1];
						for (int j1 = 0; j1 < rank; ++j1) {
							gfm_temp0[i1 * rank + j1] ^=
								q2 * sum_t1s[(mi_prime * l_SNOVA + a1) * l_SNOVA * SNOVA_r2 +
								               b1 * SNOVA_r2 + i1 * rank + j1];
						}
					}
				}
				direct_mtk2_compress_array(gfm_temp0, SNOVA_r2);
				const uint16_t q1 = q1x[(mi * alpha_SNOVA + alpha) * l_SNOVA + a1];
				for (int idx = 0; idx < SNOVA_r2; ++idx) {
					gfm_temp1[idx] ^= gfm_temp0[idx] * q1;
				}
			}
			direct_mtk2_compress_array(gfm_temp1, SNOVA_r2);
			for (int i1 = 0; i1 < rank; ++i1) {
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					for (int k1 = 0; k1 < rank; ++k1) {
						gfm_temp2[i1 * l_SNOVA + j1] ^=
							gfm_temp1[i1 * rank + k1] *
							Bmx[(mi * alpha_SNOVA + alpha) * SNOVA_lr + k1 * l_SNOVA + j1];
					}
				}
			}
			direct_mtk2_compress_array(gfm_temp2, SNOVA_lr);
			for (int i1 = 0; i1 < rank; ++i1) {
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					for (int k1 = 0; k1 < rank; ++k1) {
						hash_in_GF[mi * SNOVA_lr + i1 * l_SNOVA + j1] ^=
							Amx[(mi * alpha_SNOVA + alpha) * SNOVA_r2 + i1 * rank + k1] *
							gfm_temp2[k1 * l_SNOVA + j1];
					}
				}
			}
		}
	}
	direct_mtk2_compress_array(hash_in_GF, o_SNOVA * SNOVA_lr);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_DIRECT_MTK2_SECONDARY, verify_direct_mtk2_secondary_t);

	SNOVA_PROFILE_BEGIN(verify_compare_t);
	for (int idx = 0; idx < GF16_HASH; ++idx) {
		if ((uint8_t)hash_in_GF[idx] != signed_gf[idx]) {
			SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_COMPARE, verify_compare_t);
			return -1;
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_COMPARE, verify_compare_t);
	return 0;
}
#endif

static void evaluation_vector_vtl(gf16m_t* restrict hash_in_GF16Matrix, const public_key_expand* restrict pkx,
                                  gf16m_t* restrict signature_in_GF16Matrix) {
	VECTOR S_J_256[rank][vtl_len] = {0};
	VECTOR S_tr_J_256[rank][vtl_len] = {0};
	alignas(VEC_LENGTH) verify_pair_t L_J_pair;
	alignas(VEC_LENGTH) verify_pair_t R_tr_J_pair;
	VECTOR h_256[m_SNOVA][rank][VERIFY_H_COLS] = {0};
	SNOVA_PROFILE_BEGIN(verify_jogress_t);
	snova_pack_blocks((VECTOR*)S_J_256, (const uint8_t*)signature_in_GF16Matrix, 1, n_SNOVA);
	snova_pack_blocks_tr((VECTOR*)S_tr_J_256, (const uint8_t*)signature_in_GF16Matrix, 1, n_SNOVA);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_JOGRESS, verify_jogress_t);
	SNOVA_PROFILE_BEGIN(verify_lr_t);
	verify_eval_lr(pkx, S_J_256, S_tr_J_256, L_J_pair, R_tr_J_pair);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_LR, verify_lr_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		VECTOR P_J_256_local[n_SNOVA * rank][vtl_len];
		VECTOR (*P_J_256)[vtl_len] = P_J_256_local;
		SNOVA_PROFILE_BEGIN(verify_p_layout_t);
		snova_fill_public_layout_mi(P_J_256_local, pkx, mi);
		SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_P_LAYOUT, verify_p_layout_t);
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const int mi_prime = i_prime_inv(mi, alpha);
			VECTOR LJxPJ_256_nibble[rank_next32][vtl_len];
			SNOVA_PROFILE_BEGIN(verify_ljpj_t);
			verify_eval_ljpj((const uint8_t*)L_J_pair[mi][alpha], (const VECTOR(*)[vtl_len])P_J_256, LJxPJ_256_nibble);
			SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_LJPJ, verify_ljpj_t);
			SNOVA_PROFILE_BEGIN(verify_secondary_t);
			verify_eval_secondary(LJxPJ_256_nibble, (const uint8_t*)R_tr_J_pair[mi][alpha], h_256[mi_prime]);
			SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_SECONDARY, verify_secondary_t);
		}
	}
	SNOVA_PROFILE_BEGIN(verify_reduce_t);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < rank; ++ni) {
#if SNOVA_VERIFY_PACKED_SECONDARY && rank == 4
			for (int p = 0; p < rank_next32; ++p) {
				uint8_t folded = horizontal_xor_256(h_256[mi][ni][p]);
				set_gf16m(hash_in_GF16Matrix[mi], ni, 2 * p, folded & 0x0f);
				if (2 * p + 1 < rank) {
					set_gf16m(hash_in_GF16Matrix[mi], ni, 2 * p + 1, folded >> 4);
				}
			}
#else
			for (int nj = 0; nj < rank; ++nj) {
				set_gf16m(hash_in_GF16Matrix[mi], ni, nj, horizontal_xor_256(h_256[mi][ni][nj]));
			}
#endif
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_EVAL_REDUCE, verify_reduce_t);
}
int snova_verify_core(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature,
                       const public_key_expand* pkx) {
	uint8_t hash_in_bytes[bytes_hash];
	uint8_t signed_hash[bytes_hash];
	const uint8_t* pt_salt = pt_signature + bytes_signature;
	gf16m_t hash_in_GF16Matrix[m_SNOVA];
	gf16m_t signature_in_GF16Matrix[n_SNOVA];
	SNOVA_PROFILE_BEGIN(verify_hash_t);
	snova_signed_hash(pt_digest, bytes_digest, pkx->pt_public_key_seed, pt_salt, signed_hash);
#if (o_SNOVA * l_SNOVA) & 0x1 == 1
	signed_hash[bytes_hash - 1] &= 0x0f;
#endif
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_HASH, verify_hash_t);
	SNOVA_PROFILE_BEGIN(verify_unpack_t);
	if (snova_bytes_to_gf16(pt_signature, (gf16_t*)signature_in_GF16Matrix, GF16s_signature)) {
		SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);
		return -1;
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_UNPACK, verify_unpack_t);
#if SNOVA_VERIFY_DIRECT_OFFICIAL && !FIXED_ABQ
	evaluation_vector_direct_official(hash_in_GF16Matrix, pkx, signature_in_GF16Matrix);
#else
	evaluation_vector_vtl(hash_in_GF16Matrix, pkx, signature_in_GF16Matrix);
#endif
	SNOVA_PROFILE_BEGIN(verify_compress_t);
	snova_gf16_to_bytes(hash_in_bytes, (const gf16_t*)hash_in_GF16Matrix, m_SNOVA * lsq_SNOVA);
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_COMPRESS, verify_compress_t);
	int result = 0;
	SNOVA_PROFILE_BEGIN(verify_compare_t);
	for (int i = 0; i < bytes_hash; ++i) {
		if (hash_in_bytes[i] != signed_hash[i]) {
			result = -1;
			break;
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_VERIFY_COMPARE, verify_compare_t);
	return result;
}

#if SNOVA_q != 16
#error "SNOVA_q != 16"
#endif
static void snova_init(void) {
	static int first_time = 1;
	if (first_time) {
		first_time = 0;
		snova_init_gf16_tables();
		snova_gen_S_array();
#if FIXED_ABQ
		snova_gen_fixed_ABQ("SNOVA_ABQ");
#endif
		snova_init_vector_table();
	}
}
static void generate_secret_key_core(snova_key_elems* key_elems, const uint8_t* pk_seed, const uint8_t* sk_seed) {
	SNOVA_PROFILE_BEGIN(genkeys_map_expand_t);
	snova_gen_seeds_and_T12(key_elems->T12, sk_seed);
	memcpy(key_elems->pk.pt_public_key_seed, pk_seed, seed_length_public);
	snova_gen_ABQ(&(key_elems->map1), pk_seed);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_MAP_EXPAND, genkeys_map_expand_t);
	SNOVA_PROFILE_BEGIN(genkeys_gen_f_t);
	snova_gen_F_vtl(&(key_elems->map2), &(key_elems->map1), key_elems->T12);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_GEN_F, genkeys_gen_f_t);
}
/* Keypair only needs P11/P12/P21; ABQ is regenerated for sign and verify. */
static void generate_keypair_core(snova_key_elems* key_elems, const uint8_t* pk_seed, const uint8_t* sk_seed) {
	SNOVA_PROFILE_BEGIN(genkeys_map_expand_t);
	snova_gen_seeds_and_T12(key_elems->T12, sk_seed);
	memcpy(key_elems->pk.pt_public_key_seed, pk_seed, seed_length_public);
	snova_gen_public_P_only(&(key_elems->map1), pk_seed);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_MAP_EXPAND, genkeys_map_expand_t);
	SNOVA_PROFILE_BEGIN(genkeys_gen_f_t);
	snova_gen_F_vtl(&(key_elems->map2), &(key_elems->map1), key_elems->T12);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_GEN_F, genkeys_gen_f_t);
}
int SNOVA_NAMESPACE(genkeys)(uint8_t* pk, uint8_t* ssk, const uint8_t* seed) {
	const uint8_t* pkseed = seed;
	const uint8_t* skseed = seed + seed_length_public;
	snova_key_elems key_elems;
	snova_init();
	generate_keypair_core(&key_elems, pkseed, skseed);
	snova_gen_P22_vtl(key_elems.pk.P22, key_elems.T12, key_elems.map1.P21, key_elems.map2.F12);
	SNOVA_PROFILE_BEGIN(genkeys_pk_pack_t);
	snova_pk_pack(pk, &key_elems);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_PK_PACK, genkeys_pk_pack_t);
	SNOVA_PROFILE_BEGIN(genkeys_sk_pack_t);
	memcpy(ssk, pkseed, seed_length_public);
	memcpy(ssk + seed_length_public, skseed, seed_length_private);
	SNOVA_PROFILE_END(SNOVA_PROFILE_GENKEYS_SK_PACK, genkeys_sk_pack_t);
	return 0;
}
int SNOVA_NAMESPACE(sk_expand)(expanded_SK* skx, const uint8_t* sk) {
	const uint8_t* pk_seed = sk;
	const uint8_t* sk_seed = sk + seed_length_public;
	snova_key_elems key_elems;
	snova_init();
	SNOVA_PROFILE_BEGIN(sk_expand_core_t);
	generate_secret_key_core(&key_elems, pk_seed, sk_seed);
	SNOVA_PROFILE_END(SNOVA_PROFILE_SK_EXPAND_CORE, sk_expand_core_t);
	SNOVA_PROFILE_BEGIN(sk_expand_pack_t);
	snova_sk_pack((uint8_t*)skx, &key_elems, sk_seed);
	SNOVA_PROFILE_END(SNOVA_PROFILE_SK_EXPAND_PACK, sk_expand_pack_t);
	return 0;
}
int SNOVA_NAMESPACE(sign)(const expanded_SK* esk, uint8_t* pt_signature, const uint8_t* digest, const size_t bytes_digest,
                          const uint8_t* array_salt) {
	sk_gf16 sk_upk;
	snova_init();
	snova_sk_unpack(&sk_upk, (const uint8_t*)esk);
	return snova_sign_core(pt_signature, digest, bytes_digest, array_salt, &sk_upk);
}
int SNOVA_NAMESPACE(pk_expand)(expanded_PK* pkx, const uint8_t* pk) {
	const public_key* pk_stru = (const public_key*)pk;
	snova_init();
#if SNOVA_VERIFY_DIRECT_MTK2 && !FIXED_ABQ
	uint8_t* const direct = (uint8_t*)(void*)pkx;
	uint8_t prng_output_public[bytes_prng_public];
	const size_t p11_offset = 0;
	const size_t p12_offset = p11_offset + sizeof(P11_t);
	const size_t p21_offset = p12_offset + sizeof(P12_t);
	const size_t a_offset = p21_offset + sizeof(P21_t);
	const size_t b_offset = a_offset + sizeof(Aalpha_t);
	const size_t q1_offset = b_offset + sizeof(Balpha_t);
	const size_t q2_offset = q1_offset + sizeof(Qalpha_coeff_t);

	SNOVA_PROFILE_BEGIN(pk_expand_map_t);
	snova_pk_expand_seed(pk_stru->pt_public_key_seed, prng_output_public);
	direct_mtk2_unpack_gf16_run(direct + DIRECT_MTK2_AM_OFFSET,
	                             prng_output_public + a_offset / 2, DIRECT_MTK2_AM_BYTES);
	direct_mtk2_unpack_gf16_run(direct + DIRECT_MTK2_BM_OFFSET,
	                             prng_output_public + b_offset / 2, DIRECT_MTK2_BM_BYTES);
	direct_mtk2_unpack_gf16_run(direct + DIRECT_MTK2_Q1_OFFSET,
	                             prng_output_public + q1_offset / 2, DIRECT_MTK2_Q1_BYTES);
	direct_mtk2_unpack_gf16_run(direct + DIRECT_MTK2_Q2_OFFSET,
	                             prng_output_public + q2_offset / 2, DIRECT_MTK2_Q2_BYTES);
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			const size_t matrix = (size_t)mi * alpha_SNOVA + alpha;
			snova_be_invertible_by_add_aS(direct + DIRECT_MTK2_AM_OFFSET + matrix * SNOVA_r2);
			snova_be_invertible_by_add_aS(direct + DIRECT_MTK2_BM_OFFSET + matrix * SNOVA_lr);
			gf16_t* const q1 = direct + DIRECT_MTK2_Q1_OFFSET + matrix * l_SNOVA;
			gf16_t* const q2 = direct + DIRECT_MTK2_Q2_OFFSET + matrix * l_SNOVA;
			q1[l_SNOVA - 1] = snova_adjust_q_coeff(q1, l_SNOVA - 1);
			q2[l_SNOVA - 1] = snova_adjust_q_coeff(q2, l_SNOVA - 1);
		}
	}
	SNOVA_PROFILE_END(SNOVA_PROFILE_PK_EXPAND_MAP, pk_expand_map_t);

	SNOVA_PROFILE_BEGIN(pk_expand_p22_t);
	uint8_t* const direct_P = direct + DIRECT_MTK2_P_OFFSET;
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < v_SNOVA; ++ni) {
			uint8_t* const dst = direct_P + ((size_t)(mi * n_SNOVA + ni) * n_SNOVA) * lsq_SNOVA;
			const size_t p11_row = p11_offset + ((size_t)mi * v_SNOVA + ni) * v_SNOVA * lsq_SNOVA;
			const size_t p12_row = p12_offset + ((size_t)mi * v_SNOVA + ni) * o_SNOVA * lsq_SNOVA;
			direct_mtk2_unpack_gf16_run(dst, prng_output_public + p11_row / 2,
			                               v_SNOVA * lsq_SNOVA);
			direct_mtk2_unpack_gf16_run(dst + v_SNOVA * lsq_SNOVA,
			                               prng_output_public + p12_row / 2, o_SNOVA * lsq_SNOVA);
		}
		for (int ni = 0; ni < o_SNOVA; ++ni) {
			uint8_t* const dst = direct_P + ((size_t)(mi * n_SNOVA + v_SNOVA + ni) * n_SNOVA) * lsq_SNOVA;
			const size_t p21_row = p21_offset + ((size_t)mi * o_SNOVA + ni) * v_SNOVA * lsq_SNOVA;
			const size_t p22_row = ((size_t)mi * o_SNOVA + ni) * o_SNOVA * lsq_SNOVA;
			direct_mtk2_unpack_gf16_run(dst, prng_output_public + p21_row / 2,
			                               v_SNOVA * lsq_SNOVA);
			direct_mtk2_unpack_gf16_run(dst + v_SNOVA * lsq_SNOVA,
			                               (const uint8_t*)pk_stru->P22 + p22_row / 2,
			                               o_SNOVA * lsq_SNOVA);
		}
	}
	memcpy(direct + DIRECT_MTK2_PK_SEED_OFFSET, pk_stru->pt_public_key_seed, SEED_LENGTH_PUBLIC);
	SNOVA_PROFILE_END(SNOVA_PROFILE_PK_EXPAND_P22, pk_expand_p22_t);
	return 0;
#else
	public_key_expand tmp;
	memset(pkx, 0, sizeof(*pkx));
	memset(&tmp, 0, sizeof(tmp));
	memcpy(tmp.pt_public_key_seed, pk_stru->pt_public_key_seed, sizeof(pk_stru->pt_public_key_seed));
#if !SNOVA_VERIFY_VTL && !SNOVA_VERIFY_DIRECT_MTK2
	memcpy(pkx->pk_seed, pk_stru->pt_public_key_seed, sizeof(pk_stru->pt_public_key_seed));
#endif
	SNOVA_PROFILE_BEGIN(pk_expand_map_t);
#if !FIXED_ABQ
	snova_gen_ABQ_ex(&(tmp.map1), pk_stru->pt_public_key_seed, tmp.q1coeff, tmp.q2coeff);
#else
	snova_gen_ABQ(&(tmp.map1), pk_stru->pt_public_key_seed);
#endif
	SNOVA_PROFILE_END(SNOVA_PROFILE_PK_EXPAND_MAP, pk_expand_map_t);
	SNOVA_PROFILE_BEGIN(pk_expand_p22_t);
	snova_input_P22((uint8_t*)tmp.P22, (const uint8_t*)pk_stru->P22);
	#if SNOVA_VERIFY_VTL
	/* expanded_PK is an opaque API buffer; keep the VTL-native expansion in it. */
	memcpy(pkx, &tmp, sizeof(tmp));
	#elif SNOVA_VERIFY_DIRECT_MTK2
	uint8_t* const direct = (uint8_t*)(void*)pkx;
	verify_copy_public_blocks_mtk2(direct + DIRECT_MTK2_P_OFFSET, &tmp);
	memcpy(direct + DIRECT_MTK2_AM_OFFSET, tmp.map1.Aalpha, DIRECT_MTK2_AM_BYTES);
	memcpy(direct + DIRECT_MTK2_BM_OFFSET, tmp.map1.Balpha, DIRECT_MTK2_BM_BYTES);
#if !FIXED_ABQ
	memcpy(direct + DIRECT_MTK2_Q1_OFFSET, tmp.q1coeff, DIRECT_MTK2_Q1_BYTES);
	memcpy(direct + DIRECT_MTK2_Q2_OFFSET, tmp.q2coeff, DIRECT_MTK2_Q2_BYTES);
#else
	memset(direct + DIRECT_MTK2_Q1_OFFSET, 0, DIRECT_MTK2_Q1_BYTES);
	memset(direct + DIRECT_MTK2_Q2_OFFSET, 0, DIRECT_MTK2_Q2_BYTES);
#endif
	memcpy(direct + DIRECT_MTK2_PK_SEED_OFFSET, pk_stru->pt_public_key_seed, SEED_LENGTH_PUBLIC);
	#else
	verify_copy_public_blocks(pkx, &tmp);
	memcpy(pkx->Am, tmp.map1.Aalpha, sizeof(pkx->Am));
	memcpy(pkx->Bm, tmp.map1.Balpha, sizeof(pkx->Bm));
#if !FIXED_ABQ
	memcpy(pkx->q1, tmp.q1coeff, sizeof(pkx->q1));
	memcpy(pkx->q2, tmp.q2coeff, sizeof(pkx->q2));
#else
	for (int idx = 0; idx < o_SNOVA * alpha_SNOVA * l_SNOVA; ++idx) {
		pkx->q1[idx] = 0;
		pkx->q2[idx] = 0;
	}
#endif
	#endif
	SNOVA_PROFILE_END(SNOVA_PROFILE_PK_EXPAND_P22, pk_expand_p22_t);
	return 0;
#endif
}
int SNOVA_NAMESPACE(verify)(const expanded_PK* pk, const uint8_t* pt_signature, const uint8_t* pt_digest,
	const size_t bytes_digest) {
	snova_init();
	#if SNOVA_VERIFY_DIRECT_MTK2
	return snova_verify_direct_mtk2_after_pk_expand((const uint8_t*)(const void*)pk,
	                                               pt_signature, pt_digest, bytes_digest);
	#elif !SNOVA_VERIFY_VTL
	return snova_verify_official_after_pk_expand(pk, pt_signature, pt_digest, bytes_digest);
	#else
	return snova_verify_core(pt_digest, bytes_digest, pt_signature,
	                         (const public_key_expand*)(const void*)pk);
#endif
}

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk) {
	uint8_t seed[SEED_LENGTH];
	SNOVA_PROFILE_BEGIN(api_keypair_rng_t);
	randombytes(seed, SEED_LENGTH);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_KEYPAIR_RNG, api_keypair_rng_t);
#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_UNDEFINED(seed, SEED_LENGTH);
#endif
	return SNOVA_NAMESPACE(genkeys)(pk, sk, seed);
}
int crypto_sign(unsigned char* sm, unsigned long long* smlen, const unsigned char* m, unsigned long long mlen,
                const unsigned char* sk) {
	uint8_t digest[BYTES_DIGEST];
	uint8_t salt[BYTES_SALT];
	expanded_SK skx;
	SNOVA_PROFILE_BEGIN(api_sign_sk_expand_t);
	int res = SNOVA_NAMESPACE(sk_expand)(&skx, sk);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_SIGN_SK_EXPAND, api_sign_sk_expand_t);
	if (res) {
		return res;
	}
	SNOVA_PROFILE_BEGIN(api_sign_rng_salt_t);
	randombytes(salt, BYTES_SALT);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_SIGN_RNG_SALT, api_sign_rng_salt_t);
#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_UNDEFINED(salt, BYTES_SALT);
#endif
	SNOVA_PROFILE_BEGIN(api_sign_hash_msg_t);
	shake256(digest, BYTES_DIGEST, m, mlen);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_SIGN_HASH_MSG, api_sign_hash_msg_t);
	SNOVA_PROFILE_BEGIN(api_sign_core_t);
	res = SNOVA_NAMESPACE(sign)(&skx, sm, digest, BYTES_DIGEST, salt);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_SIGN_CORE, api_sign_core_t);
	if (!res) {
		SNOVA_PROFILE_BEGIN(api_sign_append_t);
		memcpy(sm + CRYPTO_BYTES, m, mlen);
		*smlen = mlen + CRYPTO_BYTES;
		SNOVA_PROFILE_END(SNOVA_PROFILE_API_SIGN_APPEND, api_sign_append_t);
	}
	return res;
}
int crypto_sign_open(unsigned char* m, unsigned long long* mlen, const unsigned char* sm, unsigned long long smlen,
                     const unsigned char* pk) {
	uint8_t digest[BYTES_DIGEST];
	expanded_PK pkx;
	if (smlen < CRYPTO_BYTES) {
		return -1;
	}
	SNOVA_PROFILE_BEGIN(api_open_pk_expand_t);
	int res = SNOVA_NAMESPACE(pk_expand)(&pkx, pk);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_OPEN_PK_EXPAND, api_open_pk_expand_t);
	if (res) {
		return -1;
	}
	SNOVA_PROFILE_BEGIN(api_open_hash_msg_t);
	shake256(digest, BYTES_DIGEST, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_OPEN_HASH_MSG, api_open_hash_msg_t);
	SNOVA_PROFILE_BEGIN(api_open_verify_t);
	res = SNOVA_NAMESPACE(verify)(&pkx, sm, digest, BYTES_DIGEST);
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_OPEN_VERIFY, api_open_verify_t);
	if (res) {
		return -1;
	}
	SNOVA_PROFILE_BEGIN(api_open_copy_t);
	memcpy(m, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
	*mlen = smlen - CRYPTO_BYTES;
	SNOVA_PROFILE_END(SNOVA_PROFILE_API_OPEN_COPY, api_open_copy_t);
	return 0;
}
