#ifndef SNOVA_VECOPS_H
#define SNOVA_VECOPS_H

#include <stdalign.h>
#include <stdint.h>
#include <string.h>

#include "snova.h"

#define v_SNOVA SNOVA_v
#define o_SNOVA SNOVA_o
#define m_SNOVA SNOVA_o
#define l_SNOVA SNOVA_l
#define n_SNOVA SNOVA_n
#define lsq_SNOVA SNOVA_l2
#define alpha_SNOVA SNOVA_alpha

#define rank SNOVA_l
#define sq_rank SNOVA_l2

#define seed_length_public SEED_LENGTH_PUBLIC
#define seed_length_private SEED_LENGTH_PRIVATE
#define seed_length SEED_LENGTH
#define bytes_salt BYTES_SALT

#define i_prime_inv(mi, alpha) ((o_SNOVA * alpha_SNOVA - (alpha) + (mi)) % o_SNOVA)

typedef uint8_t gf16_t;
typedef gf16_t gf16m_t[sq_rank];

#define get_gf16m(gf16m, x, y) (gf16m[(((x) * rank) + (y))])
#define set_gf16m(gf16m, x, y, value) (gf16m[(((x) * rank) + (y))] = value)
#define gf16_get_add(a, b) ((a) ^ (b))

extern uint8_t mt4b[256];
extern uint8_t inv4b[16];

#define mt(p, q) mt4b[((p) << 4) ^ (q)]
#define inv(gf16) inv4b[(gf16)]

#if defined(__AVX2__)
#include <immintrin.h>
typedef __m256i VECTOR;
#define VEC_LENGTH 32
#define VEC_BYTE(vect, idx) (((uint8_t*)vect)[idx])
#define VEC_SHUFFLE _mm256_shuffle_epi8
#define VEC_CMP_GT _mm256_cmpgt_epi32
#define VEC_CMP_EQ _mm256_cmpeq_epi8
static inline VECTOR vec_set1_u8(uint8_t value) {
	return _mm256_set1_epi8((char)value);
}
static inline VECTOR vec_test_bits(VECTOR value, VECTOR mask) {
	return VEC_CMP_EQ(value & mask, mask);
}
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
typedef uint8x16_t VECTOR;
#define VEC_LENGTH 16
#define VEC_BYTE(vect, idx) vect[(idx) / 16][(idx) % 16]
#define VEC_SHUFFLE vqtbl1q_u8
#define VEC_CMP_GT vcgtq_u8
#define VEC_CMP_EQ vceqq_u8
static inline VECTOR vec_set1_u8(uint8_t value) {
	return vdupq_n_u8(value);
}
static inline VECTOR vec_test_bits(VECTOR value, VECTOR mask) {
	return vtstq_u8(value, mask);
}
#else
#error "Vectorization not supported"
#endif

#if (defined(__ARM_NEON) || defined(__ARM_NEON__)) && (rank == 4)
#define SNOVA_NEON_RANK4 1
#else
#define SNOVA_NEON_RANK4 0
#endif

#if SNOVA_NEON_RANK4 || (defined(__AVX2__) && rank == 4)
#define SNOVA_RANK4_FAST 1
#else
#define SNOVA_RANK4_FAST 0
#endif

#if defined(__AVX2__) && rank == 4
#define SNOVA_AVX2_RANK4 1
#else
#define SNOVA_AVX2_RANK4 0
#endif

#define vtl_mainRow_x_rank(mainRow) (((mainRow) * rank + VEC_LENGTH - 1) / VEC_LENGTH * VEC_LENGTH)
#define vtl_mainRow_x_rank32(mainRow) (((mainRow) * rank + VEC_LENGTH - 1) / VEC_LENGTH)

#define vtl_len vtl_mainRow_x_rank32(SNOVA_n)
#define vtl_v_len vtl_mainRow_x_rank32(SNOVA_v)
#define vtl_o_len vtl_mainRow_x_rank32(SNOVA_o)

#define rank_next32 ((rank + 1) / 2)
#define rank_next2 (rank_next32 * 2)
#define rank_floor2 ((rank) / 2 * 2)
#define alpha_SNOVA_next2 ((alpha_SNOVA + 1) / 2 * 2)

extern VECTOR *mtk2_16;
extern VECTOR vector_inv_table;
extern VECTOR l_mask;
extern VECTOR vtl_multmask1, vtl_multmask2, vtl_multmask4, vtl_multmask8;
extern VECTOR vtl_mult_table1, vtl_mult_table2, vtl_mult_table4, vtl_mult_table8;
extern VECTOR zero256;

void snova_init_gf16_tables(void);
void snova_init_vector_table(void);
void snova_gen_S_array(void);

static inline VECTOR vtl_ct_multtab(uint8_t val) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	VECTOR val256 = vld1q_dup_u8(&val);
#else
	VECTOR val256 = _mm256_set1_epi8(val);
#endif

	return (vtl_mult_table1 & VEC_CMP_GT(val256 & vtl_multmask1, zero256)) ^
	       (vtl_mult_table2 & VEC_CMP_GT(val256 & vtl_multmask2, zero256)) ^
	       (vtl_mult_table4 & VEC_CMP_GT(val256 & vtl_multmask4, zero256)) ^
	       (vtl_mult_table8 & VEC_CMP_GT(val256 & vtl_multmask8, zero256));
}

#ifndef SNOVA_VTL_MULTTAB_LOOKUP
#define SNOVA_VTL_MULTTAB_LOOKUP 1
#endif

#if SNOVA_VTL_MULTTAB_LOOKUP
#define SNOVA_MULTTAB(k) (mtk2_16[(k)])
#else
#define SNOVA_MULTTAB(k) vtl_ct_multtab((uint8_t)(k))
#endif

static inline void gf16_32_mul_k(VECTOR* a_256, gf16_t k, VECTOR* ak) {
	*ak = VEC_SHUFFLE(SNOVA_MULTTAB(k), *a_256);
}

static inline void gf16_32_mul_k_add(VECTOR* a_256, gf16_t k, VECTOR* ak) {
	*ak ^= VEC_SHUFFLE(SNOVA_MULTTAB(k), *a_256);
}

static inline void gf16_32_mul_tab(VECTOR* a_256, VECTOR tab, VECTOR* ak) {
	*ak = VEC_SHUFFLE(tab, *a_256);
}

static inline void gf16_32_mul_tab_add(VECTOR* a_256, VECTOR tab, VECTOR* ak) {
	*ak ^= VEC_SHUFFLE(tab, *a_256);
}

static inline void gf16_32_mul_32(VECTOR* a_256, VECTOR* b_256, VECTOR* sum) {
	memset(sum, 0, sizeof(VECTOR));

	*sum ^= *b_256 & VEC_CMP_EQ(*a_256 & vtl_multmask1, vtl_multmask1);
	*sum ^= VEC_SHUFFLE(vtl_mult_table2, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask2, vtl_multmask2);
	*sum ^= VEC_SHUFFLE(vtl_mult_table4, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask4, vtl_multmask4);
	*sum ^= VEC_SHUFFLE(vtl_mult_table8, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask8, vtl_multmask8);
}

static inline void gf16_32_mul_32_add(VECTOR* a_256, VECTOR* b_256, VECTOR* sum) {
	*sum ^= *b_256 & VEC_CMP_EQ(*a_256 & vtl_multmask1, vtl_multmask1);
	*sum ^= VEC_SHUFFLE(vtl_mult_table2, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask2, vtl_multmask2);
	*sum ^= VEC_SHUFFLE(vtl_mult_table4, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask4, vtl_multmask4);
	*sum ^= VEC_SHUFFLE(vtl_mult_table8, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask8, vtl_multmask8);
}

static inline VECTOR gf16_32_mul_32_precomputed_b(VECTOR a, VECTOR b1, VECTOR b2, VECTOR b4, VECTOR b8) {
	return (b1 & VEC_CMP_EQ(a & vtl_multmask1, vtl_multmask1)) ^
	       (b2 & VEC_CMP_EQ(a & vtl_multmask2, vtl_multmask2)) ^
	       (b4 & VEC_CMP_EQ(a & vtl_multmask4, vtl_multmask4)) ^
	       (b8 & VEC_CMP_EQ(a & vtl_multmask8, vtl_multmask8));
}

#if rank == 4
#ifndef SNOVA_PACKED_MUL2_ARITH
#define SNOVA_PACKED_MUL2_ARITH 0
#endif

static inline VECTOR gf16_32_packed_mul2_arith(VECTOR packed) {
	VECTOR keep = vec_set1_u8(0x77);
	VECTOR carry = packed & vec_set1_u8(0x88);
	return ((packed & keep) << 1) ^ (carry >> 2) ^ (carry >> 3);
}

static inline VECTOR gf16_32_packed_mul_k(VECTOR packed, gf16_t k) {
#if SNOVA_PACKED_MUL2_ARITH
	if (k == 2) {
		return gf16_32_packed_mul2_arith(packed);
	}
	if (k == 4) {
		return gf16_32_packed_mul2_arith(gf16_32_packed_mul2_arith(packed));
	}
	if (k == 8) {
		return gf16_32_packed_mul2_arith(gf16_32_packed_mul2_arith(gf16_32_packed_mul2_arith(packed)));
	}
#endif
	VECTOR lo = VEC_SHUFFLE(SNOVA_MULTTAB(k), packed & l_mask);
	VECTOR hi = VEC_SHUFFLE(SNOVA_MULTTAB(k), (packed >> 4) & l_mask);
	return lo ^ (hi << 4);
}

static inline VECTOR gf16_32_mul_32_precomputed_packed_b(VECTOR a, VECTOR b1, VECTOR b2, VECTOR b4, VECTOR b8) {
	return (b1 & VEC_CMP_EQ(a & vtl_multmask1, vtl_multmask1)) ^
	       (b2 & VEC_CMP_EQ(a & vtl_multmask2, vtl_multmask2)) ^
	       (b4 & VEC_CMP_EQ(a & vtl_multmask4, vtl_multmask4)) ^
	       (b8 & VEC_CMP_EQ(a & vtl_multmask8, vtl_multmask8));
}

#ifndef SNOVA_NEON_PACKED_A_VTST
#define SNOVA_NEON_PACKED_A_VTST 1
#endif

static inline void gf16_32_mul_32_split_packed_a(VECTOR a_packed, VECTOR b1, VECTOR b2, VECTOR b4, VECTOR b8,
                                                 VECTOR* out_lo, VECTOR* out_hi) {
#if (defined(__ARM_NEON) || defined(__ARM_NEON__)) && SNOVA_NEON_PACKED_A_VTST
	const VECTOR m1_lo = vec_set1_u8(0x01);
	const VECTOR m2_lo = vec_set1_u8(0x02);
	const VECTOR m4_lo = vec_set1_u8(0x04);
	const VECTOR m8_lo = vec_set1_u8(0x08);
	const VECTOR m1_hi = vec_set1_u8(0x10);
	const VECTOR m2_hi = vec_set1_u8(0x20);
	const VECTOR m4_hi = vec_set1_u8(0x40);
	const VECTOR m8_hi = vec_set1_u8(0x80);

	*out_lo = ((b1 & vec_test_bits(a_packed, m1_lo)) ^
	           (b2 & vec_test_bits(a_packed, m2_lo))) ^
	          ((b4 & vec_test_bits(a_packed, m4_lo)) ^
	           (b8 & vec_test_bits(a_packed, m8_lo)));
	*out_hi = ((b1 & vec_test_bits(a_packed, m1_hi)) ^
	           (b2 & vec_test_bits(a_packed, m2_hi))) ^
	          ((b4 & vec_test_bits(a_packed, m4_hi)) ^
	           (b8 & vec_test_bits(a_packed, m8_hi)));
#else
	VECTOR a_lo = a_packed & l_mask;
	VECTOR a_hi = (a_packed >> 4) & l_mask;

	*out_lo = (b1 & VEC_CMP_EQ(a_lo & vtl_multmask1, vtl_multmask1)) ^
	          (b2 & VEC_CMP_EQ(a_lo & vtl_multmask2, vtl_multmask2)) ^
	          (b4 & VEC_CMP_EQ(a_lo & vtl_multmask4, vtl_multmask4)) ^
	          (b8 & VEC_CMP_EQ(a_lo & vtl_multmask8, vtl_multmask8));
	*out_hi = (b1 & VEC_CMP_EQ(a_hi & vtl_multmask1, vtl_multmask1)) ^
	          (b2 & VEC_CMP_EQ(a_hi & vtl_multmask2, vtl_multmask2)) ^
	          (b4 & VEC_CMP_EQ(a_hi & vtl_multmask4, vtl_multmask4)) ^
	          (b8 & VEC_CMP_EQ(a_hi & vtl_multmask8, vtl_multmask8));
#endif
}

static inline void gf16_32_prepare_packed_a_masks(
    VECTOR a_packed,
    VECTOR* m1_lo, VECTOR* m2_lo, VECTOR* m4_lo, VECTOR* m8_lo,
    VECTOR* m1_hi, VECTOR* m2_hi, VECTOR* m4_hi, VECTOR* m8_hi) {
#if (defined(__ARM_NEON) || defined(__ARM_NEON__)) && SNOVA_NEON_PACKED_A_VTST
	*m1_lo = vec_test_bits(a_packed, vec_set1_u8(0x01));
	*m2_lo = vec_test_bits(a_packed, vec_set1_u8(0x02));
	*m4_lo = vec_test_bits(a_packed, vec_set1_u8(0x04));
	*m8_lo = vec_test_bits(a_packed, vec_set1_u8(0x08));
	*m1_hi = vec_test_bits(a_packed, vec_set1_u8(0x10));
	*m2_hi = vec_test_bits(a_packed, vec_set1_u8(0x20));
	*m4_hi = vec_test_bits(a_packed, vec_set1_u8(0x40));
	*m8_hi = vec_test_bits(a_packed, vec_set1_u8(0x80));
#else
	const VECTOR a_lo = a_packed & l_mask;
	const VECTOR a_hi = (a_packed >> 4) & l_mask;
	*m1_lo = VEC_CMP_EQ(a_lo & vtl_multmask1, vtl_multmask1);
	*m2_lo = VEC_CMP_EQ(a_lo & vtl_multmask2, vtl_multmask2);
	*m4_lo = VEC_CMP_EQ(a_lo & vtl_multmask4, vtl_multmask4);
	*m8_lo = VEC_CMP_EQ(a_lo & vtl_multmask8, vtl_multmask8);
	*m1_hi = VEC_CMP_EQ(a_hi & vtl_multmask1, vtl_multmask1);
	*m2_hi = VEC_CMP_EQ(a_hi & vtl_multmask2, vtl_multmask2);
	*m4_hi = VEC_CMP_EQ(a_hi & vtl_multmask4, vtl_multmask4);
	*m8_hi = VEC_CMP_EQ(a_hi & vtl_multmask8, vtl_multmask8);
#endif
}

static inline VECTOR gf16_32_apply_packed_a_masks(
    VECTOR b1, VECTOR b2, VECTOR b4, VECTOR b8,
    VECTOR m1, VECTOR m2, VECTOR m4, VECTOR m8) {
	const VECTOR x12 = (b1 & m1) ^ (b2 & m2);
	const VECTOR x48 = (b4 & m4) ^ (b8 & m8);
	return x12 ^ x48;
}

static inline void gf16_32_mul_32_add_4_same_packed_a_values_masked(
    const VECTOR* restrict b_packed,
    VECTOR m01_1, VECTOR m01_2, VECTOR m01_4, VECTOR m01_8,
    VECTOR m01_1_hi, VECTOR m01_2_hi, VECTOR m01_4_hi, VECTOR m01_8_hi,
    VECTOR m23_1, VECTOR m23_2, VECTOR m23_4, VECTOR m23_8,
    VECTOR m23_1_hi, VECTOR m23_2_hi, VECTOR m23_4_hi, VECTOR m23_8_hi,
    VECTOR* restrict sum0, VECTOR* restrict sum1,
    VECTOR* restrict sum2, VECTOR* restrict sum3) {
	const VECTOR b1 = *b_packed;
	const VECTOR b2 = gf16_32_packed_mul_k(b1, 2);
	const VECTOR b4 = gf16_32_packed_mul_k(b1, 4);
	const VECTOR b8 = gf16_32_packed_mul_k(b1, 8);
	const VECTOR lo0 = gf16_32_apply_packed_a_masks(b1, b2, b4, b8,
	                                                  m01_1, m01_2, m01_4, m01_8);
	const VECTOR hi0 = gf16_32_apply_packed_a_masks(b1, b2, b4, b8,
	                                                  m01_1_hi, m01_2_hi, m01_4_hi, m01_8_hi);
	const VECTOR lo1 = gf16_32_apply_packed_a_masks(b1, b2, b4, b8,
	                                                  m23_1, m23_2, m23_4, m23_8);
	const VECTOR hi1 = gf16_32_apply_packed_a_masks(b1, b2, b4, b8,
	                                                  m23_1_hi, m23_2_hi, m23_4_hi, m23_8_hi);
	*sum0 ^= lo0;
	*sum1 ^= hi0;
	*sum2 ^= lo1;
	*sum3 ^= hi1;
}

static inline VECTOR gf16_32_mul_32_precomputed_packed_a_cmp(VECTOR a_packed, VECTOR b1, VECTOR b2, VECTOR b4,
                                                             VECTOR b8) {
	VECTOR out_lo;
	VECTOR out_hi;
	gf16_32_mul_32_split_packed_a(a_packed, b1, b2, b4, b8, &out_lo, &out_hi);
	return out_lo ^ (out_hi << 4);
}

static inline void gf16_32_mul_32_add_4_same_packed_a_values(VECTOR a01, VECTOR a23,
                                                             const VECTOR* restrict b_packed,
                                                             VECTOR* restrict sum0, VECTOR* restrict sum1,
                                                             VECTOR* restrict sum2, VECTOR* restrict sum3) {
	VECTOR b1 = *b_packed;
	VECTOR b2 = gf16_32_packed_mul_k(b1, 2);
	VECTOR b4 = gf16_32_packed_mul_k(b1, 4);
	VECTOR b8 = gf16_32_packed_mul_k(b1, 8);

	VECTOR lo0;
	VECTOR hi0;
	VECTOR lo1;
	VECTOR hi1;

	gf16_32_mul_32_split_packed_a(a01, b1, b2, b4, b8, &lo0, &hi0);
	gf16_32_mul_32_split_packed_a(a23, b1, b2, b4, b8, &lo1, &hi1);

	*sum0 ^= lo0;
	*sum1 ^= hi0;
	*sum2 ^= lo1;
	*sum3 ^= hi1;
}

static inline void gf16_32_mul_32_add_4_same_packed_b_values(VECTOR a0, VECTOR a1, VECTOR a2, VECTOR a3,
                                                             const VECTOR* restrict b_packed,
                                                             VECTOR* restrict sum0, VECTOR* restrict sum1,
                                                             VECTOR* restrict sum2, VECTOR* restrict sum3) {
	VECTOR b1 = *b_packed;
	VECTOR b2 = gf16_32_packed_mul_k(b1, 2);
	VECTOR b4 = gf16_32_packed_mul_k(b1, 4);
	VECTOR b8 = gf16_32_packed_mul_k(b1, 8);

	*sum0 ^= gf16_32_mul_32_precomputed_packed_b(a0, b1, b2, b4, b8);
	*sum1 ^= gf16_32_mul_32_precomputed_packed_b(a1, b1, b2, b4, b8);
	*sum2 ^= gf16_32_mul_32_precomputed_packed_b(a2, b1, b2, b4, b8);
	*sum3 ^= gf16_32_mul_32_precomputed_packed_b(a3, b1, b2, b4, b8);
}

static inline void gf16_32_mul_32_add_2_same_b(const VECTOR* restrict a0_256, const VECTOR* restrict a1_256,
                                               const VECTOR* restrict b_256, VECTOR* restrict sum0,
                                               VECTOR* restrict sum1) {
	VECTOR b1 = *b_256;
	VECTOR b2 = VEC_SHUFFLE(vtl_mult_table2, b1);
	VECTOR b4 = VEC_SHUFFLE(vtl_mult_table4, b1);
	VECTOR b8 = VEC_SHUFFLE(vtl_mult_table8, b1);

	*sum0 ^= gf16_32_mul_32_precomputed_b(*a0_256, b1, b2, b4, b8);
	*sum1 ^= gf16_32_mul_32_precomputed_b(*a1_256, b1, b2, b4, b8);
}

static inline void gf16_32_mul_32_add_4_same_b_values(VECTOR a0, VECTOR a1, VECTOR a2, VECTOR a3,
                                                      const VECTOR* restrict b_256, VECTOR* restrict sum0,
                                                      VECTOR* restrict sum1, VECTOR* restrict sum2,
                                                      VECTOR* restrict sum3) {
	VECTOR b1 = *b_256;
	VECTOR b2 = VEC_SHUFFLE(vtl_mult_table2, b1);
	VECTOR b4 = VEC_SHUFFLE(vtl_mult_table4, b1);
	VECTOR b8 = VEC_SHUFFLE(vtl_mult_table8, b1);

	*sum0 ^= gf16_32_mul_32_precomputed_b(a0, b1, b2, b4, b8);
	*sum1 ^= gf16_32_mul_32_precomputed_b(a1, b1, b2, b4, b8);
	*sum2 ^= gf16_32_mul_32_precomputed_b(a2, b1, b2, b4, b8);
	*sum3 ^= gf16_32_mul_32_precomputed_b(a3, b1, b2, b4, b8);
}
#endif

static inline uint8_t horizontal_xor_256(VECTOR vec) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	int64x2_t combined = (int64x2_t)vec;
#elif VEC_LENGTH > 16
	__m128i low = _mm256_castsi256_si128(vec);
	__m128i high = _mm256_extracti128_si256(vec, 1);
	__m128i combined = _mm_xor_si128(low, high);
#else
	__m128i combined = vec;
#endif

	combined ^= combined >> 32;
	combined ^= combined >> 16;
	combined ^= combined >> 8;

	return combined[0] ^ combined[1];
}

#ifndef SNOVA_GF16_SCALAR_MODEXP
#define SNOVA_GF16_SCALAR_MODEXP 1
#endif

#ifndef SNOVA_GF16_MODEXP_NEON
#define SNOVA_GF16_MODEXP_NEON 1
#endif

#ifndef SNOVA_GF16_MODEXP_AVX2
#define SNOVA_GF16_MODEXP_AVX2 1
#endif

#if SNOVA_GF16_SCALAR_MODEXP

alignas(32) static const uint8_t gf16_modexp_log16[16] = {
	0, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12
};

alignas(32) static const uint8_t gf16_modexp_exp16[15] = {
	1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9
};

static inline gf16_t gf16_mul_modexp_scalar(gf16_t a, gf16_t b) {
	a &= 0x0f;
	b &= 0x0f;

	uint8_t e = (uint8_t)(gf16_modexp_log16[a] + gf16_modexp_log16[b]);
	uint8_t ge15_mask = (uint8_t)-(uint8_t)(e >= 15);
	e = (uint8_t)(e - (uint8_t)(15u & ge15_mask));

	uint8_t nonzero_mask = (uint8_t)-(uint8_t)((a != 0) & (b != 0));
	return (gf16_t)(gf16_modexp_exp16[e] & nonzero_mask);
}

#if defined(__ARM_NEON) && SNOVA_GF16_MODEXP_NEON
static inline uint8x16_t gf16_mul_vec_modexp_scalar_neon(uint8x16_t values, uint8_t factor) {
	static const uint8_t log_lut_bytes[16] = {
		0, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12
	};
	static const uint8_t exp_lut_bytes[16] = {
		1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9, 0
	};
	const uint8x16_t log_lut = vld1q_u8(log_lut_bytes);
	const uint8x16_t exp_lut = vld1q_u8(exp_lut_bytes);
	const uint8x16_t zero = vdupq_n_u8(0);
	const uint8x16_t low_mask = vdupq_n_u8(0x0f);

	values = vandq_u8(values, low_mask);
	factor &= 0x0f;
	if (factor == 0) {
		return zero;
	}

	uint8x16_t logs = vqtbl1q_u8(log_lut, values);
	uint8x16_t sum = vaddq_u8(logs, vdupq_n_u8(log_lut_bytes[factor]));
	uint8x16_t ge15 = vcgeq_u8(sum, vdupq_n_u8(15));
	sum = vsubq_u8(sum, vandq_u8(ge15, vdupq_n_u8(15)));

	uint8x16_t res = vqtbl1q_u8(exp_lut, sum);
	uint8x16_t nonzero = vmvnq_u8(vceqq_u8(values, zero));
	return vandq_u8(res, nonzero);
}

static inline uint8x16_t gf16_vec_from_u32(uint32_t word) {
	return vreinterpretq_u8_u32(vsetq_lane_u32(word, vdupq_n_u32(0), 0));
}

static inline uint32_t gf16_vec_low_u32(uint8x16_t vec) {
	return vget_lane_u32(vreinterpret_u32_u8(vget_low_u8(vec)), 0);
}
#endif

#if defined(__AVX2__) && SNOVA_GF16_MODEXP_AVX2 && rank == 4
static inline __m128i gf16_mul_vec_modexp_scalar_avx2(__m128i values, uint8_t factor) {
	static const uint8_t log_lut_bytes[16] = {
		0, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12
	};
	static const uint8_t exp_lut_bytes[16] = {
		1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9, 0
	};
	const __m128i log_lut = _mm_loadu_si128((const __m128i*)log_lut_bytes);
	const __m128i exp_lut = _mm_loadu_si128((const __m128i*)exp_lut_bytes);
	const __m128i zero = _mm_setzero_si128();
	const __m128i low_mask = _mm_set1_epi8(0x0f);

	values = _mm_and_si128(values, low_mask);
	factor &= 0x0f;
	if (factor == 0) {
		return zero;
	}

	__m128i logs = _mm_shuffle_epi8(log_lut, values);
	__m128i sum = _mm_add_epi8(logs, _mm_set1_epi8((char)log_lut_bytes[factor]));
	__m128i ge15 = _mm_cmpgt_epi8(sum, _mm_set1_epi8(14));
	sum = _mm_sub_epi8(sum, _mm_and_si128(ge15, _mm_set1_epi8(15)));
	__m128i result = _mm_shuffle_epi8(exp_lut, sum);
	__m128i nonzero = _mm_andnot_si128(_mm_cmpeq_epi8(values, zero), _mm_set1_epi8((char)0xff));
	return _mm_and_si128(result, nonzero);
}
#endif

#define gf16_get_mul(a, b) gf16_mul_modexp_scalar((a), (b))
#else
#define gf16_get_mul(a, b) (mt((a), (b)))
#endif

static inline uint32_t gf16_reduce(uint32_t idx) {
	uint32_t res, upper;

	res = idx & 0x49249249;
	upper = idx >> 12;
	res = res ^ upper ^ (upper << 3);
	upper = res >> 12;
	res = res ^ upper ^ (upper << 3);
	upper = res >> 12;
	res = res ^ upper ^ (upper << 3);

	return res & 0x249;
}

static inline uint32_t gf16_from_nibble(uint8_t idx) {
	uint32_t middle = idx | idx << 4;
	return (middle & 0x41) | ((middle << 2) & 0x208);
}

static inline uint8_t gf16_to_nibble(uint32_t idx) {
	uint32_t res = gf16_reduce(idx);
	res = res | (res >> 4);
	return (res & 0x5) | ((res >> 2) & 0xa);
}

static inline uint8_t xgf16_to_nibble(uint32_t res) {
	res = res | (res >> 4);
	return (res & 0x5) | ((res >> 2) & 0xa);
}

static inline uint32_t gf16_inv(uint32_t val) {
	val = gf16_reduce(val);
	uint32_t res2 = gf16_reduce(val * val);
	uint32_t res4 = gf16_reduce(res2 * res2);
	uint32_t res8 = gf16_reduce(res4 * res4);

	return gf16_reduce(res2 * ((res4 * res8) & 0x49249249));
}

static inline uint32_t ct_gf16_is_not_zero(uint8_t val) {
	return val != 0;
}

static inline void gf16m_set_zero(gf16m_t a) {
	memset(a, 0, sq_rank);
}

static inline void gf16m_clone(gf16m_t target, const gf16m_t source) {
	memcpy(target, source, sq_rank);
}

static inline void gf16m_add(const gf16m_t a, const gf16m_t b, gf16m_t c) {
#if defined(__ARM_NEON) && SNOVA_GF16_MODEXP_NEON && rank == 4
	uint8x16_t va = vld1q_u8((const uint8_t*)a);
	uint8x16_t vb = vld1q_u8((const uint8_t*)b);
	vst1q_u8((uint8_t*)c, veorq_u8(va, vb));
#elif defined(__AVX2__) && SNOVA_GF16_MODEXP_AVX2 && rank == 4
	__m128i va = _mm_loadu_si128((const __m128i*)a);
	__m128i vb = _mm_loadu_si128((const __m128i*)b);
	_mm_storeu_si128((__m128i*)c, _mm_xor_si128(va, vb));
#else
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(c, i, j, gf16_get_add(get_gf16m(a, i, j), get_gf16m(b, i, j)));
		}
	}
#endif
}

static inline void gf16m_mul(const gf16m_t a, const gf16m_t b, gf16m_t c) {
#if defined(__ARM_NEON) && SNOVA_GF16_MODEXP_NEON && rank == 4
	uint8x16_t b_row_vec[4];
	for (int k = 0; k < 4; ++k) {
		uint32_t row_word = 0;
		memcpy(&row_word, (const uint8_t*)b + k * 4, sizeof(row_word));
		b_row_vec[k] = gf16_vec_from_u32(row_word);
	}

	for (int i = 0; i < 4; ++i) {
		uint32_t acc_word = 0;
		for (int k = 0; k < 4; ++k) {
			uint8x16_t prod = gf16_mul_vec_modexp_scalar_neon(b_row_vec[k], get_gf16m(a, i, k));
			acc_word ^= gf16_vec_low_u32(prod);
		}
		memcpy((uint8_t*)c + i * 4, &acc_word, sizeof(acc_word));
	}
#elif defined(__AVX2__) && SNOVA_GF16_MODEXP_AVX2 && rank == 4
	__m128i b_row_vec[4];
	for (int k = 0; k < 4; ++k) {
		uint32_t row_word;
		memcpy(&row_word, (const uint8_t*)b + k * 4, sizeof(row_word));
		b_row_vec[k] = _mm_cvtsi32_si128((int)row_word);
	}
	for (int i = 0; i < 4; ++i) {
		__m128i acc = _mm_setzero_si128();
		for (int k = 0; k < 4; ++k) {
			acc = _mm_xor_si128(acc,
			                    gf16_mul_vec_modexp_scalar_avx2(b_row_vec[k], get_gf16m(a, i, k)));
		}
		const uint32_t row_word = (uint32_t)_mm_cvtsi128_si32(acc);
		memcpy((uint8_t*)c + i * 4, &row_word, sizeof(row_word));
	}
#else
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(c, i, j, gf16_get_mul(get_gf16m(a, i, 0), get_gf16m(b, 0, j)));
			for (int k = 1; k < rank; ++k) {
				set_gf16m(c, i, j, gf16_get_add(get_gf16m(c, i, j), gf16_get_mul(get_gf16m(a, i, k), get_gf16m(b, k, j))));
			}
		}
	}
#endif
}

static inline void gf16m_scale(const gf16m_t a, gf16_t k, gf16m_t c) {
#if defined(__ARM_NEON) && SNOVA_GF16_MODEXP_NEON && rank == 4
	uint8x16_t vec = vld1q_u8((const uint8_t*)a);
	vec = gf16_mul_vec_modexp_scalar_neon(vec, k);
	vst1q_u8((uint8_t*)c, vec);
#elif defined(__AVX2__) && SNOVA_GF16_MODEXP_AVX2 && rank == 4
	__m128i vec = _mm_loadu_si128((const __m128i*)a);
	vec = gf16_mul_vec_modexp_scalar_avx2(vec, k);
	_mm_storeu_si128((__m128i*)c, vec);
#else
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(c, i, j, gf16_get_mul(get_gf16m(a, i, j), k));
		}
	}
#endif
}

static inline void gf16m_transpose(const gf16m_t a, gf16m_t ap) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(ap, i, j, get_gf16m(a, j, i));
		}
	}
}

static inline void be_aI(gf16m_t target, gf16_t a) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(target, i, j, (i == j) ? a : 0);
		}
	}
}

static inline void be_the_S(gf16m_t target) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(target, i, j, (8 - (i + j)));
		}
	}
#if rank == 5
	set_gf16m(target, 4, 4, 9);
#endif
}

static inline void be_aS(gf16m_t target, gf16_t a) {
#if defined(__ARM_NEON) && SNOVA_GF16_MODEXP_NEON && rank == 4
	static const uint8_t coeff_bytes[16] = {
		8, 7, 6, 5,
		7, 6, 5, 4,
		6, 5, 4, 3,
		5, 4, 3, 2
	};
	uint8x16_t coeff = vld1q_u8(coeff_bytes);
	coeff = gf16_mul_vec_modexp_scalar_neon(coeff, a);
	vst1q_u8((uint8_t*)target, coeff);
#elif defined(__AVX2__) && SNOVA_GF16_MODEXP_AVX2 && rank == 4
	static const uint8_t coeff_bytes[16] = {
		8, 7, 6, 5,
		7, 6, 5, 4,
		6, 5, 4, 3,
		5, 4, 3, 2
	};
	__m128i coeff = _mm_loadu_si128((const __m128i*)coeff_bytes);
	coeff = gf16_mul_vec_modexp_scalar_avx2(coeff, a);
	_mm_storeu_si128((__m128i*)target, coeff);
#else
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(target, i, j, gf16_get_mul((8 - (i + j)), a));
		}
	}
#if rank == 5
	set_gf16m(target, 4, 4, gf16_get_mul(9, a));
#endif
#endif
}

extern gf16m_t snova_S[l_SNOVA];
extern uint32_t snova_xS[l_SNOVA][lsq_SNOVA];

gf16_t snova_gf16m_det(const gf16m_t entry);
void snova_be_invertible_by_add_aS(gf16m_t source);
void snova_gen_a_FqS(const gf16_t* c, gf16m_t pt_matrix);
void snova_gen_a_FqS_ct(const gf16_t* c, gf16m_t pt_matrix);

#endif
