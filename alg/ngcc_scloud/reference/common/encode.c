/**
 * @file encode.c
 * @brief Message encoder/decoder selected by compile-time BW code parameters.
 *
 * The level parameter header selects BW32, BW128 or RBW128 through compact
 * numeric macros.  This file keeps one implementation shape and lets the
 * compiler fold the few public-parameter choices at build time.
 */

#include "encode.h"
#include "gf/mod_pow2.h"
#include "scloudplus_param_common.h"
#include "util/bits.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/*
 * Helpers below are independent of the active BW/RBW message layout. Keeping
 * one copy avoids the old per-parameter source-file shape while preserving the
 * compile-time-selected branches for layout-dependent code.
 */
/** @brief Add two integer complex coordinates. */
static inline Complex complex_add(Complex a, Complex b)
{
	return (Complex){a.real + b.real, a.imag + b.imag};
}

/** @brief Subtract two integer complex coordinates. */
static inline Complex complex_sub(Complex a, Complex b)
{
	return (Complex){a.real - b.real, a.imag - b.imag};
}

/**
 * @brief Multiply by phi = 1+i in integer complex coordinates.
 *
 * (a+bi)(1+i) = (a-b) + (a+b)i.  No division is involved, so this operation is
 * exact in the integer representation.
 */
static inline Complex complex_mul_phi(Complex z)
{
	return (Complex){z.real - z.imag, z.real + z.imag};
}

/* BDD rounding relies on arithmetic right shift for signed powers of two. */
#if ((-1 >> 1) != -1)
#error "Scloud+ codec requires arithmetic right shift for signed integers"
#endif

/**
 * @brief Multiply by phi^{-1} = (1-i)/2 on integral lattice coordinates.
 *
 * This helper is used after a BW lattice point has already been recovered.
 * At that point the numerator is known to be divisible by two in both
 * coordinates, so the integer division is exact.
 */
static inline Complex complex_mul_phi_inv(Complex z)
{
	return (Complex){(int32_t)(((int64_t)z.real + z.imag) >> 1U),
					 (int32_t)(((int64_t)z.imag - z.real) >> 1U)};
}

/**
 * @brief Multiply by 1-i without dividing by two.
 *
 * Algorithm 20 represents a target as integer numerators with an implicit
 * denominator 2^s.  Applying phi^{-1}=(1-i)/2 is therefore implemented as
 * multiplication by 1-i and an increment of the public scale s.  This keeps
 * all arithmetic integral and avoids the old fixed-point over-scaling.
 */
static inline Complex complex_mul_one_minus_i(Complex z)
{
	return (Complex){z.real + z.imag, z.imag - z.real};
}

/**
 * @brief Round a scaled integer numerator to the nearest integer.
 *
 * The represented real value is value / 2^scale_bits.  The rounding rule is
 * floor(x+1/2), matching the notation used by the specification.
 */
static inline int32_t scaled_round_to_int(int32_t value, unsigned int scale_bits)
{
	return (int32_t)(((int64_t)value +
					  ((int64_t)1 << (scale_bits - 1U))) >> scale_bits);
}

/**
 * @brief Return 2^scale_bits * round(value / 2^scale_bits).
 */
static inline int32_t round_to_scaled_integer(int32_t value,
											  unsigned int scale_bits)
{
	const int64_t rounded = scaled_round_to_int(value, scale_bits);

	return (int32_t)(rounded * ((int64_t)1 << scale_bits));
}

static inline uint32_t mask_bits(unsigned int bits)
{
	return (1U << bits) - 1U;
}

/** @brief Reduce a signed integer modulo 2^bits, returned as an unsigned word. */
static inline uint32_t mod_pow2(int32_t value, unsigned int bits)
{
	return ((uint32_t)value) & mask_bits(bits);
}

/**
 * @brief Return all-ones if lhs < rhs and zero otherwise.
 *
 * This mask is used for branch-free candidate selection in BDD decoding.  It
 * keeps the selected output independent of branch predictor behavior.
 */
static inline uint32_t ct_mask_u64_lt(uint64_t lhs, uint64_t rhs)
{
	const uint64_t diff = lhs - rhs;
	const uint64_t borrow = ((~lhs & rhs) | (~(lhs ^ rhs) & diff)) >> 63U;

	return 0U - (uint32_t)borrow;
}

static inline int32_t ct_select_i32(int32_t lhs, int32_t rhs, uint32_t mask)
{
	return (int32_t)((((uint32_t)lhs) & mask) | (((uint32_t)rhs) & ~mask));
}

static inline Complex ct_select_complex(Complex lhs, Complex rhs, uint32_t mask)
{
	return (Complex){ct_select_i32(lhs.real, rhs.real, mask),
					 ct_select_i32(lhs.imag, rhs.imag, mask)};
}

static inline unsigned int popcount_u32(uint32_t value)
{
#if defined(__GNUC__) || defined(__clang__)
	return (unsigned int)__builtin_popcount(value);
#else
	value = value - ((value >> 1U) & 0x55555555U);
	value = (value & 0x33333333U) + ((value >> 2U) & 0x33333333U);
	value = (value + (value >> 4U)) & 0x0f0f0f0fU;
	return (unsigned int)((value * 0x01010101U) >> 24U);
#endif
}

static inline unsigned int real_bits(unsigned int bits)
{
	return (bits + 1U) >> 1U;
}

#define SCLOUDPLUS_MESSAGE_CODE_BW32 1
#define SCLOUDPLUS_MESSAGE_CODE_BW128 2
#define SCLOUDPLUS_MESSAGE_CODE_RBW128 3

/*
 * The pair (bw_n, tau) identifies the message code.  tau alone is ambiguous:
 * BW32 and BW128 both use tau=3.
 */
#if (scloudplus_bw_n == 32U) && (scloudplus_tau == 3)
#define SCLOUDPLUS_MESSAGE_CODE SCLOUDPLUS_MESSAGE_CODE_BW32
#elif (scloudplus_bw_n == 128U) && (scloudplus_tau == 3)
#define SCLOUDPLUS_MESSAGE_CODE SCLOUDPLUS_MESSAGE_CODE_BW128
#elif (scloudplus_bw_n == 128U) && (scloudplus_tau == 4)
#define SCLOUDPLUS_MESSAGE_CODE SCLOUDPLUS_MESSAGE_CODE_RBW128
#else
#error "Unsupported Scloud+ message-code parameter pair"
#endif

#if (mat_subm * scloudplus_mu) != scloudplus_l
#error "Scloud+ message blocks must exactly cover the message bit length"
#endif

#if (mat_subm * scloudplus_bw_n) > (scloudplus_mbar * scloudplus_nbar)
#error "Scloud+ message code does not fit in the message matrix"
#endif

#if (scloudplus_msg_scale_bits == 0) || \
	(scloudplus_msg_scale_bits >= scloudplus_logq)
#error "Scloud+ message scale must produce a positive BDD denominator"
#endif

/** @brief Squared Euclidean distance in integer complex coordinates. */
static uint64_t squared_distance(const Complex *lhs, const Complex *rhs, int size)
{
	uint64_t sum = 0;

	for (int i = 0; i < size; i++)
	{
		const int64_t dr = (int64_t)lhs[i].real - (int64_t)rhs[i].real;
		const int64_t di = (int64_t)lhs[i].imag - (int64_t)rhs[i].imag;

		sum += (uint64_t)(dr * dr + di * di);
	}

	return sum;
}

/**
 * @brief Recursive integer bounded-distance decoder for BW_n.
 *
 * The input vector t stores integer numerators with public implicit
 * denominator 2^scale_bits.  The output vector y is a BW lattice point stored
 * with the same denominator, i.e. y belongs to 2^scale_bits * BW_n.  This is
 * the integer BDD formulation from Algorithm 20.
 *
 * The phi^{-1} recursive calls are implemented as multiplication by 1-i and
 * scale_bits+1, so no fractional complex coordinate is ever materialized.
 *
 * The selection itself is constant-time with respect to the candidate
 * distances.  The recursion depth is fixed by the public parameter bw_n.
 */
static void bddbwn(const Complex *t, Complex *y, int n,
				   unsigned int scale_bits)
{
	const int half = n >> 1;
	const int quarter = half >> 1;

	if (n == 2)
	{
		y[0] = (Complex){round_to_scaled_integer(t[0].real, scale_bits),
						 round_to_scaled_integer(t[0].imag, scale_bits)};
		return;
	}

	Complex t1[quarter];
	Complex t2[quarter];
	Complex y1[quarter];
	Complex y2[quarter];
	Complex z1[quarter];
	Complex z2[quarter];
	Complex z1in[quarter];
	Complex z2in[quarter];
	Complex out1[half];
	Complex out2[half];

	for (int i = 0; i < quarter; i++)
	{
		t1[i] = t[i];
		t2[i] = t[i + quarter];
	}

	bddbwn(t1, y1, half, scale_bits);
	bddbwn(t2, y2, half, scale_bits);

	for (int i = 0; i < quarter; i++)
	{
		z1in[i] = complex_mul_one_minus_i(complex_sub(t2[i], y1[i]));
		z2in[i] = complex_mul_one_minus_i(complex_sub(t1[i], y2[i]));
	}

	bddbwn(z1in, z1, half, scale_bits + 1U);
	bddbwn(z2in, z2, half, scale_bits + 1U);

	for (int i = 0; i < quarter; i++)
	{
		const Complex z1_aligned =
			{(int32_t)((int64_t)z1[i].real >> 1U),
			 (int32_t)((int64_t)z1[i].imag >> 1U)};
		const Complex z2_aligned =
			{(int32_t)((int64_t)z2[i].real >> 1U),
			 (int32_t)((int64_t)z2[i].imag >> 1U)};
		const Complex phi_z1 = complex_mul_phi(z1_aligned);
		const Complex phi_z2 = complex_mul_phi(z2_aligned);

		out1[i] = y1[i];
		out1[quarter + i] = complex_add(y1[i], phi_z1);
		out2[i] = complex_add(y2[i], phi_z2);
		out2[quarter + i] = y2[i];
	}

	{
		const uint64_t dist1 = squared_distance(out1, t, half);
		const uint64_t dist2 = squared_distance(out2, t, half);
		const uint32_t use_out1 = ct_mask_u64_lt(dist1, dist2);

		for (int i = 0; i < half; i++)
		{
			y[i] = ct_select_complex(out1[i], out2[i], use_out1);
		}
	}
}

#if (SCLOUDPLUS_MESSAGE_CODE == SCLOUDPLUS_MESSAGE_CODE_BW32) || \
	(SCLOUDPLUS_MESSAGE_CODE == SCLOUDPLUS_MESSAGE_CODE_BW128) || \
	(SCLOUDPLUS_MESSAGE_CODE == SCLOUDPLUS_MESSAGE_CODE_RBW128)
/* Begin compile-time codec branch. */
/* BW/RBW codec.  The active parameters decide BW32/BW128/RBW128 and blocks. */


static inline unsigned int chunk_bits(size_t index)
{
	unsigned int bits =
		(unsigned int)(2 * scloudplus_tau - popcount_u32((uint32_t)index));

#if (SCLOUDPLUS_MESSAGE_CODE == SCLOUDPLUS_MESSAGE_CODE_RBW128)
	bits -= 1U;
#endif
	return bits;
}

static void apply_bw_layers(Complex *v)
{
	for (size_t stride = 1; stride < scloudplus_bw_complex; stride <<= 1U)
	{
		for (size_t base = 0; base < scloudplus_bw_complex; base += (stride << 1U))
		{
			for (size_t offset = 0; offset < stride; offset++)
			{
				v[base + stride + offset] =
					complex_add(v[base + offset],
								complex_mul_phi(v[base + stride + offset]));
			}
		}
	}
}

static void undo_bw_layers(Complex *w)
{
	for (size_t stride = scloudplus_bw_complex >> 1U; stride > 0U;
		 stride >>= 1U)
	{
		for (size_t base = 0; base < scloudplus_bw_complex; base += (stride << 1U))
		{
			for (size_t offset = 0; offset < stride; offset++)
			{
				w[base + stride + offset] =
					complex_mul_phi_inv(complex_sub(w[base + stride + offset],
														 w[base + offset]));
			}
		}
	}
}

static void reduce_labels_mod_2tau(Complex *v)
{
	for (size_t i = 0; i < scloudplus_bw_complex; i++)
	{
		v[i].real = (int32_t)mod_pow2(v[i].real, scloudplus_tau);
		v[i].imag = (int32_t)mod_pow2(v[i].imag, scloudplus_tau);
	}
}

static void adjust_delabel_output(Complex *v, unsigned int scale_bits)
{
	for (size_t j = 0; j < scloudplus_bw_complex; j++)
	{
		const unsigned int bits = chunk_bits(j);
		const unsigned int real_mod_bits = real_bits(bits);
		const unsigned int imag_mod_bits = bits - real_mod_bits;
		const int32_t imag_value =
			scaled_round_to_int(v[j].imag, scale_bits);
		const int32_t imag_reduced =
			(int32_t)mod_pow2(imag_value, imag_mod_bits);
		const int32_t real_reduced =
			(int32_t)mod_pow2(scaled_round_to_int(v[j].real, scale_bits) -
								  (imag_value - imag_reduced),
							  real_mod_bits);

		v[j].real = real_reduced;
		v[j].imag = imag_reduced;
	}
}

static void load_message_block(const uint8_t *msg, size_t start_bit, Complex *v)
{
	size_t bit_offset = start_bit;

	for (size_t j = 0; j < scloudplus_bw_complex; j++)
	{
		const unsigned int bits = chunk_bits(j);
		const unsigned int re_bits = real_bits(bits);
		const uint32_t packed =
			pqc_load_bits_le(msg, bit_offset, bits);

		v[j].real = (int32_t)(packed & mask_bits(re_bits));
		v[j].imag = (int32_t)(packed >> re_bits);
		bit_offset += bits;
	}
}

static void store_message_block(const Complex *v, uint8_t *msg, size_t start_bit)
{
	size_t bit_offset = start_bit;

	for (size_t j = 0; j < scloudplus_bw_complex; j++)
	{
		const unsigned int bits = chunk_bits(j);
		const unsigned int re_bits = real_bits(bits);
		const uint32_t packed =
			((uint32_t)v[j].imag << re_bits) | (uint32_t)v[j].real;

		pqc_store_bits_le(msg, bit_offset, packed, bits);
		bit_offset += bits;
	}
}

static void write_matrix_block(const Complex *w, uint16_t *matrix,
							   size_t position_offset)
{
	const unsigned int shift =
		(unsigned int)(scloudplus_logq - scloudplus_msg_scale_bits);

	/*
	 * Message-matrix layout: store each complex coordinate as an adjacent
	 * public pair (Re_0, Im_0, Re_1, Im_1, ...).  MsgDec uses the exact same
	 * mapping before invoking BDD.
	 */
	for (size_t i = 0; i < scloudplus_bw_complex; i++)
	{
		const size_t real_idx = position_offset + (i << 1U);
		const size_t imag_idx = real_idx + 1U;

		matrix[real_idx] =
			pqc_mod_pow2(((uint32_t)w[i].real) << shift, scloudplus_logq);
		matrix[imag_idx] =
			pqc_mod_pow2(((uint32_t)w[i].imag) << shift, scloudplus_logq);
	}
}

static inline int32_t center_q_coefficient(uint16_t value)
{
	const uint32_t reduced = (uint32_t)(value & scloudplus_q_mask);
	const uint32_t subtract =
		(reduced >> (scloudplus_logq - 1U)) << scloudplus_logq;

	return (int32_t)reduced - (int32_t)subtract;
}

static void read_matrix_block(const uint16_t *matrix, Complex *w,
							  size_t position_offset)
{
	for (size_t i = 0; i < scloudplus_bw_complex; i++)
	{
		const size_t real_idx = position_offset + (i << 1U);
		const size_t imag_idx = real_idx + 1U;

		w[i] = (Complex){center_q_coefficient(matrix[real_idx]),
						  center_q_coefficient(matrix[imag_idx])};
	}
}

void msg_encode(const uint8_t *msg, uint16_t *matrix)
{
	memset(matrix, 0,
		   (size_t)scloudplus_mbar * scloudplus_nbar * sizeof(uint16_t));

	for (size_t block = 0; block < mat_subm; block++)
	{
		Complex v[scloudplus_bw_complex];

		load_message_block(msg, block * scloudplus_mu, v);
		apply_bw_layers(v);
		reduce_labels_mod_2tau(v);
#if (SCLOUDPLUS_MESSAGE_CODE == SCLOUDPLUS_MESSAGE_CODE_RBW128)
		for (size_t i = 0; i < scloudplus_bw_complex; i++)
		{
			v[i] = complex_mul_phi(v[i]);
		}
#endif
		write_matrix_block(v, matrix, block * scloudplus_bw_n);
	}
}

void msg_decode(const uint16_t *matrix, uint8_t *msg)
{
	memset(msg, 0, scloudplus_ss);

	for (size_t block = 0; block < mat_subm; block++)
	{
		Complex target[scloudplus_bw_complex];
		Complex lattice[scloudplus_bw_complex];
		unsigned int scale_bits =
			(unsigned int)(scloudplus_logq - scloudplus_msg_scale_bits);

		read_matrix_block(matrix, target, block * scloudplus_bw_n);
#if (SCLOUDPLUS_MESSAGE_CODE == SCLOUDPLUS_MESSAGE_CODE_RBW128)
		for (size_t i = 0; i < scloudplus_bw_complex; i++)
		{
			target[i] = complex_mul_one_minus_i(target[i]);
		}
		scale_bits += 1U;
#endif
		bddbwn(target, lattice, (int)scloudplus_bw_n, scale_bits);
		undo_bw_layers(lattice);
		adjust_delabel_output(lattice, scale_bits);
		store_message_block(lattice, msg, block * scloudplus_mu);
	}
}
/* End compile-time codec branch. */
#else
#error "No implementation variant selected for common/encode.c"
#endif
