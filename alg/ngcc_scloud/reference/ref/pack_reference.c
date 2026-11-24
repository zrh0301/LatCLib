/**
 * @file pack_reference.c
 * @brief Scalar byte serialization for Scloud+ keys and ciphertext components.
 *
 * The public key matrix B and ciphertext matrices C1/C2 are elements of
 * Z/qZ with q = 2^10.  They are serialized in row-major order using the
 * algorithm-level 10-bit little-endian packing: every four coefficients
 * occupy five bytes.  This file is therefore part of the byte-level
 * interoperability contract for pk and ct.
 *
 * The secret matrix S is different: its coefficients are sampled from
 * {-1,0,1}.  The implementation stores S with ShortPack, a lossless two-bit
 * two's-complement encoding: 0 -> 00, 1 -> 01, and -1 -> 11.  This keeps the
 * KEM private key short while preserving the exact mathematical secret matrix.
 *
 * The PACK_LINEAR_Q10_FIXED and UNPACK_LINEAR_Q10_FIXED macros deliberately
 * receive compile-time counts at their call sites.  Earlier generic helpers
 * used a run-time count parameter and compiled poorly on some platforms; the
 * fixed-count form gives the compiler enough information to specialize the
 * pk, C1, and C2 loops without changing the serialized format.
 */

#include "pack.h"
#include "scloudplus_param_common.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/**
 * @brief Pack four q-coefficients into five little-endian bytes.
 *
 * Coefficient a occupies bits [0,9], b occupies [10,19], c occupies [20,29],
 * and d occupies [30,39].  Only the low 10 bits of each input are serialized.
 */
static void pack4_q10(uint16_t a, uint16_t b, uint16_t c, uint16_t d,
					  uint8_t out[5])
{
	const uint64_t w = ((uint64_t)(a & 0x03ffU) << 0U) |
					   ((uint64_t)(b & 0x03ffU) << 10U) |
					   ((uint64_t)(c & 0x03ffU) << 20U) |
					   ((uint64_t)(d & 0x03ffU) << 30U);

	out[0] = (uint8_t)(w >> 0U);
	out[1] = (uint8_t)(w >> 8U);
	out[2] = (uint8_t)(w >> 16U);
	out[3] = (uint8_t)(w >> 24U);
	out[4] = (uint8_t)(w >> 32U);
}

/**
 * @brief Inverse of pack4_q10 for one five-byte group.
 */
static void unpack4_q10(const uint8_t in[5], uint16_t *a, uint16_t *b,
						uint16_t *c, uint16_t *d)
{
	const uint64_t w = ((uint64_t)in[0] << 0U) |
					   ((uint64_t)in[1] << 8U) |
					   ((uint64_t)in[2] << 16U) |
					   ((uint64_t)in[3] << 24U) |
					   ((uint64_t)in[4] << 32U);

	*a = (uint16_t)((w >> 0U) & 0x03ffU);
	*b = (uint16_t)((w >> 10U) & 0x03ffU);
	*c = (uint16_t)((w >> 20U) & 0x03ffU);
	*d = (uint16_t)((w >> 30U) & 0x03ffU);
}

/**
 * @brief Fixed-count row-major Pack10 loop.
 *
 * COUNT is always a public parameter expression such as m*nbar or mbar*n.
 * All supported dimensions are multiples of four, so the loop has no tail
 * case and the byte length is exactly 5*COUNT/4.
 */
#define PACK_LINEAR_Q10_FIXED(IN, COUNT, OUT)                           \
	do                                                                  \
	{                                                                   \
		const uint16_t *pack_in_ = (IN);                                \
		uint8_t *pack_out_ = (OUT);                                     \
		for (size_t pack_i_ = 0U; pack_i_ < (size_t)(COUNT);            \
			 pack_i_ += 4U)                                             \
		{                                                               \
			pack4_q10(pack_in_[pack_i_ + 0U], pack_in_[pack_i_ + 1U],   \
					  pack_in_[pack_i_ + 2U], pack_in_[pack_i_ + 3U],   \
					  pack_out_);                                       \
			pack_out_ += 5U;                                            \
		}                                                               \
	} while (0)

/**
 * @brief Fixed-count row-major Unpack10 loop.
 */
#define UNPACK_LINEAR_Q10_FIXED(IN, COUNT, OUT)                         \
	do                                                                  \
	{                                                                   \
		const uint8_t *unpack_in_ = (IN);                               \
		uint16_t *unpack_out_ = (OUT);                                  \
		for (size_t unpack_i_ = 0U; unpack_i_ < (size_t)(COUNT);        \
			 unpack_i_ += 4U)                                           \
		{                                                               \
			unpack4_q10(unpack_in_, &unpack_out_[unpack_i_ + 0U],       \
						&unpack_out_[unpack_i_ + 1U],                  \
						&unpack_out_[unpack_i_ + 2U],                  \
						&unpack_out_[unpack_i_ + 3U]);                 \
			unpack_in_ += 5U;                                           \
		}                                                               \
	} while (0)

/**
 * @brief Reduce a public matrix buffer to q by keeping the low 10 bits.
 */
static void mod_q_vec(const uint16_t *in, size_t count, uint16_t *out)
{
	for (size_t i = 0U; i < count; i++)
	{
		out[i] = (uint16_t)(in[i] & scloudplus_q_mask);
	}
}

/**
 * @brief ShortPack for the secret matrix S.
 *
 * S is stored row-major.  Each coefficient is interpreted as a signed ternary
 * integer and serialized by its two low bits.  Thus -1 is encoded as 3
 * because the in-memory uint16_t representation of -1 is 0xffff.  Four
 * coefficients occupy one byte.  This encoding is private-key serialization
 * only; public q-elements still use Pack10.
 */
void pack_sk(const uint16_t *S, uint8_t *sk)
{
	size_t i = 0U;
	const size_t count = (size_t)scloudplus_n * scloudplus_nbar;

	for (; i < count; i += 4U)
	{
		const uint8_t s0 = (uint8_t)(S[i + 0U] & 0x03U);
		const uint8_t s1 = (uint8_t)(S[i + 1U] & 0x03U);
		const uint8_t s2 = (uint8_t)(S[i + 2U] & 0x03U);
		const uint8_t s3 = (uint8_t)(S[i + 3U] & 0x03U);

		sk[i >> 2U] = (uint8_t)(s0 | (uint8_t)(s1 << 2U) |
						  (uint8_t)(s2 << 4U) | (uint8_t)(s3 << 6U));
	}
}

static inline uint16_t unpack_sk_coeff(uint8_t value)
{
	const uint16_t x = (uint16_t)(value & 0x03U);
	const uint16_t sign = (uint16_t)(x >> 1U);
	const uint16_t sign_mask = (uint16_t)(0U - sign);

	return (uint16_t)(x | (sign_mask & 0xfffcU));
}

/**
 * @brief Inverse of pack_sk.
 */
void unpack_sk(const uint8_t *sk, uint16_t *S)
{
	size_t i = 0U;
	const size_t count = (size_t)scloudplus_n * scloudplus_nbar;

	for (; i < count; i += 4U)
	{
		const uint8_t packed = sk[i >> 2U];

		S[i + 0U] = unpack_sk_coeff(packed);
		S[i + 1U] = unpack_sk_coeff((uint8_t)(packed >> 2U));
		S[i + 2U] = unpack_sk_coeff((uint8_t)(packed >> 4U));
		S[i + 3U] = unpack_sk_coeff((uint8_t)(packed >> 6U));
	}
}

/**
 * @brief Serialize the public matrix B in pk after seedA.
 */
void pack_pk(const uint16_t *B, uint8_t *pk)
{
	PACK_LINEAR_Q10_FIXED(B, (size_t)scloudplus_m * scloudplus_nbar, pk);
}

/**
 * @brief Parse the public matrix B from pk.
 */
void unpack_pk(const uint8_t *pk, uint16_t *B)
{
	UNPACK_LINEAR_Q10_FIXED(pk, (size_t)scloudplus_m * scloudplus_nbar, B);
}

/**
 * @brief Canonicalize C1 coefficients before ciphertext packing.
 */
void reduce_c1(const uint16_t *C, uint16_t *out)
{
	mod_q_vec(C, (size_t)scloudplus_mbar * scloudplus_n, out);
}

/**
 * @brief Canonicalize C2 coefficients before ciphertext packing.
 */
void reduce_c2(const uint16_t *C, uint16_t *out)
{
	mod_q_vec(C, (size_t)scloudplus_mbar * scloudplus_nbar, out);
}


/**
 * @brief Serialize C1 in row-major Pack10 order.
 */
void pack_c1(const uint16_t *C, uint8_t *out)
{
	PACK_LINEAR_Q10_FIXED(C, (size_t)scloudplus_mbar * scloudplus_n, out);
}

/**
 * @brief Parse C1 from row-major Pack10 order.
 */
void unpack_c1(const uint8_t *in, uint16_t *C)
{
	UNPACK_LINEAR_Q10_FIXED(in, (size_t)scloudplus_mbar * scloudplus_n, C);
}

/**
 * @brief Serialize C2 in row-major Pack10 order immediately after C1.
 */
void pack_c2(const uint16_t *C, uint8_t *out)
{
	PACK_LINEAR_Q10_FIXED(C, (size_t)scloudplus_mbar * scloudplus_nbar, out);
}

/**
 * @brief Parse C2 from row-major Pack10 order.
 */
void unpack_c2(const uint8_t *in, uint16_t *C)
{
	UNPACK_LINEAR_Q10_FIXED(in, (size_t)scloudplus_mbar * scloudplus_nbar, C);
}
