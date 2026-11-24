/**
 * @file sample.c
 * @brief Unified short-noise sampler for all Scloud+ families and backends.
 *
 * The sampler has two independent choices:
 *
 * - the XOF that supplies the bit stream: SHAKE for AES/SHAKE instances, SM3
 *   bounded SM3 expansion interface for SM3 instances;
 * - the public BD distribution parameters selected by the concrete instance.
 *
 * Both choices are compile-time decisions.  The implementation below keeps one
 * reader abstraction and one implementation of each BD sampler, so the source is
 * reviewed as one algorithm instead of as one file per parameter set.  The
 * preprocessor only removes unused primitive adapters and unused BD kernels; it
 * does not introduce runtime dispatch or change the serialized outputs.
 */

#include "hash.h"
#include "scloudplus_param_common.h"
#include "sample.h"
#include <stdint.h>
#include <string.h>

#if defined(SCLOUDPLUS_SAMPLE_XOF_SM3)
#define SAMPLE_XOF_SM3 1
#else
#define SAMPLE_XOF_SHAKE 1
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD != 2) && \
	(SCLOUDPLUS_SAMPLE_SECRET_BD != 4) && \
	(SCLOUDPLUS_SAMPLE_SECRET_BD != 6) && \
	(SCLOUDPLUS_SAMPLE_SECRET_BD != 12)
#error "Unsupported secret sample distribution"
#endif

#if (SCLOUDPLUS_SAMPLE_ERROR_BD != 2) && \
	(SCLOUDPLUS_SAMPLE_ERROR_BD != 6) && \
	(SCLOUDPLUS_SAMPLE_ERROR_BD != 12)
#error "Unsupported error sample distribution"
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 6) || \
	(SCLOUDPLUS_SAMPLE_ERROR_BD == 6) || \
	(SCLOUDPLUS_SAMPLE_SECRET_BD == 12) || \
	(SCLOUDPLUS_SAMPLE_ERROR_BD == 12)
#define SAMPLE_USE_SHAKE256 1
#else
#define SAMPLE_USE_SHAKE128 1
#endif

#if defined(SAMPLE_XOF_SM3)
#define SAMPLE_SM3_BLOCKS 4U
#define SAMPLE_BLOCK_BYTES (SM3_BOUNDED_XOF_RATE * SAMPLE_SM3_BLOCKS)
#else
#define SAMPLE_SHAKE_BLOCKS 4U
#if defined(SAMPLE_USE_SHAKE256)
#define SAMPLE_BLOCK_BYTES (SHAKE256_RATE * SAMPLE_SHAKE_BLOCKS)
#else
#define SAMPLE_BLOCK_BYTES (SHAKE128_RATE * SAMPLE_SHAKE_BLOCKS)
#endif
#endif

typedef struct
{
#if defined(SAMPLE_XOF_SM3)
	sm3_bounded_xof_state state;
#else
	keccak_state state;
#endif
	uint8_t block[SAMPLE_BLOCK_BYTES];
	size_t block_pos;
	size_t block_len;
} sample_reader;

#if defined(SCLOUDPLUS_BACKEND_AVX2) || defined(SCLOUDPLUS_BACKEND_NEON)
#define SAMPLE_HAS_SIMD_BD 1
#define SAMPLE_SIMD_INPUT_BYTES 32U
void scloudplus_sample_bd2_backend(const uint8_t *bytes, size_t byte_count,
								   uint16_t *out);
void scloudplus_sample_bd4_backend(const uint8_t *bytes, size_t byte_count,
								   uint16_t *out);
#endif

#if defined(SCLOUDPLUS_BACKEND_AVX2) || defined(SCLOUDPLUS_BACKEND_NEON)
#define SAMPLE_HAS_REJECT_BD_BACKEND 1
size_t scloudplus_sample_bd6_bits_backend(const uint8_t *bytes,
										  size_t byte_count, uint8_t *bits);
size_t scloudplus_sample_bd12_bits_backend(const uint8_t *bytes,
										   size_t byte_count, uint8_t *bits);
#endif

static void reader_refill(sample_reader *r)
{
	r->block_len = sizeof(r->block);
#if defined(SAMPLE_XOF_SM3)
	sm3_bounded_xof_squeezeblocks(r->block, SAMPLE_SM3_BLOCKS, &r->state);
#elif defined(SAMPLE_USE_SHAKE256)
	shake256_squeezeblocks(r->block, SAMPLE_SHAKE_BLOCKS, &r->state);
#else
	shake128_squeezeblocks(r->block, SAMPLE_SHAKE_BLOCKS, &r->state);
#endif
	r->block_pos = 0U;
}

static void reader_init(sample_reader *r, const uint8_t *seed, size_t seedlen)
{
	memset(r, 0, sizeof(*r));
#if defined(SAMPLE_XOF_SM3)
	sm3_bounded_xof_absorb_once(&r->state, seed, seedlen);
#elif defined(SAMPLE_USE_SHAKE256)
	shake256_absorb_once(&r->state, seed, seedlen);
#else
	shake128_absorb_once(&r->state, seed, seedlen);
#endif
	reader_refill(r);
}

static inline uint8_t reader_byte(sample_reader *r)
{
	if (r->block_pos == r->block_len)
	{
		reader_refill(r);
	}
	return r->block[r->block_pos++];
}

#if defined(SAMPLE_HAS_SIMD_BD) || \
	(SCLOUDPLUS_SAMPLE_SECRET_BD == 6) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 6) || \
	(SCLOUDPLUS_SAMPLE_SECRET_BD == 12) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 12)
static inline const uint8_t *reader_contiguous(sample_reader *r, size_t bytes)
{
	if ((r->block_len - r->block_pos) >= bytes)
	{
		const uint8_t *out = r->block + r->block_pos;

		r->block_pos += bytes;
		return out;
	}
	return 0;
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 6) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 6) || \
	(SCLOUDPLUS_SAMPLE_SECRET_BD == 12) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 12)
static inline const uint8_t *reader_bytes(sample_reader *r, size_t bytes,
										  uint8_t *tmp)
{
	const uint8_t *buf = reader_contiguous(r, bytes);

	if (buf != 0)
	{
		return buf;
	}
	for (size_t i = 0U; i < bytes; i++)
	{
		tmp[i] = reader_byte(r);
	}
	return tmp;
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 4) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 4) || \
	(!defined(SAMPLE_HAS_REJECT_BD_BACKEND) && \
	 ((SCLOUDPLUS_SAMPLE_SECRET_BD == 6) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 6) || \
	  (SCLOUDPLUS_SAMPLE_SECRET_BD == 12) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 12)))
static inline uint16_t ct_is_zero_u32(uint32_t value)
{
	return (uint16_t)((((value | (0U - value)) >> 31U) ^ 1U) & 1U);
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 2) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 2)
#if ((scloudplus_n * scloudplus_nbar) % 4) != 0 || \
	((scloudplus_mbar * scloudplus_m) % 4) != 0 || \
	((scloudplus_m * scloudplus_nbar) % 4) != 0 || \
	((scloudplus_mbar * scloudplus_n) % 4) != 0 || \
	((scloudplus_mbar * scloudplus_nbar) % 4) != 0
#error "BD2 sampler assumes current Scloud+ sample sizes are multiples of 4"
#endif

static inline uint16_t bd2_coeff_from_pair(uint8_t pair)
{
	const uint16_t lhs = (uint16_t)((pair >> 1U) & 1U);
	const uint16_t rhs = (uint16_t)(pair & 1U);

	return (uint16_t)(lhs - rhs);
}

static void sample_bd2_from_reader(sample_reader *r, size_t coeffs,
								   uint16_t *out)
{
	size_t i = 0U;

#if defined(SAMPLE_HAS_SIMD_BD)
	while ((i + 4U * SAMPLE_SIMD_INPUT_BYTES) <= coeffs)
	{
		const uint8_t *bytes = reader_contiguous(r, SAMPLE_SIMD_INPUT_BYTES);

		if (bytes != 0)
		{
			scloudplus_sample_bd2_backend(bytes, SAMPLE_SIMD_INPUT_BYTES, out + i);
			i += 4U * SAMPLE_SIMD_INPUT_BYTES;
			continue;
		}

		const uint8_t byte = reader_byte(r);

		out[i + 0U] = bd2_coeff_from_pair((uint8_t)(byte & 0x03U));
		out[i + 1U] = bd2_coeff_from_pair((uint8_t)((byte >> 2U) & 0x03U));
		out[i + 2U] = bd2_coeff_from_pair((uint8_t)((byte >> 4U) & 0x03U));
		out[i + 3U] = bd2_coeff_from_pair((uint8_t)(byte >> 6U));
		i += 4U;
	}
#endif

	for (; i < coeffs; i += 4U)
	{
		const uint8_t byte = reader_byte(r);

		out[i + 0U] = bd2_coeff_from_pair((uint8_t)(byte & 0x03U));
		out[i + 1U] = bd2_coeff_from_pair((uint8_t)((byte >> 2U) & 0x03U));
		out[i + 2U] = bd2_coeff_from_pair((uint8_t)((byte >> 4U) & 0x03U));
		out[i + 3U] = bd2_coeff_from_pair((uint8_t)(byte >> 6U));
	}
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 4) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 4)
#if ((scloudplus_n * scloudplus_nbar) % 2) != 0 || \
	((scloudplus_mbar * scloudplus_m) % 2) != 0 || \
	((scloudplus_m * scloudplus_nbar) % 2) != 0 || \
	((scloudplus_mbar * scloudplus_n) % 2) != 0 || \
	((scloudplus_mbar * scloudplus_nbar) % 2) != 0
#error "BD4 sampler assumes current Scloud+ sample sizes are multiples of 2"
#endif

static inline uint16_t bd4_coeff_from_nibble(uint8_t nibble)
{
	const uint16_t lhs = ct_is_zero_u32(nibble & 0x03U);
	const uint16_t rhs = ct_is_zero_u32((uint32_t)(nibble >> 2U) & 0x03U);

	return (uint16_t)(lhs - rhs);
}

static void sample_bd4_from_reader(sample_reader *r, size_t coeffs,
								   uint16_t *out)
{
	size_t i = 0U;

#if defined(SAMPLE_HAS_SIMD_BD)
	while ((i + 2U * SAMPLE_SIMD_INPUT_BYTES) <= coeffs)
	{
		const uint8_t *bytes = reader_contiguous(r, SAMPLE_SIMD_INPUT_BYTES);

		if (bytes != 0)
		{
			scloudplus_sample_bd4_backend(bytes, SAMPLE_SIMD_INPUT_BYTES, out + i);
			i += 2U * SAMPLE_SIMD_INPUT_BYTES;
			continue;
		}

		const uint8_t byte = reader_byte(r);

		out[i + 0U] = bd4_coeff_from_nibble((uint8_t)(byte & 0x0FU));
		out[i + 1U] = bd4_coeff_from_nibble((uint8_t)(byte >> 4U));
		i += 2U;
	}
#endif

	for (; i < coeffs; i += 2U)
	{
		const uint8_t byte = reader_byte(r);

		out[i + 0U] = bd4_coeff_from_nibble((uint8_t)(byte & 0x0FU));
		out[i + 1U] = bd4_coeff_from_nibble((uint8_t)(byte >> 4U));
	}
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 6) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 6) || \
	(SCLOUDPLUS_SAMPLE_SECRET_BD == 12) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 12)
#define SAMPLE_REJECT_BITS_CAPACITY 256U

typedef struct
{
	sample_reader reader;
	uint8_t bits[SAMPLE_REJECT_BITS_CAPACITY];
	size_t bit_pos;
	size_t bit_len;
} sample_reject_reader;

#if !defined(SAMPLE_HAS_REJECT_BD_BACKEND)
static inline uint8_t sample_accept_lt_u32(uint32_t candidate, uint32_t limit)
{
	return (uint8_t)((candidate - limit) >> 31U);
}
#endif

static inline uint16_t sample_coeff_from_bits(uint8_t lhs, uint8_t rhs)
{
	return (uint16_t)((int)lhs - (int)rhs);
}

static void sample_reject_reader_init(sample_reject_reader *reader,
									  const uint8_t *seed, size_t seedlen)
{
	reader_init(&reader->reader, seed, seedlen);
	reader->bit_pos = 0U;
	reader->bit_len = 0U;
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 6) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 6)
#define SAMPLE_BD6_BATCH_BYTES 96U
#if ((SAMPLE_BD6_BATCH_BYTES * 8U) / 3U) > SAMPLE_REJECT_BITS_CAPACITY
#error "BD6 reject batch exceeds accepted-bit queue"
#endif

#if !defined(SAMPLE_HAS_REJECT_BD_BACKEND)
static size_t sample_bd6_bits_scalar(const uint8_t *bytes, size_t byte_count,
									 uint8_t *bits);
#endif

static void sample_bd6_refill(sample_reject_reader *reader)
{
	do
	{
		uint8_t tmp[SAMPLE_BD6_BATCH_BYTES];
		const uint8_t *bytes = reader_bytes(&reader->reader,
											SAMPLE_BD6_BATCH_BYTES, tmp);

#if defined(SAMPLE_HAS_REJECT_BD_BACKEND)
		reader->bit_len = scloudplus_sample_bd6_bits_backend(
			bytes, SAMPLE_BD6_BATCH_BYTES, reader->bits);
#else
		reader->bit_len = sample_bd6_bits_scalar(bytes, SAMPLE_BD6_BATCH_BYTES,
												 reader->bits);
#endif
		reader->bit_pos = 0U;
	} while (reader->bit_len == 0U);
}

#if !defined(SAMPLE_HAS_REJECT_BD_BACKEND)
static size_t sample_bd6_bits_scalar(const uint8_t *bytes, size_t byte_count,
									 uint8_t *bits)
{
	size_t count = 0U;

	for (size_t i = 0U; i < byte_count; i += 3U)
	{
		const uint32_t chunk = (uint32_t)bytes[i + 0U] |
							   ((uint32_t)bytes[i + 1U] << 8U) |
							   ((uint32_t)bytes[i + 2U] << 16U);

		for (unsigned int lane = 0U; lane < 8U; lane++)
		{
			const uint32_t candidate = (chunk >> (3U * lane)) & 0x07U;

			bits[count] = (uint8_t)ct_is_zero_u32(candidate);
			count += sample_accept_lt_u32(candidate, 6U);
		}
	}
	return count;
}
#endif

static inline uint8_t sample_bd6_next_bit(sample_reject_reader *reader)
{
	if (reader->bit_pos == reader->bit_len)
	{
		sample_bd6_refill(reader);
	}
	return reader->bits[reader->bit_pos++];
}

static void sample_bd6_from_reader(sample_reject_reader *reader, size_t coeffs,
								   uint16_t *out)
{
	for (size_t i = 0U; i < coeffs; i++)
	{
		const uint8_t lhs = sample_bd6_next_bit(reader);
		const uint8_t rhs = sample_bd6_next_bit(reader);

		out[i] = sample_coeff_from_bits(lhs, rhs);
	}
}
#endif

#if (SCLOUDPLUS_SAMPLE_SECRET_BD == 12) || (SCLOUDPLUS_SAMPLE_ERROR_BD == 12)
#define SAMPLE_BD12_BATCH_BYTES 128U
#if (SAMPLE_BD12_BATCH_BYTES * 2U) > SAMPLE_REJECT_BITS_CAPACITY
#error "BD12 reject batch exceeds accepted-bit queue"
#endif

#if !defined(SAMPLE_HAS_REJECT_BD_BACKEND)
static size_t sample_bd12_bits_scalar(const uint8_t *bytes, size_t byte_count,
									  uint8_t *bits);
#endif

static void sample_bd12_refill(sample_reject_reader *reader)
{
	do
	{
		uint8_t tmp[SAMPLE_BD12_BATCH_BYTES];
		const uint8_t *bytes = reader_bytes(&reader->reader,
											SAMPLE_BD12_BATCH_BYTES, tmp);

#if defined(SAMPLE_HAS_REJECT_BD_BACKEND)
		reader->bit_len = scloudplus_sample_bd12_bits_backend(
			bytes, SAMPLE_BD12_BATCH_BYTES, reader->bits);
#else
		reader->bit_len = sample_bd12_bits_scalar(bytes, SAMPLE_BD12_BATCH_BYTES,
												  reader->bits);
#endif
		reader->bit_pos = 0U;
	} while (reader->bit_len == 0U);
}

#if !defined(SAMPLE_HAS_REJECT_BD_BACKEND)
static size_t sample_bd12_bits_scalar(const uint8_t *bytes, size_t byte_count,
									  uint8_t *bits)
{
	size_t count = 0U;

	for (size_t i = 0U; i < byte_count; i++)
	{
		const uint32_t lo = bytes[i] & 0x0FU;
		const uint32_t hi = bytes[i] >> 4U;

		bits[count] = (uint8_t)ct_is_zero_u32(lo);
		count += sample_accept_lt_u32(lo, 12U);
		bits[count] = (uint8_t)ct_is_zero_u32(hi);
		count += sample_accept_lt_u32(hi, 12U);
	}
	return count;
}
#endif

static inline uint8_t sample_bd12_next_bit(sample_reject_reader *reader)
{
	if (reader->bit_pos == reader->bit_len)
	{
		sample_bd12_refill(reader);
	}
	return reader->bits[reader->bit_pos++];
}

static void sample_bd12_from_reader(sample_reject_reader *reader, size_t coeffs,
									uint16_t *out)
{
	for (size_t i = 0U; i < coeffs; i++)
	{
		const uint8_t lhs = sample_bd12_next_bit(reader);
		const uint8_t rhs = sample_bd12_next_bit(reader);

		out[i] = sample_coeff_from_bits(lhs, rhs);
	}
}
#endif

static void sample_secret_from_reader(const uint8_t *seed, size_t seedlen,
									  size_t coeffs, uint16_t *out)
{
#if SCLOUDPLUS_SAMPLE_SECRET_BD == 2
	sample_reader r;
	reader_init(&r, seed, seedlen);
	sample_bd2_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#elif SCLOUDPLUS_SAMPLE_SECRET_BD == 4
	sample_reader r;
	reader_init(&r, seed, seedlen);
	sample_bd4_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#elif SCLOUDPLUS_SAMPLE_SECRET_BD == 6
	sample_reject_reader r;
	sample_reject_reader_init(&r, seed, seedlen);
	sample_bd6_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#else
	sample_reject_reader r;
	sample_reject_reader_init(&r, seed, seedlen);
	sample_bd12_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#endif
}

static void sample_error_from_reader(const uint8_t *seed, size_t seedlen,
									 size_t coeffs, uint16_t *out)
{
#if SCLOUDPLUS_SAMPLE_ERROR_BD == 2
	sample_reader r;
	reader_init(&r, seed, seedlen);
	sample_bd2_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#elif SCLOUDPLUS_SAMPLE_ERROR_BD == 6
	sample_reject_reader r;
	sample_reject_reader_init(&r, seed, seedlen);
	sample_bd6_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#else
	sample_reject_reader r;
	sample_reject_reader_init(&r, seed, seedlen);
	sample_bd12_from_reader(&r, coeffs, out);
	memset(&r, 0, sizeof(r));
#endif
}

void sample_s(const uint8_t *seed, size_t seedlen, uint16_t *s)
{
	sample_secret_from_reader(seed, seedlen,
							  (size_t)scloudplus_n * scloudplus_nbar, s);
}

void sample_sp(const uint8_t *seed, size_t seedlen, uint16_t *sp)
{
	sample_secret_from_reader(seed, seedlen,
							  (size_t)scloudplus_mbar * scloudplus_m, sp);
}

void sample_e(const uint8_t *seed, size_t seedlen, uint16_t *e)
{
	sample_error_from_reader(seed, seedlen,
							 (size_t)scloudplus_m * scloudplus_nbar, e);
}

void sample_e12(const uint8_t *seed, size_t seedlen, uint16_t *e1, uint16_t *e2)
{
#if SCLOUDPLUS_SAMPLE_ERROR_BD == 2
	sample_reader r;
	reader_init(&r, seed, seedlen);
	sample_bd2_from_reader(&r, (size_t)scloudplus_mbar * scloudplus_n, e1);
	sample_bd2_from_reader(&r, (size_t)scloudplus_mbar * scloudplus_nbar, e2);
	memset(&r, 0, sizeof(r));
#elif SCLOUDPLUS_SAMPLE_ERROR_BD == 6
	sample_reject_reader r;
	sample_reject_reader_init(&r, seed, seedlen);
	sample_bd6_from_reader(&r, (size_t)scloudplus_mbar * scloudplus_n, e1);
	sample_bd6_from_reader(&r, (size_t)scloudplus_mbar * scloudplus_nbar, e2);
	memset(&r, 0, sizeof(r));
#else
	sample_reject_reader r;
	sample_reject_reader_init(&r, seed, seedlen);
	sample_bd12_from_reader(&r, (size_t)scloudplus_mbar * scloudplus_n, e1);
	sample_bd12_from_reader(&r, (size_t)scloudplus_mbar * scloudplus_nbar, e2);
	memset(&r, 0, sizeof(r));
#endif
}
