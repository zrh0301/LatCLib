/**
 * @file matrix_reference.c
 * @brief Direct reference matrix arithmetic for Scloud+.
 *
 * The reference backend deliberately uses a direct row-at-a-time structure:
 * generate one public row of A, use it in the visible matrix loops, then move
 * to the next row. Optimized builds keep their own batched/vectorized kernels.
 */

#include "hash.h"
#include "matrix.h"
#include "gf/mod_pow2.h"
#include "scloudplus_param_common.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#if !defined(SCLOUDPLUS_REF_FAMILY_AES) && \
	!defined(SCLOUDPLUS_REF_FAMILY_SHAKE) && \
	!defined(SCLOUDPLUS_REF_FAMILY_SM3)
#error "matrix_reference.c requires one SCLOUDPLUS_REF_FAMILY_* macro"
#endif

#if (defined(SCLOUDPLUS_REF_FAMILY_AES) + \
	 defined(SCLOUDPLUS_REF_FAMILY_SHAKE) + \
	 defined(SCLOUDPLUS_REF_FAMILY_SM3)) != 1
#error "matrix_reference.c can select only one family"
#endif

#if defined(SCLOUDPLUS_REF_FAMILY_AES)
#include "aes.h"
#endif

typedef struct
{
#if defined(SCLOUDPLUS_REF_FAMILY_AES)
	uint8_t aes_schedule[16 * 11];
#else
	uint8_t seedA[scloudplus_seedA_bytes];
#endif
} ref_a_generator;

static uint16_t mul16(uint16_t a, uint16_t b)
{
	return (uint16_t)((uint32_t)a * (uint32_t)b);
}

static void ref_a_generator_init(ref_a_generator *gen, const uint8_t *seedA)
{
#if defined(SCLOUDPLUS_REF_FAMILY_AES)
	AES128_load_schedule(seedA, gen->aes_schedule);
#else
	memcpy(gen->seedA, seedA, scloudplus_seedA_bytes);
#endif
}

static void ref_a_generator_clear(ref_a_generator *gen)
{
#if defined(SCLOUDPLUS_REF_FAMILY_AES)
	AES128_free_schedule(gen->aes_schedule);
#else
	memset(gen, 0, sizeof(*gen));
#endif
}

static void unpack10_x4(const uint8_t in[5], uint16_t out[4])
{
	out[0] = (uint16_t)(in[0] | ((uint16_t)(in[1] & 0x03U) << 8U));
	out[1] = (uint16_t)((in[1] >> 2U) | ((uint16_t)(in[2] & 0x0FU) << 6U));
	out[2] = (uint16_t)((in[2] >> 4U) | ((uint16_t)(in[3] & 0x3FU) << 4U));
	out[3] = (uint16_t)((in[3] >> 6U) | ((uint16_t)in[4] << 2U));
}

static void unpack_a_row(const uint8_t *packed, uint16_t *a_row)
{
	for (size_t i = 0U, j = 0U; j < (size_t)scloudplus_n; i += 5U, j += 4U)
	{
		unpack10_x4(packed + i, a_row + j);
	}
}

static void generate_a_row(ref_a_generator *gen, size_t row, uint16_t *a_row)
{
	uint8_t packed[scloudplus_a_padded_bytes_per_row];

#if defined(SCLOUDPLUS_REF_FAMILY_AES)
	AES128_CTR_zero_sch((uint32_t)(row * (size_t)scloudplus_a_blocks_per_row),
						(size_t)scloudplus_a_blocks_per_row,
						gen->aes_schedule, packed);
#elif defined(SCLOUDPLUS_REF_FAMILY_SHAKE)
	uint8_t input[scloudplus_seedA_bytes + 4U];

	memcpy(input, gen->seedA, scloudplus_seedA_bytes);
	input[scloudplus_seedA_bytes + 0U] = (uint8_t)(row >> 0U);
	input[scloudplus_seedA_bytes + 1U] = (uint8_t)(row >> 8U);
	input[scloudplus_seedA_bytes + 2U] = (uint8_t)(row >> 16U);
	input[scloudplus_seedA_bytes + 3U] = (uint8_t)(row >> 24U);
	shake128(packed, scloudplus_a_bytes_per_row, input, sizeof(input));
	memset(input, 0, sizeof(input));
#else
	sm3_bounded_xof_seed_row(packed, scloudplus_a_bytes_per_row,
							 gen->seedA, (uint32_t)row);
#endif
	unpack_a_row(packed, a_row);
	memset(packed, 0, sizeof(packed));
}

void mat_add(uint16_t *lhs, uint16_t *rhs, int len, uint16_t *out)
{
	for (int i = 0; i < len; i++)
	{
		out[i] = pqc_mod_pow2((uint32_t)lhs[i] + (uint32_t)rhs[i], scloudplus_logq);
	}
}

void mat_sub(uint16_t *lhs, uint16_t *rhs, int len, uint16_t *out)
{
	for (int i = 0; i < len; i++)
	{
		out[i] = pqc_mod_pow2((uint32_t)lhs[i] - (uint32_t)rhs[i], scloudplus_logq);
	}
}

void mul_as_e(const uint8_t *seedA, const uint16_t *S, const uint16_t *E,
			  uint16_t *B)
{
	ref_a_generator gen;
	uint16_t a_row[scloudplus_n];

	memcpy(B, E, (size_t)scloudplus_m * scloudplus_nbar * sizeof(uint16_t));
	ref_a_generator_init(&gen, seedA);
	for (int i = 0; i < scloudplus_m; i++)
	{
		generate_a_row(&gen, (size_t)i, a_row);
		for (int k = 0; k < scloudplus_nbar; k++)
		{
			uint16_t sum = 0;
			const uint16_t *s_col = S + (size_t)k * scloudplus_n;

			for (int j = 0; j < scloudplus_n; j++)
			{
				sum = (uint16_t)(sum + mul16(a_row[j], s_col[j]));
			}
			B[(size_t)i * scloudplus_nbar + (size_t)k] =
				(uint16_t)(B[(size_t)i * scloudplus_nbar + (size_t)k] + sum);
		}
	}
	ref_a_generator_clear(&gen);
	memset(a_row, 0, sizeof(a_row));
}

void mul_sa_e(const uint8_t *seedA, const uint16_t *S, uint16_t *E,
			  uint16_t *C)
{
	ref_a_generator gen;
	uint16_t a_row[scloudplus_n];

	memcpy(C, E, (size_t)scloudplus_mbar * scloudplus_n * sizeof(uint16_t));
	ref_a_generator_init(&gen, seedA);
	for (int i = 0; i < scloudplus_m; i++)
	{
		generate_a_row(&gen, (size_t)i, a_row);
		for (int p = 0; p < scloudplus_mbar; p++)
		{
			const uint16_t s = S[(size_t)p * scloudplus_m + (size_t)i];
			uint16_t *c_row = C + (size_t)p * scloudplus_n;

			for (int q = 0; q < scloudplus_n; q++)
			{
				c_row[q] = (uint16_t)(c_row[q] + mul16(s, a_row[q]));
			}
		}
	}
	ref_a_generator_clear(&gen);
	memset(a_row, 0, sizeof(a_row));
}

void mul_sb_e(const uint16_t *S, const uint16_t *B, const uint16_t *E,
			  uint16_t *out)
{
	memcpy(out, E,
		   (size_t)scloudplus_mbar * scloudplus_nbar * sizeof(uint16_t));
	for (int i = 0; i < scloudplus_mbar; i++)
	{
		const uint16_t *s_row = S + (size_t)i * scloudplus_m;
		uint16_t *out_row = out + (size_t)i * scloudplus_nbar;

		for (int k = 0; k < scloudplus_m; k++)
		{
			const uint16_t s = s_row[k];
			const uint16_t *b_row = B + (size_t)k * scloudplus_nbar;

			for (int j = 0; j < scloudplus_nbar; j++)
			{
				out_row[j] = (uint16_t)(out_row[j] + mul16(s, b_row[j]));
			}
		}
	}
}

void mul_cs(uint16_t *C, uint16_t *S, uint16_t *out)
{
	memset(out, 0,
		   (size_t)scloudplus_mbar * scloudplus_nbar * sizeof(uint16_t));
	for (int i = 0; i < scloudplus_mbar; i++)
	{
		uint16_t *out_row = out + (size_t)i * scloudplus_nbar;

		for (int k = 0; k < scloudplus_n; k++)
		{
			const uint16_t c = C[(size_t)i * scloudplus_n + (size_t)k];

			for (int j = 0; j < scloudplus_nbar; j++)
			{
				out_row[j] = (uint16_t)(out_row[j] +
					mul16(c, S[(size_t)j * scloudplus_n + (size_t)k]));
			}
		}
	}
}
