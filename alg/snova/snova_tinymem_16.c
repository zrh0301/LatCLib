// SPDX-License-Identifier: MIT

/**
 * Reduced memory implementation.
 *
 * Copyright (c) 2026 SNOVA TEAM
 */

#include <stdint.h>
#include <string.h>

#include "snova.h"
#include "symmetric_tiny.h"
#include "abq_snova.h"

// Flag to control use of the public XOF bytes. Setting it results in KAT files that are different from the Round 3 specification.
#ifndef PUB_REORDER
#define PUB_REORDER 0
#endif

// Flag to use nibbled whipped_sig. Not using it is much faster.
#ifndef NIBBLE_UINT8_T
#define NIBBLE_UINT8_T 1
#endif

#if NIBBLE_UINT8_T
// Ensure that SNOVA_olr is even.
#define SNOVA_olr (((SNOVA_o * SNOVA_lr) / 2 + 1) * 2)
#else
#define SNOVA_olr (SNOVA_o * SNOVA_lr + 1)
#endif

/**
 * Utilities
 */
#define gf16_expand(a) (((a) | ((a) << 3) | ((a) << 6) | ((a) << 9)) & 0x1111)
#define gf16_compress(a) (((a) & 0xf) ^ (((a) >> 3) & 0xe) ^ (((a) >> 6) & 0xe) ^ (((a) >> 7) & 0x1) ^ (((a) >> 9) & 0xe) ^ (((a) >> 10) & 0x3) ^ (((a) >> 12) & 0xe) ^ (((a) >> 13) & 0x7))
#define ct_is_not_zero(val) ((val) != 0)

static inline void convert_bytes_to_GF(uint8_t *gf_array, const uint8_t *byte_array, size_t num) {
	for (size_t idx = 0; idx < num / 2; idx++) {
		gf_array[2 * idx] = byte_array[idx] & 0xf;
		gf_array[2 * idx + 1] = byte_array[idx] >> 4;
	}
	if (num & 1) {
		gf_array[num - 1] = byte_array[num / 2] & 0xf;
	}
}

// Used to compress SIG(sign)
static inline void compress_gf(uint8_t *byte_array, const uint8_t *gf_array, size_t num) {
	for (size_t idx = 0; idx < num / 2; idx++) {
		byte_array[idx] = gf_array[2 * idx] ^ (gf_array[2 * idx + 1] << 4);
	}
	if (num & 1) {
		byte_array[num / 2] = gf_array[num - 1];
	}
}

static inline void hash_combined(uint8_t *hash_out, const uint8_t *m, size_t mlen, const uint8_t *pk_seed,
                                 const uint8_t *salt) {
	shake_t state;
	snova_shake256_init(&state);
#if HASH_PK
	snova_shake_absorb(&state, pk_seed, BYTES_PK_HASH);
#else
	snova_shake_absorb(&state, pk_seed, SEED_LENGTH_PUBLIC);
#endif
	snova_shake_absorb(&state, m, mlen);
	snova_shake_absorb(&state, salt, BYTES_SALT);
	snova_shake_finalize(&state);
	snova_shake_squeeze(&state, hash_out, BYTES_HASH);
}

/**
 * Expand T12 matrix and coefficients. Shared by genkey and sign
 */
static inline void expand_T12(uint8_t *T12, const uint8_t *seed) {
	uint8_t T12coef[SNOVA_o * SNOVA_v * SNOVA_l];
	uint8_t sk_data[(SNOVA_o * SNOVA_v * SNOVA_l + 1) / 2];
	uint16_t sum;

	shake256(sk_data, (SNOVA_o * SNOVA_v * SNOVA_l + 1) / 2, seed, SEED_LENGTH_PRIVATE);
	convert_bytes_to_GF(T12coef, sk_data, SNOVA_o * SNOVA_v * SNOVA_l);

	for (size_t i1 = 0; i1 < SNOVA_o * SNOVA_v; i1++) {
		for (int j1 = 0; j1 < SNOVA_l2; j1++) {
			sum = 0;
			for (int k1 = 0; k1 < SNOVA_l; k1++) {
				sum ^= T12coef[i1 * SNOVA_l + k1] * gf_Sx[k1 * SNOVA_l2 + j1];
			}
			T12[i1 * SNOVA_l2 + j1] = gf16_compress(sum);
		}
	}

	SNOVA_SECURE_CLEAR(sk_data);
	SNOVA_SECURE_CLEAR(T12coef);
	SNOVA_SECURE_CLEARI(sum);
}

/**
 * Optimized version of genkey.
 */
int SNOVA_NAMESPACE(genkeys)(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
	/**
	 * Gen T12 matrix
	 */
	uint8_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];
	expand_T12(T12, seed + SEED_LENGTH_PUBLIC);

	/**
	 * Gen P22
	 */
	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, seed, SEED_LENGTH_PUBLIC);

#if !PUB_REORDER
	snova_pk_expander_t p12_instance;
	snova_pk_expander_init(&p12_instance, seed, SEED_LENGTH_PUBLIC);
	snova_pk_expander_goto(&p12_instance, SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2);

	snova_pk_expander_t p21_instance;
	snova_pk_expander_init(&p21_instance, seed, SEED_LENGTH_PUBLIC);
	snova_pk_expander_goto(&p21_instance, SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2);
#endif

	int idx = 0;
	uint8_t *pp22 = pk + SEED_LENGTH_PUBLIC;
	uint16_t F11i[SNOVA_o * SNOVA_l2];

	for (int mi = 0; mi < SNOVA_m1; mi++) {
		uint16_t P22x[SNOVA_o * SNOVA_o * SNOVA_l2] = {0};

		for (int nj = 0; nj < SNOVA_v; nj++) {
			memset(F11i, 0, sizeof(F11i));
			for (int nk = 0; nk < SNOVA_v; nk++) {
				uint16_t P11x[SNOVA_l2];
				snova_pk_expander_squeeze(P11x, SNOVA_l2, &instance);
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					P11x[i1] = gf16_expand(P11x[i1]);
				}
				for (int ni = 0; ni < SNOVA_o; ni++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								F11i[ni * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    T12[(nk * SNOVA_o + ni) * SNOVA_l2 + k1 * SNOVA_l + j1] * P11x[i1 * SNOVA_l + k1];
			}

			for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
				F11i[i1] = gf16_expand(gf16_compress(F11i[i1]));
			}

			for (int ni = 0; ni < SNOVA_o; ni++)
				for (int nk = 0; nk < SNOVA_o; nk++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								P22x[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    T12[(nj * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    F11i[ni * SNOVA_l2 + k1 * SNOVA_l + j1];
		}

		/**
		 * Calculate F12 and F21 contributions
		 */
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int ni = 0; ni < SNOVA_o; ni++) {
				uint16_t P12x[SNOVA_l2];
#if !PUB_REORDER
				snova_pk_expander_squeeze(P12x, SNOVA_l2, &p12_instance);
#else
				snova_pk_expander_squeeze(P12x, SNOVA_l2, &instance);
#endif
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					P12x[i1] = gf16_expand(P12x[i1]);
				}
				for (int nk = 0; nk < SNOVA_o; nk++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								P22x[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    T12[(nj * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    P12x[k1 * SNOVA_l + j1];
			}

		for (int nk = 0; nk < SNOVA_o; nk++)
			for (int nj = 0; nj < SNOVA_v; nj++) {
				uint16_t P21x[SNOVA_l2];
#if !PUB_REORDER
				snova_pk_expander_squeeze(P21x, SNOVA_l2, &p21_instance);
#else
				snova_pk_expander_squeeze(P21x, SNOVA_l2, &instance);
#endif
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					P21x[i1] = gf16_expand(P21x[i1]);
				}
				for (int ni = 0; ni < SNOVA_o; ni++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								P22x[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    P21x[i1 * SNOVA_l + k1] * T12[(nj * SNOVA_o + ni) * SNOVA_l2 + k1 * SNOVA_l + j1];
			}

		for (int i1 = 0; i1 < SNOVA_o * SNOVA_o * SNOVA_l2; i1++) {
			if (idx & 1) {
				pp22[idx / 2] ^= (gf16_compress(P22x[i1]) << 4) & 0xf0;
				idx++;
			} else {
				pp22[idx / 2] = gf16_compress(P22x[i1]) & 0xf;
				idx++;
			}
		}
	}

	snova_pk_expander_free(&instance);
#if !PUB_REORDER
	snova_pk_expander_free(&p12_instance);
	snova_pk_expander_free(&p21_instance);
#endif

	/**
	 * Output public and secret keys
	 */
	memcpy(pk, seed, SEED_LENGTH_PUBLIC);
	memcpy(sk, seed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);

#if HASH_PK
	shake256(sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, BYTES_PK_HASH, pk, BYTES_PK);
#endif

	SNOVA_SECURE_CLEAR(T12);
	SNOVA_SECURE_CLEAR(F11i);

	return 0;
}

#if NIBBLE_UINT8_T

#define get_whipped_sig(x) (((x) & 1) ? whipped_sig[(x) / 2] >> 4 : whipped_sig[(x) / 2] & 0xf)

#define set_gauss16(x, y, val)                                          \
    if ((y) & 1)                                                        \
        gauss16[x][(y) / 2] = (val << 4) ^ (gauss16[x][(y) / 2] & 0xf); \
    else                                                                \
        gauss16[x][(y) / 2] = val ^ (gauss16[x][(y) / 2] & 0xf0);
#define addto_gauss16(x, y, val)         \
    if ((y) & 1)                         \
        gauss16[x][(y) / 2] ^= val << 4; \
    else                                 \
        gauss16[x][(y) / 2] ^= val;
#define get_gauss16(x, y) (((y) & 1) ? gauss16[x][(y) / 2] >> 4 : gauss16[x][(y) / 2] & 0xf)

#else

#define get_whipped_sig(x) (whipped_sig[x])

#define set_gauss16(x, y, val) gauss16[x][y] = val
#define addto_gauss16(x, y, val) gauss16[x][y] ^= val
#define get_gauss16(x, y) (gauss16[x][y])

#endif

/**
 * Optimized version of Sign. Deterministic using the salt provided
 */
int SNOVA_NAMESPACE(sign)(const expanded_SK *skx_arg, uint8_t *sig, const uint8_t *digest, const size_t len_digest,
                          const uint8_t *salt) {
	uint8_t *sk = (uint8_t *)skx_arg;
	uint8_t *s_pk_seed = sk;
	uint8_t *s_sk_seed = sk + SEED_LENGTH_PUBLIC;
#if HASH_PK
	uint8_t *s_pk_hash = sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE;
#endif

	uint8_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];
	expand_T12(T12, s_sk_seed);

	// Calculate message has of size l^2o
	uint8_t hash_in_GF16[GF16_HASH];

	{
		uint8_t sign_hashb[BYTES_HASH];
#if HASH_PK
		hash_combined(sign_hashb, digest, len_digest, s_pk_hash, salt);
#else
		hash_combined(sign_hashb, digest, len_digest, s_pk_seed, salt);
#endif
		convert_bytes_to_GF(hash_in_GF16, sign_hashb, GF16_HASH);
	}

	// Find a solution for T.X
	uint8_t signature_in_GF[SNOVA_n * SNOVA_lr] = {0};
	int flag_redo;
	uint8_t num_sign = 0;

	while (1) {
		num_sign++;
		if (num_sign == 255) {
			// Probability of getting here is about q^{-255}
			memset(sig, 0, BYTES_SIGNATURE);
			return -1;
		}
		flag_redo = 0;

#if NIBBLE_UINT8_T
		uint8_t gauss16[SNOVA_o * SNOVA_lr][SNOVA_olr / 2];
#else
		uint8_t gauss16[SNOVA_o * SNOVA_lr][SNOVA_olr];
#endif

		// generate the vinegar value
#if NIBBLE_UINT8_T
		uint8_t whipped_sig[(SNOVA_l * SNOVA_v * SNOVA_lr + 1) / 2] = {0};
#else
		uint8_t whipped_sig[SNOVA_l * SNOVA_v * SNOVA_lr] = {0};
#endif

		{
			uint8_t vinegar_in_byte[NUM_GEN_SEC_BYTES];
			shake_t v_instance;

			snova_shake256_init(&v_instance);
			snova_shake_absorb(&v_instance, s_sk_seed, SEED_LENGTH_PRIVATE);
			snova_shake_absorb(&v_instance, digest, BYTES_DIGEST);
			snova_shake_absorb(&v_instance, salt, BYTES_SALT);
			snova_shake_absorb(&v_instance, &num_sign, 1);
			snova_shake_finalize(&v_instance);
			snova_shake_squeeze(&v_instance, vinegar_in_byte, NUM_GEN_SEC_BYTES);

			convert_bytes_to_GF(signature_in_GF, vinegar_in_byte, SNOVA_v * SNOVA_lr);

			SNOVA_SECURE_CLEARI(v_instance);
			SNOVA_SECURE_CLEAR(vinegar_in_byte);
		}

		{
			/**
			 * Whip signature
			 */
			uint16_t whipped_sig1[SNOVA_lr] = {0};

			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int i1 = 0; i1 < SNOVA_l; i1++) {
					memset(whipped_sig1, 0, sizeof(whipped_sig1));
					for (int ab = 0; ab < SNOVA_l; ++ab)
						for (int j1 = 0; j1 < SNOVA_r; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								whipped_sig1[ab * SNOVA_r + j1] ^=
								    gf_Sx[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r + j1];

#if NIBBLE_UINT8_T
					for (int j1 = 0; j1 < SNOVA_lr; j1++) {
						int idx = (i1 * SNOVA_v + ni) * SNOVA_lr + j1;
						if (idx & 1) {
							whipped_sig[idx / 2] ^= gf16_compress(whipped_sig1[j1]) << 4;
						} else {
							whipped_sig[idx / 2] ^= gf16_compress(whipped_sig1[j1]);
						}
					}
#else
					for (int j1 = 0; j1 < SNOVA_lr; j1++) {
						whipped_sig[(i1 * SNOVA_v + ni) * SNOVA_lr + j1] = gf16_compress(whipped_sig1[j1]);
					}
#endif
				}

			SNOVA_SECURE_CLEAR(whipped_sig1);
		}

		// Calculate Fvv
		memset(gauss16, 0, sizeof(gauss16));

		// Set the last column of gauss matrix
		for (int mi = 0; mi < SNOVA_o; ++mi)
			for (int i1 = 0; i1 < SNOVA_lr; i1++) {
				addto_gauss16(mi * SNOVA_lr + i1, SNOVA_o * SNOVA_lr, hash_in_GF16[mi * SNOVA_lr + i1]);
			}

		{
			/**
			 * Public matrix
			 */
			snova_pk_expander_t instance;
			snova_pk_expander_init(&instance, s_pk_seed, SEED_LENGTH_PUBLIC);

#if !PUB_REORDER
			snova_pk_expander_t p12_instance;
			snova_pk_expander_init(&p12_instance, s_pk_seed, SEED_LENGTH_PUBLIC);
			snova_pk_expander_goto(&p12_instance, SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2);

			snova_pk_expander_t p21_instance;
			snova_pk_expander_init(&p21_instance, s_pk_seed, SEED_LENGTH_PUBLIC);
			snova_pk_expander_goto(&p21_instance, SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2);
#endif

			/**
			 * ****************************************************************************************************************************************************
			 */

			uint16_t sum_t1[SNOVA_l * SNOVA_r * SNOVA_lr];
			uint8_t whipped_F12[SNOVA_l * SNOVA_o * SNOVA_lr];
			uint8_t whipped_F21[SNOVA_l * SNOVA_o * SNOVA_lr];
			uint16_t sum;

			for (int mi_prime = 0; mi_prime < SNOVA_m1; ++mi_prime) {
				memset(sum_t1, 0, sizeof(sum_t1));
				memset(whipped_F12, 0, sizeof(whipped_F12));
				memset(whipped_F21, 0, sizeof(whipped_F21));

				{
					uint16_t sum_t0[SNOVA_l * SNOVA_lr];
					uint16_t F12[SNOVA_o * SNOVA_l2];

					for (int ni = 0; ni < SNOVA_v; ++ni) {
						// P11
						memset(sum_t0, 0, sizeof(sum_t0));
						memset(F12, 0, sizeof(F12));

						for (int nj = 0; nj < SNOVA_v; ++nj) {
							uint16_t P11x[SNOVA_l2];
							snova_pk_expander_squeeze(P11x, SNOVA_l2, &instance);
							for (int i1 = 0; i1 < SNOVA_l2; i1++) {
								P11x[i1] = gf16_expand(P11x[i1]);
							}

							for (int k1 = 0; k1 < SNOVA_l; k1++)
								for (int i1 = 0; i1 < SNOVA_l; i1++)
									for (int b1 = 0; b1 < SNOVA_lr; ++b1)
										sum_t0[i1 * SNOVA_lr + b1] ^=
										    P11x[i1 * SNOVA_l + k1] *
										    get_whipped_sig((k1 * SNOVA_v + nj) * SNOVA_lr + b1);

							for (int nk = 0; nk < SNOVA_o; nk++)
								for (int i1 = 0; i1 < SNOVA_l; i1++)
									for (int j1 = 0; j1 < SNOVA_l; j1++)
										for (int k1 = 0; k1 < SNOVA_l; k1++)
											F12[nk * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
											    T12[(nj * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1] *
											    P11x[i1 * SNOVA_l + k1];
						}

						for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; ++i1) {
							sum_t0[i1] = gf16_expand(gf16_compress(sum_t0[i1]));
						}

						// Left, transposed whipped_sig
						for (int a1 = 0; a1 < SNOVA_l; ++a1)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								for (int i1 = 0; i1 < SNOVA_r; i1++)
									for (int b1 = 0; b1 < SNOVA_lr; ++b1)
										sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1] ^=
										    get_whipped_sig((k1 * SNOVA_v + ni) * SNOVA_lr + a1 * SNOVA_r + i1) *
										    sum_t0[k1 * SNOVA_lr + b1];

						/**
						 * Start F12
						 */
						for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
							F12[i1] = gf16_expand(gf16_compress(F12[i1]));
						}

						for (int nk = 0; nk < SNOVA_o; nk++)
							for (int i1 = 0; i1 < SNOVA_l; i1++)
								for (int b1 = 0; b1 < SNOVA_lr; ++b1) {
									sum = 0;
									for (int k1 = 0; k1 < SNOVA_l; k1++) {
										sum ^= F12[nk * SNOVA_l2 + k1 * SNOVA_l + i1] *
										       get_whipped_sig((k1 * SNOVA_v + ni) * SNOVA_lr + b1);
									}
									whipped_F12[(i1 * SNOVA_o + nk) * SNOVA_lr + b1] ^= gf16_compress(sum);
								}

						/**
						 * Start whipped_F21
						 */
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int b1 = 0; b1 < SNOVA_l; ++b1)
								for (int nk = 0; nk < SNOVA_o; nk++)
									for (int i1 = 0; i1 < SNOVA_l; i1++)
										for (int j1 = 0; j1 < SNOVA_r; j1++)
											whipped_F21[(b1 * SNOVA_o + nk) * SNOVA_lr + i1 * SNOVA_r + j1] ^=
											    gf16_compress(T12[(ni * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + k1] * sum_t0[k1 * SNOVA_lr + b1 * SNOVA_r + j1]);
					}

					SNOVA_SECURE_CLEAR(sum_t0);
					SNOVA_SECURE_CLEAR(F12);
				}

				for (int i1 = 0; i1 < SNOVA_l * SNOVA_r * SNOVA_lr; i1++) {
					sum_t1[i1] = gf16_compress(sum_t1[i1]);
				}

				// Right
				for (int nj = 0; nj < SNOVA_v; ++nj)
					for (int idx = 0; idx < SNOVA_o; idx++) {
						uint16_t P12x[SNOVA_l2];
#if !PUB_REORDER
						snova_pk_expander_squeeze(P12x, SNOVA_l2, &p12_instance);
#else
						snova_pk_expander_squeeze(P12x, SNOVA_l2, &instance);
#endif
						for (int i1 = 0; i1 < SNOVA_l2; i1++) {
							P12x[i1] = gf16_expand(P12x[i1]);
						}

						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1) {
								sum = 0;
								for (int k1 = 0; k1 < SNOVA_l; k1++) {
									sum ^= P12x[k1 * SNOVA_l + i1] * get_whipped_sig((k1 * SNOVA_v + nj) * SNOVA_lr + b1);
								}
								whipped_F12[(i1 * SNOVA_o + idx) * SNOVA_lr + b1] ^= gf16_compress(sum);
							}
					}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int nj = 0; nj < SNOVA_v; ++nj) {
						// Right
						uint16_t P21x[SNOVA_l2];
#if !PUB_REORDER
						snova_pk_expander_squeeze(P21x, SNOVA_l2, &p21_instance);
#else
						snova_pk_expander_squeeze(P21x, SNOVA_l2, &instance);
#endif
						for (int i1 = 0; i1 < SNOVA_l2; i1++) {
							P21x[i1] = gf16_expand(P21x[i1]);
						}

						for (int b1 = 0; b1 < SNOVA_l; ++b1)
							for (int i1 = 0; i1 < SNOVA_l; i1++)
								for (int j1 = 0; j1 < SNOVA_r; j1++) {
									sum = 0;
									for (int k1 = 0; k1 < SNOVA_l; k1++) {
										sum ^= P21x[i1 * SNOVA_l + k1] *
										       get_whipped_sig((k1 * SNOVA_v + nj) * SNOVA_lr + b1 * SNOVA_r + j1);
									}
									whipped_F21[(b1 * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1] ^= gf16_compress(sum);
								}
					}

				/**
				 * Alpha loop **********************************************************************************************************
				 *
				 * Apply A, B, q1 and q2, aka E matrix
				 */

				uint16_t gfm_temp0[SNOVA_olr];
				uint16_t gfm_temp1[SNOVA_olr];
				uint16_t gfm_temp2[SNOVA_olr];
				uint16_t Q12[SNOVA_l2];

				for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
					int mi = (mi_prime + SNOVA_m1 * SNOVA_alpha - alpha) % SNOVA_m1;

					if (mi >= SNOVA_o) {
						continue;
					}

					{
						/**
						 * F_vv
						 */
						memset(gfm_temp0, 0, sizeof(gfm_temp0));
						memset(gfm_temp1, 0, sizeof(gfm_temp1));
						memset(gfm_temp2, 0, sizeof(gfm_temp2));

						// apply q1 and q2
						for (int a1 = 0; a1 < SNOVA_l; ++a1) {
							memset(gfm_temp0, 0, sizeof(gfm_temp0));

							for (int b1 = 0; b1 < SNOVA_l; ++b1)
								for (int i1 = 0; i1 < SNOVA_r; i1++)
									for (int j1 = 0; j1 < SNOVA_r; j1++)
										gfm_temp0[i1 * SNOVA_r + j1] ^=
										    sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1 * SNOVA_r + j1] *
										    fixedq2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

							for (int i1 = 0; i1 < SNOVA_r2; i1++) {
								gfm_temp0[i1] = gf16_compress(gfm_temp0[i1]);
							}

							for (int i1 = 0; i1 < SNOVA_r2; i1++) {
								gfm_temp1[i1] ^= gfm_temp0[i1] * fixedq1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
							}
						}

						for (int i1 = 0; i1 < SNOVA_r2; i1++) {
							gfm_temp1[i1] = gf16_compress(gfm_temp1[i1]);
						}

						// A and B
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								for (int k1 = 0; k1 < SNOVA_r; k1++) {
									gfm_temp2[i1 * SNOVA_l + j1] ^=
									    gfm_temp1[i1 * SNOVA_r + k1] * fixedBm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];
								}

						for (int i1 = 0; i1 < SNOVA_lr; i1++) {
							gfm_temp2[i1] = gf16_compress(gfm_temp2[i1]);
						}

						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								for (int k1 = 0; k1 < SNOVA_r; k1++) {
									addto_gauss16(mi * SNOVA_lr + i1 * SNOVA_l + j1, SNOVA_o * SNOVA_lr,
									              gf16_compress(fixedAm[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] * gfm_temp2[k1 * SNOVA_l + j1]));
								}
					}
				}

				// compute the coefficients of Xo and put into gauss matrix and compute
				// the coefficients of Xo^t and add into gauss matrix

				for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
					int mi = (mi_prime + SNOVA_m1 * SNOVA_alpha - alpha) % SNOVA_m1;

					if (mi >= SNOVA_o) {
						continue;
					}

					{
						memset(gfm_temp0, 0, sizeof(gfm_temp0));
						memset(gfm_temp1, 0, sizeof(gfm_temp1));
						memset(gfm_temp2, 0, sizeof(gfm_temp2));

						// Transpose
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int idx = 0; idx < SNOVA_o; idx++)
								for (int b1 = 0; b1 < SNOVA_l; ++b1)
									for (int j1 = 0; j1 < SNOVA_r; j1++)
										gfm_temp0[idx * SNOVA_lr + i1 * SNOVA_r + j1] ^=
										    whipped_F12[(i1 * SNOVA_o + idx) * SNOVA_lr + b1 * SNOVA_r + j1] *
										    fixedq1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

						for (int i1 = 0; i1 < SNOVA_olr; i1++) {
							gfm_temp0[i1] = gf16_compress(gfm_temp0[i1]);
						}

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int i1 = 0; i1 < SNOVA_r; i1++)
								for (int j1 = 0; j1 < SNOVA_l; j1++)
									for (int k1 = 0; k1 < SNOVA_r; k1++)
										gfm_temp1[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
										    fixedAm[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] *
										    gfm_temp0[idx * SNOVA_lr + j1 * SNOVA_r + k1];

						for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
							gfm_temp1[i1] = gf16_compress(gfm_temp1[i1]);
						}

						for (int i1 = 0; i1 < SNOVA_l2; i1++) {
							sum = 0;
							for (int j1 = 0; j1 < SNOVA_l; j1++) {
								sum ^= fixedq2[(mi * SNOVA_alpha + alpha) * SNOVA_l + j1] * gf_S[j1 * SNOVA_l2 + i1];
							}
							Q12[i1] = gf16_expand(gf16_compress(sum));
						}

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int i1 = 0; i1 < SNOVA_r; i1++)
								for (int j1 = 0; j1 < SNOVA_l; j1++)
									for (int k1 = 0; k1 < SNOVA_l; k1++) {
										gfm_temp2[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^= gfm_temp1[idx * SNOVA_lr + i1 * SNOVA_l + k1] * Q12[k1 * SNOVA_l + j1];
									}

						for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
							gfm_temp2[i1] = gf16_compress(gfm_temp2[i1]);
						}

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
								for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
									for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
										for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
											addto_gauss16(mi * SNOVA_lr + ti1 * SNOVA_l + ti2, idx * SNOVA_lr + tj1 * SNOVA_r + tj2,
											              gf16_compress(gfm_temp2[idx * SNOVA_lr + ti1 * SNOVA_l + tj1] * fixedBm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2]));
					}

					/**
					 * End of F12
					 * *****************************************
					 * Start F21
					 */

					{
						memset(gfm_temp0, 0, sizeof(gfm_temp0));
						memset(gfm_temp1, 0, sizeof(gfm_temp1));
						memset(gfm_temp2, 0, sizeof(gfm_temp2));

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int b1 = 0; b1 < SNOVA_l; ++b1)
								for (int i1 = 0; i1 < SNOVA_l; i1++)
									for (int j1 = 0; j1 < SNOVA_r; j1++)
										gfm_temp0[idx * SNOVA_lr + i1 * SNOVA_r + j1] ^=
										    whipped_F21[(b1 * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1] *
										    fixedq2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

						for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
							gfm_temp0[i1] = gf16_compress(gfm_temp0[i1]);
						}

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int i1 = 0; i1 < SNOVA_l; i1++)
								for (int j1 = 0; j1 < SNOVA_l; j1++)
									for (int k1 = 0; k1 < SNOVA_r; k1++) {
										gfm_temp1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
										    gfm_temp0[idx * SNOVA_lr + i1 * SNOVA_r + k1] *
										    fixedBm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];
									}

						for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
							gfm_temp1[i1] = gf16_compress(gfm_temp1[i1]);
						}

						for (int i1 = 0; i1 < SNOVA_l2; i1++) {
							sum = 0;
							for (int j1 = 0; j1 < SNOVA_l; j1++) {
								sum ^= fixedq1[(mi * SNOVA_alpha + alpha) * SNOVA_l + j1] * gf_S[j1 * SNOVA_l2 + i1];
							}
							Q12[i1] = gf16_expand(gf16_compress(sum));
						}

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int i1 = 0; i1 < SNOVA_l; i1++)
								for (int j1 = 0; j1 < SNOVA_l; j1++)
									for (int k1 = 0; k1 < SNOVA_l; k1++) {
										gfm_temp2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^= Q12[i1 * SNOVA_l + k1] * gfm_temp1[idx * SNOVA_l2 + k1 * SNOVA_l + j1];
									}

						for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
							gfm_temp2[i1] = gf16_compress(gfm_temp2[i1]);
						}

						for (int idx = 0; idx < SNOVA_o; idx++)
							for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
								for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
									for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
										for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
											addto_gauss16(mi * SNOVA_lr + ti1 * SNOVA_l + ti2, idx * SNOVA_lr + tj1 * SNOVA_r + tj2,
											              gf16_compress(gfm_temp2[idx * SNOVA_l2 + tj1 * SNOVA_l + ti2] * fixedAm[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r + tj2]));
					}
				}

				SNOVA_SECURE_CLEAR(gfm_temp0);
				SNOVA_SECURE_CLEAR(gfm_temp1);
				SNOVA_SECURE_CLEAR(gfm_temp2);
				SNOVA_SECURE_CLEAR(Q12);
			}

			SNOVA_SECURE_CLEAR(sum_t1);
			SNOVA_SECURE_CLEAR(whipped_F12);
			SNOVA_SECURE_CLEAR(whipped_F21);
			SNOVA_SECURE_CLEARI(sum);

			snova_pk_expander_free(&instance);
#if !PUB_REORDER
			snova_pk_expander_free(&p12_instance);
			snova_pk_expander_free(&p21_instance);
#endif
		}

		/**
		 * Gaussian elimination
		 */

		// Gaussian elimination in constant time
		uint16_t t_GF16;
		uint16_t fact;

		for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) {
			for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
				int16_t mask = ct_is_not_zero(get_gauss16(i, i)) - 1;
				for (int k = 0; k < SNOVA_olr; ++k) {
					addto_gauss16(i, k, mask & get_gauss16(j, k));
				}
			}

			flag_redo |= 1 - ct_is_not_zero(get_gauss16(i, i));

			/**
			 * Constant time GF(16) inverse
			 *
			 * Use that x^q = x and therefore x^(q-2) = x^-1
			 *
			 * t_GF16 = ct_gf_inverse(get_gauss16(i, i));
			 */
			fact = get_gauss16(i, i) * gf16_expand(get_gauss16(i, i));
			fact = gf16_compress(fact);
			t_GF16 = fact;

			fact = fact * gf16_expand(fact);
			fact = gf16_compress(fact);
			t_GF16 = t_GF16 * gf16_expand(fact);
			t_GF16 = gf16_compress(t_GF16);

			fact = fact * gf16_expand(fact);
			fact = gf16_compress(fact);
			t_GF16 = t_GF16 * gf16_expand(fact);
			t_GF16 = gf16_compress(t_GF16);
			t_GF16 = gf16_expand(t_GF16);

			for (int k = 0; k < SNOVA_olr; ++k) {
				set_gauss16(i, k, gf16_compress(get_gauss16(i, k) * t_GF16));
			}

			for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
				uint16_t gji = gf16_expand(get_gauss16(j, i));
				for (int k = 0; k < SNOVA_olr; ++k) {
					addto_gauss16(j, k, gf16_compress(get_gauss16(i, k) * gji));
				}
			}
		}

		SNOVA_SECURE_CLEARI(t_GF16);
		SNOVA_SECURE_CLEARI(fact);

		if (flag_redo) {
			continue;
		}

		// Last step of Gaussian elimination
		uint8_t solution[SNOVA_olr] = {0};
		uint16_t solution16[SNOVA_o * SNOVA_lr] = {0};
		uint16_t sum;

		for (int i = SNOVA_o * SNOVA_lr - 1; i >= 0; --i) {
			sum = 0;
			for (int k = i + 1; k < SNOVA_o * SNOVA_lr; ++k) {
				sum ^= get_gauss16(i, k) * solution16[k];
			}
			solution16[i] = gf16_expand(gf16_compress(sum ^ get_gauss16(i, SNOVA_o * SNOVA_lr)));
		}
		for (int i = 0; i < SNOVA_o * SNOVA_lr; i++) {
			solution[i] = gf16_compress(solution16[i]);
		}

		memcpy(signature_in_GF + SNOVA_v * SNOVA_lr, solution, SNOVA_o * SNOVA_lr);

		// Establish signature using T12
		for (int index = 0; index < SNOVA_v; ++index)
			for (int mi = 0; mi < SNOVA_o; ++mi)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_r; j1++) {
						sum = 0;
						for (int k1 = 0; k1 < SNOVA_l; k1++) {
							sum ^= gf16_expand(T12[(index * SNOVA_o + mi) * SNOVA_l2 + i1 * SNOVA_l + k1]) *
							       solution[mi * SNOVA_lr + k1 * SNOVA_r + j1];
						}
						signature_in_GF[index * SNOVA_lr + i1 * SNOVA_r + j1] ^= gf16_compress(sum);
					}

		memcpy(signature_in_GF + SNOVA_v * SNOVA_lr, solution, SNOVA_o * SNOVA_lr);

		SNOVA_SECURE_CLEAR(gauss16);
		SNOVA_SECURE_CLEAR(whipped_sig);
		SNOVA_SECURE_CLEAR(solution);
		SNOVA_SECURE_CLEAR(solution16);
		SNOVA_SECURE_CLEARI(sum);

		break;
	}

	SNOVA_SECURE_CLEAR(T12);

	compress_gf(sig, signature_in_GF, SNOVA_n * SNOVA_lr);
	memcpy(sig + BYTES_SIGNATURE - BYTES_SALT, salt, BYTES_SALT);

	return 0;
}

/**
 * Optimized version of verify.
 */
int SNOVA_NAMESPACE(verify)(const expanded_PK *pkx, const uint8_t *sig, const uint8_t *digest, const size_t len_digest) {
	uint8_t *pk = (uint8_t *)pkx;
	uint8_t *pk_seed = pk;

	// Prepare
	uint16_t hash_in_GF[SNOVA_o * SNOVA_lr] = {0};
	{
		uint8_t *pp22 = pk + SEED_LENGTH_PUBLIC;
		int idx_p22 = 0;

		/**
		 * Whip signature
		 */
#if NIBBLE_UINT8_T
		uint8_t whipped_sig[(SNOVA_l * SNOVA_n * SNOVA_lr + 1) / 2] = {0};
#else
		uint8_t whipped_sig[SNOVA_l * SNOVA_n * SNOVA_lr] = {0};
#endif
		for (int idx = 0; idx < SNOVA_n; ++idx) {
			uint16_t whipped_sig1[SNOVA_l * SNOVA_lr] = {0};

			for (int ab = 0; ab < SNOVA_l; ++ab) {
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_r; j1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++) {
							if ((idx * SNOVA_lr + k1 * SNOVA_r + j1) & 1) {
								whipped_sig1[i1 * SNOVA_lr + ab * SNOVA_r + j1] ^=
								    gf_Sx[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * (sig[(idx * SNOVA_lr + k1 * SNOVA_r + j1) / 2] >> 4);
							} else {
								whipped_sig1[i1 * SNOVA_lr + ab * SNOVA_r + j1] ^=
								    gf_Sx[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * (sig[(idx * SNOVA_lr + k1 * SNOVA_r + j1) / 2] & 0xf);
							}
						}
			}

#if NIBBLE_UINT8_T
			for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; i1++) {
				if ((idx * SNOVA_l * SNOVA_lr + i1) & 1) {
					whipped_sig[(idx * SNOVA_l * SNOVA_lr + i1) / 2] ^= gf16_compress(whipped_sig1[i1]) << 4;
				} else {
					whipped_sig[(idx * SNOVA_l * SNOVA_lr + i1) / 2] ^= gf16_compress(whipped_sig1[i1]);
				}
			}
#else
			for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; i1++) {
				whipped_sig[idx * SNOVA_l * SNOVA_lr + i1] = gf16_compress(whipped_sig1[i1]);
			}
#endif
		}

		/**
		 * stream P matrix
		 */
		snova_pk_expander_t instance;
		snova_pk_expander_init(&instance, pk_seed, SEED_LENGTH_PUBLIC);

#if !PUB_REORDER
		snova_pk_expander_t p12_instance;
		snova_pk_expander_init(&p12_instance, pk_seed, SEED_LENGTH_PUBLIC);
		snova_pk_expander_goto(&p12_instance, SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2);

		snova_pk_expander_t p21_instance;
		snova_pk_expander_init(&p21_instance, pk_seed, SEED_LENGTH_PUBLIC);
		snova_pk_expander_goto(&p21_instance, SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2);
#endif

		for (int mi_prime = 0; mi_prime < SNOVA_m1; ++mi_prime) {
			uint16_t sum_t1[SNOVA_l * SNOVA_r * SNOVA_lr] = {0};

			for (int ni = 0; ni < SNOVA_o; ++ni) {
				// Right
				uint16_t sum_t0[SNOVA_l * SNOVA_lr] = {0};
				for (int nj = 0; nj < SNOVA_o; ++nj) {
					uint16_t P22[SNOVA_l2];

					for (int i1 = 0; i1 < SNOVA_l2; i1++) {
						if (idx_p22 & 1) {
							P22[i1] = gf16_expand(pp22[idx_p22 / 2] >> 4);
						} else {
							P22[i1] = gf16_expand(pp22[idx_p22 / 2] & 0xf);
						}
						idx_p22++;
					}

					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t0[i1 * SNOVA_lr + b1] ^=
								    P22[i1 * SNOVA_l + k1] *
								    get_whipped_sig((nj + SNOVA_v) * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + b1);
				}

				for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; ++i1) {
					sum_t0[i1] = gf16_expand(gf16_compress(sum_t0[i1]));
				}

				// Left, transposed whipped_sig
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1] ^=
								    get_whipped_sig((ni + SNOVA_v) * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + a1 * SNOVA_r + i1) *
								    sum_t0[k1 * SNOVA_lr + b1];
			}

			// P11
			for (int ni = 0; ni < SNOVA_v; ++ni) {
				// Right
				uint16_t sum_t0[SNOVA_l * SNOVA_lr] = {0};
				for (int nj = 0; nj < SNOVA_v; ++nj) {
					uint16_t P11x[SNOVA_l2];
					snova_pk_expander_squeeze(P11x, SNOVA_l2, &instance);
					for (int i1 = 0; i1 < SNOVA_l2; i1++) {
						P11x[i1] = gf16_expand(P11x[i1]);
					}

					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t0[i1 * SNOVA_lr + b1] ^=
								    P11x[i1 * SNOVA_l + k1] *
								    get_whipped_sig(nj * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + b1);
				}

				for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; ++i1) {
					sum_t0[i1] = gf16_expand(gf16_compress(sum_t0[i1]));
				}

				// Left, transposed whipped_sig
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1] ^=
								    get_whipped_sig(ni * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + a1 * SNOVA_r + i1) *
								    sum_t0[k1 * SNOVA_lr + b1];
			}

			// P12
			for (int ni = 0; ni < SNOVA_v; ++ni) {
				// Right
				uint16_t sum_t0[SNOVA_l * SNOVA_lr] = {0};
				for (int nj = 0; nj < SNOVA_o; ++nj) {
					uint16_t P12x[SNOVA_l2];
#if !PUB_REORDER
					snova_pk_expander_squeeze(P12x, SNOVA_l2, &p12_instance);
#else
					snova_pk_expander_squeeze(P12x, SNOVA_l2, &instance);
#endif
					for (int i1 = 0; i1 < SNOVA_l2; i1++) {
						P12x[i1] = gf16_expand(P12x[i1]);
					}
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t0[i1 * SNOVA_lr + b1] ^=
								    P12x[i1 * SNOVA_l + k1] *
								    get_whipped_sig((nj + SNOVA_v) * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + b1);
				}

				for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; ++i1) {
					sum_t0[i1] = gf16_expand(gf16_compress(sum_t0[i1]));
				}

				// Left, transposed whipped_sig
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1] ^=
								    get_whipped_sig(ni * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + a1 * SNOVA_r + i1) *
								    sum_t0[k1 * SNOVA_lr + b1];
			}

			// P21
			for (int ni = 0; ni < SNOVA_o; ++ni) {
				// Right
				uint16_t sum_t0[SNOVA_l * SNOVA_lr] = {0};
				for (int nj = 0; nj < SNOVA_v; ++nj) {
					uint16_t P21x[SNOVA_l2];
#if !PUB_REORDER
					snova_pk_expander_squeeze(P21x, SNOVA_l2, &p21_instance);
#else
					snova_pk_expander_squeeze(P21x, SNOVA_l2, &instance);
#endif
					for (int i1 = 0; i1 < SNOVA_l2; i1++) {
						P21x[i1] = gf16_expand(P21x[i1]);
					}

					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t0[i1 * SNOVA_lr + b1] ^=
								    P21x[i1 * SNOVA_l + k1] *
								    get_whipped_sig(nj * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + b1);
				}

				for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr; ++i1) {
					sum_t0[i1] = gf16_expand(gf16_compress(sum_t0[i1]));
				}

				// Left, transposed whipped_sig
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr; ++b1)
								sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1] ^=
								    get_whipped_sig((ni + SNOVA_v) * SNOVA_l * SNOVA_lr + k1 * SNOVA_lr + a1 * SNOVA_r + i1) *
								    sum_t0[k1 * SNOVA_lr + b1];
			}

			for (int i1 = 0; i1 < SNOVA_l * SNOVA_r * SNOVA_lr; i1++) {
				sum_t1[i1] = gf16_compress(sum_t1[i1]);
			}

			/**
			 * Apply A, B, q1 and q2, aka E matrix
			 */
			for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
				int mi = (mi_prime + SNOVA_m1 * SNOVA_alpha - alpha) % SNOVA_m1;

				if (mi >= SNOVA_o) {
					continue;
				}

				uint16_t gfm_temp1[SNOVA_r2] = {0};
				uint16_t gfm_temp2[SNOVA_lr] = {0};

				// apply q1 and q2
				for (int a1 = 0; a1 < SNOVA_l; ++a1) {
					uint16_t gfm_temp0[SNOVA_r2] = {0};

					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								gfm_temp0[i1 * SNOVA_r + j1] ^=
								    fixedq2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1] *
								    sum_t1[a1 * SNOVA_r * SNOVA_lr + i1 * SNOVA_lr + b1 * SNOVA_r + j1];

					for (int i1 = 0; i1 < SNOVA_r2; i1++) {
						gfm_temp0[i1] = gf16_compress(gfm_temp0[i1]);
					}

					for (int i1 = 0; i1 < SNOVA_r2; i1++) {
						gfm_temp1[i1] ^= gfm_temp0[i1] * fixedq1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
					}
				}

				for (int i1 = 0; i1 < SNOVA_r2; i1++) {
					gfm_temp1[i1] = gf16_compress(gfm_temp1[i1]);
				}

				// A and B
				for (int i1 = 0; i1 < SNOVA_r; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int k1 = 0; k1 < SNOVA_r; k1++)
							gfm_temp2[i1 * SNOVA_l + j1] ^=
							    gfm_temp1[i1 * SNOVA_r + k1] * fixedBm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];

				for (int i1 = 0; i1 < SNOVA_lr; i1++) {
					gfm_temp2[i1] = gf16_compress(gfm_temp2[i1]);
				}

				for (int i1 = 0; i1 < SNOVA_r; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int k1 = 0; k1 < SNOVA_r; k1++)
							hash_in_GF[mi * SNOVA_lr + i1 * SNOVA_l + j1] ^=
							    fixedAm[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] * gfm_temp2[k1 * SNOVA_l + j1];
			}
		}

		snova_pk_expander_free(&instance);
#if !PUB_REORDER
		snova_pk_expander_free(&p12_instance);
		snova_pk_expander_free(&p21_instance);
#endif

		for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
			hash_in_GF[i1] = gf16_compress(hash_in_GF[i1]);
		}
	}

	/**
	 * Check hashes
	 */
	uint8_t signed_gf[GF16_HASH] = {0};
	const uint8_t *salt = sig + BYTES_SIGNATURE - BYTES_SALT;
	uint8_t signed_bytes[BYTES_HASH];
#if HASH_PK
	uint8_t pk_hash[BYTES_PK_HASH];
	shake256(pk_hash, BYTES_PK_HASH, pk, BYTES_PK);
	hash_combined(signed_bytes, digest, len_digest, pk_hash, salt);
#else
	hash_combined(signed_bytes, digest, len_digest, pk_seed, salt);
#endif
	convert_bytes_to_GF(signed_gf, signed_bytes, GF16_HASH);

	for (int i = 0; i < GF16_HASH; ++i) {
		if (hash_in_GF[i] != signed_gf[i]) {
			return -1;
		}
	}

	return 0;
}
