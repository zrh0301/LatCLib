// SPDX-License-Identifier: MIT

/**
 * Reference implementation
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <stdint.h>
#include <string.h>

#include "snova.h"
#include "symmetric.h"
#include "gf16.h"
#include "gf16_matrix.h"

#define _str(s) #s
#define str(s) _str(s)
#define CRYPTO_ALGNAME str(SNOVA_NAME)

typedef uint8_t gf_t;

#if SNOVA_q != 16
gf_t gf_multtab[SNOVA_q * SNOVA_q] = {0};
gf_t gf_invtab[SNOVA_q] = {0};
gf_t gf_addtab[SNOVA_q * SNOVA_q] = {0};
#endif
gf_t gf_S[SNOVA_l * SNOVA_l2] = {0};

static inline gf_t gf_mult(const gf_t a, const gf_t b) {
#if SNOVA_q == 16
	return gf16_mul(a, b);
#else
	return gf_multtab[a * SNOVA_q + b];
#endif
}

static inline gf_t gf_inv(const gf_t a) {
#if SNOVA_q == 16
	return gf16_inv(a);
#else
	return gf_invtab[a];
#endif
}

static inline gf_t gf_add(const gf_t a, const gf_t b) {
#if SNOVA_q == 16
	return gf16_add(a, b);
#else
	return gf_addtab[a * SNOVA_q + b];
#endif
}

static inline void gf_set_add(gf_t* a, const gf_t b) {
#if SNOVA_q == 16
	*a = gf16_add(*a, b);
#else
	*a = gf_addtab[*a * SNOVA_q + b];
#endif
}

static inline gf_t gf_sub(const gf_t a, const gf_t b) {
#if SNOVA_q != 16
	return gf_addtab[a * SNOVA_q + (SNOVA_q - b) % SNOVA_q];
#else
	return gf16_sub(a, b);
#endif
}

static inline void gf_mat_mul(gf_t* a, const gf_t* b, const gf_t* c) {
#if SNOVA_q == 16
	(void)gf16_matrix_mul(a, b, c, SNOVA_l, SNOVA_l, SNOVA_l);
#else
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_t sum = 0;
			for (int k1 = 0; k1 < SNOVA_l; k1++) {
				gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
			}
			a[i1 * SNOVA_l + j1] = sum;
		}
#endif
}

static inline void gf_mat_mul_add(gf_t* a, const gf_t* b, const gf_t* c) {
#if SNOVA_q == 16
	(void)gf16_matrix_mul_add(a, b, c, SNOVA_l, SNOVA_l, SNOVA_l);
#else
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_t sum = 0;
			for (int k1 = 0; k1 < SNOVA_l; k1++) {
				gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
			}
			gf_set_add(&a[i1 * SNOVA_l + j1], sum);
		}
#endif
}

static inline void gf_mat_mul_add16(gf_t* a, const uint16_t *b, const gf_t* c) {
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_r; j1++) {
			gf_t sum = 0;
			for (int k1 = 0; k1 < SNOVA_l; k1++) {
				gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_r + j1]));
			}
			gf_set_add(&a[i1 * SNOVA_r + j1], sum);
		}
}

static inline void gf_mat_mul_add_lr(gf_t* a, const gf_t* b, const gf_t* c, int ad, int bd, int cd) {
#if SNOVA_q == 16
	(void)gf16_matrix_mul_add(a, b, c, (size_t)ad, (size_t)bd, (size_t)cd);
#else
	for (int i1 = 0; i1 < ad; i1++)
		for (int j1 = 0; j1 < cd; j1++) {
			gf_t sum = 0;
			for (int k1 = 0; k1 < bd; k1++) {
				gf_set_add(&sum, gf_mult(b[i1 * bd + k1], c[k1 * cd + j1]));
			}
			gf_set_add(&a[i1 * cd + j1], sum);
		}
#endif
}

static inline gf_t gf_mat_det(gf_t* a) {
#if SNOVA_q == 16
	gf_t det = 0;
	(void)gf16_matrix_determinant(&det, a, SNOVA_l);
	return det;
#else
	gf_t det = 0;
#if SNOVA_l == 1
	det = a[0];
#elif SNOVA_l == 2
	det = gf_sub(gf_mult(a[0], a[3]), gf_mult(a[1], a[2]));
#elif SNOVA_l == 3
	det = gf_mult(a[0], gf_sub(gf_mult(a[4], a[8]), gf_mult(a[5], a[7])));
	gf_set_add(&det, gf_mult(a[1], gf_sub(gf_mult(a[5], a[6]), gf_mult(a[3], a[8]))));
	gf_set_add(&det, gf_mult(a[2], gf_sub(gf_mult(a[3], a[7]), gf_mult(a[4], a[6]))));
#elif SNOVA_l == 4
	gf_t det_l;
	gf_t det_r;
#define DET_L(x, y) det_l = gf_sub(gf_mult(a[x], a[4 + y]), gf_mult(a[y], a[4 + x]))
#define DET_R(x, y) det_r = gf_sub(gf_mult(a[8 + x], a[12 + y]), gf_mult(a[8 + y], a[12 + x]))
#define DET22(x1, y1, x2, y2) \
    DET_L(x1, y1);            \
    DET_R(x2, y2);            \
    gf_set_add(&det, gf_mult(det_l, det_r))
	DET22(0, 1, 2, 3);
	DET22(0, 2, 3, 1);
	DET22(0, 3, 1, 2);
	DET22(1, 2, 0, 3);
	DET22(1, 3, 2, 0);
	DET22(2, 3, 0, 1);
#undef DET_R
#undef DET22
#undef DET_L
#elif SNOVA_l == 5
	gf_t det_l;
	gf_t det_r;
#define DET_L(x, y) det_l = gf_sub(gf_mult(a[x], a[5 + y]), gf_mult(a[y], a[5 + x]))
#define DET_R2(x, y, z) gf_mult(gf_sub(gf_mult(a[10 + x], a[15 + y]), gf_mult(a[10 + y], a[15 + x])), a[20 + z])
#define DET_R3(x, y, z) det_r = gf_add(DET_R2(x, y, z), gf_add(DET_R2(y, z, x), DET_R2(z, x, y)))
#define DET23(x1, y1, x2, y2, z2) \
    DET_L(x1, y1);                \
    DET_R3(x2, y2, z2);           \
    gf_set_add(&det, gf_mult(det_l, det_r))
	DET23(0, 1, 2, 3, 4);
	DET23(0, 2, 3, 1, 4);
	DET23(0, 3, 1, 2, 4);
	DET23(0, 4, 1, 3, 2);
	DET23(1, 2, 0, 3, 4);
	DET23(1, 3, 2, 0, 4);
	DET23(1, 4, 2, 3, 0);
	DET23(2, 3, 0, 1, 4);
	DET23(2, 4, 0, 3, 1);
	DET23(3, 4, 2, 0, 1);
#undef DET_R2
#undef DET_R3
#undef DET23
#undef DET_L
#else
#error "Unsupported rank"
#endif
	return det;
#endif
}

static void init_gf_tables(void) {
#if SNOVA_q == 16
	/* GF(16) scalar and matrix operators are provided by the shared core. */
#else
	// GF(prime)
	for (int i1 = 0; i1 < SNOVA_q; i1++)
		for (int j1 = 0; j1 < SNOVA_q; j1++) {
			gf_multtab[i1 * SNOVA_q + j1] = (i1 * j1) % SNOVA_q;
			gf_addtab[i1 * SNOVA_q + j1] = (i1 + j1) % SNOVA_q;
		}
#endif
#if SNOVA_q != 16
	// Use that x^q = x and therefore x^(q-2) = x^-1
	for (int i1 = 0; i1 < SNOVA_q; i1++) {
		gf_t val = i1;
		for (int j1 = 3; j1 < SNOVA_q; j1++) {
			val = gf_mult(val, i1);
		}
		gf_invtab[i1] = val;
	}
#endif
}

// Set the irreducible S matrix
static void set_S(gf_t* gf_S1) {
#if SNOVA_q == 16
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_S1[i1 * SNOVA_l + j1] = 8 - (i1 + j1);
		}
#if SNOVA_l == 5
	gf_S1[SNOVA_l2 - 1] = 9;
#endif
#else
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_S1[i1 * SNOVA_l + j1] = (Q_A + i1 + j1) & Q_B;
		}
	gf_S1[SNOVA_l2 - 1] = Q_C;
#endif
}

static void gen_S_array(void) {
	memset(gf_S, 0, sizeof(gf_S));

	for (int i1 = 0; i1 < SNOVA_l; i1++) {
		gf_S[i1 * SNOVA_l + i1] = 1;
	}
#if SNOVA_l > 1
	set_S(&gf_S[1 * SNOVA_l2]);

	for (int i1 = 2; i1 < SNOVA_l; i1++) {
		gf_mat_mul(&gf_S[i1 * SNOVA_l2], &gf_S[1 * SNOVA_l2], &gf_S[(i1 - 1) * SNOVA_l2]);
	}
#endif
}

static int first_time = 1;

#if FIXED_ABQ
static void gen_fixed_ABQ(const char *abq_seed);
#endif

static void snova_init(void) {
	first_time = 0;
	init_gf_tables();
	gen_S_array();
#if FIXED_ABQ
// Special values only for parameter sets that are not recommended for Round 3.
#if !ROUND2_T12 && SNOVA_l == 2
#if SNOVA_o == 17
	gen_fixed_ABQ("SNOVA_ABQ_2");
#elif SNOVA_o == 25
	gen_fixed_ABQ("SNOVA_ABQ_3");
#elif SNOVA_o == 33
	gen_fixed_ABQ("SNOVA_ABQ_4");
#else
#error "Unsupported SNOVA_o for SNOVA_l ==2"
#endif
#else
	gen_fixed_ABQ("SNOVA_ABQ");
#endif
#endif
}

/**
 * Utilities
 */
static void convert_bytes_to_GF(gf_t* gf_array, const uint8_t *byte_array, size_t num) {
#if SNOVA_q != 16
	for (size_t idx = 0; idx < num; idx++) {
		gf_array[idx] = byte_array[idx] % SNOVA_q;
	}
#else
	for (size_t idx = 0; idx < num / 2; idx++) {
		gf_array[2 * idx] = (byte_array[idx] & 0xf) % SNOVA_q;
		gf_array[2 * idx + 1] = (byte_array[idx] >> 4) % SNOVA_q;
	}
	if (num & 1) {
		gf_array[num - 1] = (byte_array[num / 2] & 0xf) % SNOVA_q;
	}
#endif
}

// Used to compress PK (genkey) and SIG(sign)
static void compress_gf(uint8_t *byte_array, const gf_t* gf_array, size_t num) {
	size_t idx = 0;
	size_t out_idx = 0;
	size_t num_bytes = BYTES_GF(num);

	do {
		uint64_t val = 0;
		uint64_t fact = 1;

		int i1 = 0;
		while (i1 < PACK_GF && idx < num) {
			val += fact * (gf_array[idx] % SNOVA_q);
			idx++;
			i1++;
			fact *= SNOVA_q;
		}

		i1 = (i1 + 1) / 2;
		int j1 = 0;
		while (j1 < PACK_BYTES && out_idx < num_bytes) {
			byte_array[out_idx] = val & 0xff;
			out_idx++;
			val = val >> 8;
			j1++;
		}
	} while (idx < num);
}

// Used to expand PK(verify) and SIG(verify)
static int expand_gf(gf_t* gf_array, const uint8_t *byte_array, size_t num) {
	size_t num_bytes = BYTES_GF(num);
	size_t idx = 0;
	size_t out_idx = 0;
	uint64_t val;
	int res = 0;

	do {
		val = 0;

		int i1 = 0;
		while (i1 < PACK_BYTES && idx < num_bytes) {
			val = val ^ ((uint64_t)(byte_array[idx]) << (8 * i1));
			idx++;
			i1++;
		}

		int j1 = 0;
		while (j1 < PACK_GF && out_idx < num) {
			gf_array[out_idx] = val % SNOVA_q;
			val = val / SNOVA_q;
			out_idx++;
			j1++;
		}

		res |= val;
	} while (out_idx < num);

#if SNOVA_q == 16
	if (num & 1) {
		return byte_array[num / 2] & 0xF0;
	}
#endif

	return res;
}

// Used to compress PK (genkey)
static void compress_pk(uint8_t *pk, gf_t* P22) {
#if defined(SYMMETRIC)

#if SNOVA_q == 16
	// Convert to Upper Triangular
	for (int mi = 0; mi < SNOVA_m1; ++mi)
		for (int ni = 0; ni < SNOVA_o; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				for (int j1 = i1 + 1; j1 < SNOVA_l; j1++) {
					P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] =
					    P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^
					    P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + j1 * SNOVA_l + i1];
				}

				for (int nj = ni + 1; nj < SNOVA_o; ++nj)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P22[((mi * SNOVA_o + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    P22[((mi * SNOVA_o + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] ^
						    P22[((mi * SNOVA_o + nj) * SNOVA_o + ni) * SNOVA_l2 + j1 * SNOVA_l + i1];
					}
			}
		}
#endif
	gf_t P22c[NUMGF_PK] = {0};
	gf_t* curval = &P22c[0];

	for (int mi = 0; mi < SNOVA_m1; ++mi)
		for (int ni = 0; ni < SNOVA_o; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
					*curval = P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1];
					curval++;
				}

				for (int nj = ni + 1; nj < SNOVA_o; ++nj)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						*curval = P22[((mi * SNOVA_o + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
						curval++;
					}
			}
		}

	compress_gf(pk, P22c, NUMGF_PK);
#else

	compress_gf(pk, P22, NUMGF_PK);
#endif
}

// Used to expand PK(verify)
static int expand_pk(gf_t* P22, const uint8_t *pk) {
#if defined(SYMMETRIC)
	gf_t P22c[NUMGF_PK] = {0};
	gf_t* curval = &P22c[0];

	int res = expand_gf(P22c, pk, NUMGF_PK);

	for (int mi = 0; mi < SNOVA_m1; ++mi)
		for (int ni = 0; ni < SNOVA_o; ++ni)
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
#if SNOVA_q == 16
					P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = 0;
#else
					P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = *curval;
#endif
					P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
					curval++;
				}

				for (int nj = ni + 1; nj < SNOVA_o; ++nj)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P22[((mi * SNOVA_o + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
#if SNOVA_q == 16
						P22[((mi * SNOVA_o + nj) * SNOVA_o + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = 0;
#else
						P22[((mi * SNOVA_o + nj) * SNOVA_o + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = *curval;
#endif
						curval++;
					}
			}

	return res;
#else

	return expand_gf(P22, pk, NUMGF_PK);
#endif
}

/**
 * Expand the public key from a seed. Make symmetric
 */
void expand_public(gf_t* P_matrix, const uint8_t *seed) {
	uint8_t pk_bytes[NUM_GEN_PUB_BYTES];

	snova_pk_expand(pk_bytes, NUM_GEN_PUB_BYTES, seed, SEED_LENGTH_PUBLIC);

#if defined(SYMMETRIC)
	gf_t pk_gf[NUM_GEN_PUB_GF];
	convert_bytes_to_GF(pk_gf, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);

	// Make symmetric
	gf_t* P11 = P_matrix;
	gf_t* P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t* P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
	gf_t* abq = P21 + SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2;

	gf_t* curval = &pk_gf[0];

	for (int mi = 0; mi < SNOVA_m1; ++mi)
		for (int ni = 0; ni < SNOVA_v; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
#if SNOVA_q == 16
					P11[((mi * SNOVA_v + ni) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = 0;
#else
					P11[((mi * SNOVA_v + ni) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = *curval;
#endif
					P11[((mi * SNOVA_v + ni) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
					curval++;
				}

			for (int nj = ni + 1; nj < SNOVA_v; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
#if SNOVA_q == 16
						P11[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = 0;
#else
						P11[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = *curval;
#endif
						P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}

			for (int nj = 0; nj < SNOVA_o; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P12[((mi * SNOVA_v + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
#if SNOVA_q == 16
						P21[((mi * SNOVA_o + nj) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = 0;
#else
						P21[((mi * SNOVA_o + nj) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] = *curval;
#endif
						curval++;
					}
		}

	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l); idx++) {
		abq[idx] = curval[idx];
	}
#elif PUB_REORDER
	uint8_t pk_gf[NUM_GEN_PUB_GF];
	convert_bytes_to_GF(pk_gf, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);

	gf_t* P11 = P_matrix;
	gf_t* P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t* P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
	gf_t* curval = &pk_gf[0];

	for (int mi = 0; mi < SNOVA_m1; ++mi) {
		for (int ni = 0; ni < SNOVA_v; ++ni)
			for (int nj = 0; nj < SNOVA_v; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}

		for (int ni = 0; ni < SNOVA_v; ++ni)
			for (int nj = 0; nj < SNOVA_o; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P12[((mi * SNOVA_v + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}

		for (int ni = 0; ni < SNOVA_o; ++ni)
			for (int nj = 0; nj < SNOVA_v; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P21[((mi * SNOVA_o + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}
	}
#else
	convert_bytes_to_GF(P_matrix, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);
#endif
}

static void hash_combined(uint8_t *hash_out, const uint8_t *digest, const size_t len_digest, const uint8_t *pk_seed,
                          const uint8_t *salt) {
	shake_t state;
	shake256_init(&state);
#if HASH_PK
	sha3_update(&state, pk_seed, BYTES_PK_HASH);
#else
	sha3_update(&state, pk_seed, SEED_LENGTH_PUBLIC);
#endif
	sha3_update(&state, digest, len_digest);
	sha3_update(&state, salt, BYTES_SALT);
	shake_xof_final(&state);
	shake_squeeze(&state, hash_out, BYTES_HASH);
}

/**
 * Improve q and calculate Q matrix
 */
static inline void gen_a_FqS(gf_t* Qm, gf_t* q) {
#if ROUND2_T12
	if (!q[SNOVA_l - 1]) {
		q[SNOVA_l - 1] = SNOVA_q - (q[0] + (q[0] == 0));
	}
#endif

	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		gf_t sum = 0;
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_set_add(&sum, gf_mult(q[j1], gf_S[j1 * SNOVA_l2 + i1]));
		}
		Qm[i1] = sum;
	}
}

static inline void gen_a_FqS_sec(gf_t* Qm, gf_t* q) {
#if ROUND2_T12
	if (!q[SNOVA_l - 1]) {
		q[SNOVA_l - 1] = SNOVA_q - (q[0] + (q[0] == 0));
	}
#endif

	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		gf_t sum = 0;
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_set_add(&sum, gf_mult(q[j1], gf_S[j1 * SNOVA_l2 + i1]));
		}
		Qm[i1] = sum;
	}
}

/**
 * Expand T12 matrix and coefficients. Shared by genkey and sign
 */
#define REJECTION_LIMIT ((256 / SNOVA_q) * SNOVA_q)
#define SK_BLOCK_SIZE 32
static void expand_T12(gf_t* T12, const uint8_t *seed) {
	gf_t T12coef[SNOVA_o * SNOVA_v * SNOVA_l];
	gf_t sk_data[SK_BLOCK_SIZE];
	shake_t state;

	shake256_init(&state);
	sha3_update(&state, seed, SEED_LENGTH_PRIVATE);
	shake_xof_final(&state);

	size_t idx = SK_BLOCK_SIZE;
	size_t t_idx = 0;

	while (t_idx < SNOVA_o * SNOVA_v * SNOVA_l) {
		if (idx >= SK_BLOCK_SIZE) {
			shake_squeeze(&state, sk_data, SK_BLOCK_SIZE);
			idx = 0;
		}

#if SNOVA_q != 16
		// Rejection sampling
		if (sk_data[idx] < REJECTION_LIMIT) {
			T12coef[t_idx] = sk_data[idx] % SNOVA_q;
			t_idx++;
		}
#else
		T12coef[t_idx] = sk_data[idx] & 0xf;
		t_idx++;
		T12coef[t_idx] = sk_data[idx] >> 4;
		t_idx++;
#endif

		idx++;
	}

	for (size_t i1 = 0; i1 < SNOVA_o * SNOVA_v; i1++) {
		gen_a_FqS_sec(&T12[i1 * SNOVA_l2], &T12coef[i1 * SNOVA_l]);
	}
}

/**
 * Ensure that a matrix is invertible by adding multiples of S
 */
static inline void be_invertible_by_add_aS(gf_t* mat, const gf_t* orig, const int l1, const int l2) {
	memcpy(mat, orig, l1 * l2);
#if ABQ_ALG2
	if ((l1 == SNOVA_l) && (l2 == SNOVA_l))
		if (gf_mat_det(mat) == 0) {
			for (gf_t f1 = 1; f1 < SNOVA_q; f1++) {
#if SNOVA_l > 1
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					gf_set_add(&mat[i1], gf_mult(f1, gf_S[SNOVA_l2 + i1]));
				}
#else
				mat[0] = 1;
#endif
				if (gf_mat_det(mat) != 0) {
					break;
				}
			}
		}
#endif
}

/**
 * Use last part of the P matrix to establish ABQ
 */
static void gen_ABQ(gf_t* A, gf_t* Am, gf_t* Bm, gf_t* Q1m, gf_t* Q2m) {
	gf_t* B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
	gf_t* q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
	gf_t* q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (size_t idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&Am[idx * SNOVA_r2], &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
		be_invertible_by_add_aS(&Bm[idx * SNOVA_lr], &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);
		gen_a_FqS(&Q1m[idx * SNOVA_l2], &q1[idx * SNOVA_l]);
		gen_a_FqS(&Q2m[idx * SNOVA_l2], &q2[idx * SNOVA_l]);
	}
}

/**
 * Fix the ABQ to constants
 */
#if FIXED_ABQ
uint8_t fixed_abq[SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l)] = {0};

static void gen_fixed_ABQ(const char *abq_seed) {
	uint8_t rng_out[SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l)] = {0};

	shake256(rng_out, SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l), (uint8_t *)abq_seed, strlen(abq_seed));
	convert_bytes_to_GF(fixed_abq, rng_out, SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l));

#if !ROUND2_T12
	// Check if q1 and q2 are always non zero
	gf_t* aptr = fixed_abq;
	gf_t* q1 = aptr + SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr);
	gf_t* q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (int i1 = 0 ; i1 < SNOVA_o * SNOVA_alpha; i1++) {
		uint32_t sum1 = 0;
		uint32_t sum2 = 0;
		for (int j1 = 0 ; j1 < SNOVA_l; j1 ++) {
			sum1 |= q1[i1 * SNOVA_l + j1];
			sum2 |= q2[i1 * SNOVA_l + j1];
		}
		if ((sum1 == 0) || (sum2 == 0)) {
			static int first = 1;
			if (first) {
				printf("Warning: Some q1,q2 are zero for '%s' (%s)\n", CRYPTO_ALGNAME, abq_seed);
				first = 0;
			}
		}
	}
#endif
}
#endif

/**
 * Reference version of genkey.
 */
int SNOVA_NAMESPACE(genkeys)(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
	if (first_time) {
		snova_init();
	}

	/**
	 * Gen T12 matrix
	 */
	gf_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];
	gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2] = {0};

	expand_T12(T12, seed + SEED_LENGTH_PUBLIC);

	/**
	 * Gen Public matrix but not ABQ
	 */
	gf_t P_matrix[NUM_PUB_GF];

	expand_public(P_matrix, seed);

	/**
	 * Calculate F12 matrix, use P11
	 */
	gf_t* P11 = P_matrix;
	gf_t* P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t* P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
	gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2] = {0};

	for (int i1 = 0; i1 < SNOVA_m1; i1++) {
		for (int j1 = 0; j1 < SNOVA_v; j1++) {
			// Squeeze P11[0..SNOVA_v * SNOVA_l2]
			for (int idx = 0; idx < SNOVA_v; idx++) {
				for (int k1 = 0; k1 < SNOVA_o; k1++) {
					// Multiply public and secret
					gf_mat_mul_add(&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
					               &P11[((i1 * SNOVA_v + j1) * SNOVA_v + idx) * SNOVA_l2],
					               &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
				}
			}
		}
	}

	// Use P12
	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		gf_set_add(&F12[i1], P12[i1]);
	}

	for (int i1 = 0; i1 < SNOVA_m1; i1++) {
		for (int idx = 0; idx < SNOVA_v; idx++) {
			for (int k1 = 0; k1 < SNOVA_o; k1++) {
				for (int j1 = 0; j1 < SNOVA_o; j1++) {
					// Multiply two secrets
					gf_mat_mul_add(&P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2], &T12[(idx * SNOVA_o + j1) * SNOVA_l2],
					               &F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);
				}
			}
		}
	}

	/**
	 * Calculate P22. Uses P21
	 */
	for (int i1 = 0; i1 < SNOVA_m1; i1++) {
		for (int j1 = 0; j1 < SNOVA_o; j1++) {
			for (int idx = 0; idx < SNOVA_v; idx++) {
				for (int k1 = 0; k1 < SNOVA_o; k1++) {
					// Multiply public and secret
					gf_mat_mul_add(&P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
					               &P21[((i1 * SNOVA_o + j1) * SNOVA_v + idx) * SNOVA_l2],
					               &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
				}
			}
		}
	}

#if SNOVA_q != 16
	// Negate P22
	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2; i1++) {
		P22[i1] = (SNOVA_q - P22[i1]) % SNOVA_q;
	}
#endif

	/**
	 * Output public and secret keys
	 */
	memcpy(pk, seed, SEED_LENGTH_PUBLIC);
	compress_pk(pk + SEED_LENGTH_PUBLIC, P22);
	memcpy(sk, seed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);

#if HASH_PK
	shake256(sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, BYTES_PK_HASH, pk, BYTES_PK);
#endif

	return 0;
}

/**
 * Dummy sk_expand.
 */
int SNOVA_NAMESPACE(sk_expand)(expanded_SK* skx, const uint8_t *sk) {
	memset(skx, 0, sizeof(expanded_SK));
	memcpy(skx, sk, BYTES_SK);

	return 0;
}

/**
 * Reference version of Sign. Deterministic using the salt provided
 */
int SNOVA_NAMESPACE(sign)(const expanded_SK* skx, uint8_t *sig, const uint8_t *digest, const size_t len_digest,
                          const uint8_t *salt) {
	if (first_time) {
		snova_init();
	}

	const uint8_t *seed = (uint8_t *)skx;

	gf_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];
	expand_T12(T12, seed + SEED_LENGTH_PUBLIC);

	gf_t P_matrix[NUM_PUB_GF];

	expand_public(P_matrix, seed);

	/**
	 * Calculate F12, F21
	 */
	gf_t* P11 = P_matrix;
	gf_t* P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t* P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;

	gf_t F21[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2] = {0};

	for (int i1 = 0; i1 < SNOVA_m1; i1++) {
		for (int j1 = 0; j1 < SNOVA_v; j1++) {
			for (int j2 = 0; j2 < SNOVA_v; j2++) {
				for (int k1 = 0; k1 < SNOVA_o; k1++) {
					// Multiply public and secret
					gf_mat_mul_add(&F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2], &T12[(j2 * SNOVA_o + k1) * SNOVA_l2],
					               &P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2]);
				}
			}
		}
	}

	gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2] = {0};

	for (int i1 = 0; i1 < SNOVA_m1; i1++) {
		for (int j1 = 0; j1 < SNOVA_v; j1++) {
			for (int j2 = 0; j2 < SNOVA_v; j2++) {
				for (int k1 = 0; k1 < SNOVA_o; k1++) {
					// Multiply public and secret
					gf_mat_mul_add(&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
					               &P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2], &T12[(j2 * SNOVA_o + k1) * SNOVA_l2]);
				}
			}
		}
	}

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		gf_set_add(&F12[i1], P12[i1]);
	}

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		gf_set_add(&F21[i1], P21[i1]);
	}

	// Generate ABQ
	gf_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
	gf_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
	gf_t Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
	gf_t Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2];

	gf_t* aptr = P_matrix + SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2;
	gf_t* q1 = aptr + SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr);
	gf_t* q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

#if FIXED_ABQ
	memcpy(aptr, fixed_abq, sizeof(fixed_abq));
#endif
	gen_ABQ(aptr, Am, Bm, Q1, Q2);

	// Calculate message has of size l^2o
	gf_t hash_in_GF16[GF16_HASH];

	uint8_t sign_hashb[BYTES_HASH];
#if HASH_PK
	hash_combined(sign_hashb, digest, len_digest, seed + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, salt);
#else
	hash_combined(sign_hashb, digest, len_digest, seed, salt);
#endif
	expand_gf(hash_in_GF16, sign_hashb, GF16_HASH);

	// Find a solution for T.X
	gf_t gauss[SNOVA_o * SNOVA_lr][SNOVA_o * SNOVA_lr + 1];
	gf_t solution[SNOVA_o * SNOVA_lr] = {0};
	gf_t signature_in_GF[SNOVA_n * SNOVA_lr] = {0};
	int flag_redo = 1;
	uint8_t num_sign = 0;

	do {
		memset(gauss, 0, sizeof(gauss));
		num_sign++;
		if (num_sign == 255) {
			// Probability of getting here is about q^{-255}
			memset(sig, 0, BYTES_SIGNATURE);
			return -1;
		}
		flag_redo = 0;

		// generate the vinegar value
		uint8_t vinegar_in_byte[NUM_GEN_SEC_BYTES];
		shake_t v_instance;

		shake256_init(&v_instance);
		sha3_update(&v_instance, seed + SEED_LENGTH_PUBLIC, SEED_LENGTH_PRIVATE);
		sha3_update(&v_instance, digest, BYTES_DIGEST);
		sha3_update(&v_instance, salt, BYTES_SALT);
		sha3_update(&v_instance, &num_sign, 1);
		shake_xof_final(&v_instance);

		size_t b_idx = SK_BLOCK_SIZE;
		size_t t_idx = 0;

		while (t_idx < SNOVA_v * SNOVA_l * SNOVA_r) {
			if (b_idx >= SK_BLOCK_SIZE) {
				shake_squeeze(&v_instance, vinegar_in_byte, SK_BLOCK_SIZE);
				b_idx = 0;
			}

#if SNOVA_q != 16
			// Rejection sampling
			if (vinegar_in_byte[b_idx] < REJECTION_LIMIT) {
				signature_in_GF[t_idx] = vinegar_in_byte[b_idx] % SNOVA_q;
				t_idx++;
			}
#else
			signature_in_GF[t_idx] = vinegar_in_byte[b_idx] & 0xf;
			t_idx++;
			signature_in_GF[t_idx] = vinegar_in_byte[b_idx] >> 4;
			t_idx++;
#endif
			b_idx++;
		}

		// Calculate Fvv
		gf_t Fvv_in_GF16Matrix[SNOVA_o * SNOVA_lr] = {0};

		/**
		 * Evaluate whipped central map
		 */
		gf_t sum_t0[SNOVA_m1 * SNOVA_l * SNOVA_v * SNOVA_lr] = {0};
		gf_t sum_t1[SNOVA_m1 * SNOVA_l2 * SNOVA_r2] = {0};

		/**
		 * Whip signature
		 */
		gf_t whipped_sig[SNOVA_l * SNOVA_v * SNOVA_lr] = {0};

		for (int ab = 0; ab < SNOVA_l; ++ab)
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_r; j1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							gf_set_add(&whipped_sig[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r + j1],
							           gf_mult(gf_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1],
							                   signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r + j1]));

		// Right
		for (int mi = 0; mi < SNOVA_m1; ++mi) {
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int nj = 0; nj < SNOVA_v; ++nj)
						gf_mat_mul_add_lr(&sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr],
						                  &P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2],
						                  &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);

			// Left, transposed whipped_sig
			for (int a1 = 0; a1 < SNOVA_l; ++a1)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int ni = 0; ni < SNOVA_v; ++ni)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								for (int k1 = 0; k1 < SNOVA_l; k1++)
									gf_set_add(
									    &sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
									    gf_mult(whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_lr + k1 * SNOVA_r + i1],
									            sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr + k1 * SNOVA_r + j1]));
		}

		/**
		 * Apply A, B, q1 and q2, aka E matrix
		 */
		for (int mi = 0; mi < SNOVA_o; ++mi)
			for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
				int mi_prime = i_prime(mi, alpha);

				gf_t gf16m_temp1[SNOVA_r2] = {0};
				gf_t gf16m_temp2[SNOVA_lr] = {0};

				// apply q1 and q2
				for (int a1 = 0; a1 < SNOVA_l; ++a1) {
					gf_t sumb[SNOVA_r2] = {0};

					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								gf_set_add(
								    &sumb[i1 * SNOVA_r + j1],
								    gf_mult(sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
								            q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));

					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int j1 = 0; j1 < SNOVA_r; j1++)
							gf_set_add(&gf16m_temp1[i1 * SNOVA_r + j1],
							           gf_mult(sumb[i1 * SNOVA_r + j1], q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]));
				}

				// A and B
				gf_mat_mul_add_lr(gf16m_temp2, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_r, SNOVA_r,
				                  SNOVA_l);
				gf_mat_mul_add_lr(&Fvv_in_GF16Matrix[mi * SNOVA_lr], &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp2,
				                  SNOVA_r, SNOVA_r, SNOVA_l);
			}

		for (int mi = 0; mi < SNOVA_o; ++mi)
			// Set the last column of gauss matrix
			for (int i1 = 0; i1 < SNOVA_lr; i1++) {
				gauss[mi * SNOVA_lr + i1][SNOVA_o * SNOVA_lr] =
				    gf_sub(hash_in_GF16[mi * SNOVA_lr + i1], Fvv_in_GF16Matrix[mi * SNOVA_lr + i1]);
			}

		// Whipped F21
		gf_t whipped_F21[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr] = {0};
		gf_t whipped_F12[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr] = {0};

		/**
		 * Whipped F21
		 */

		for (int mi = 0; mi < SNOVA_m1; mi++) {
			for (int idx = 0; idx < SNOVA_o; idx++)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int nj = 0; nj < SNOVA_v; ++nj)
						gf_mat_mul_add_lr(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr],
						                  &F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2],
						                  &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);

			// Transpose for F12
			for (int idx = 0; idx < SNOVA_o; idx++)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int nj = 0; nj < SNOVA_v; ++nj)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								for (int k1 = 0; k1 < SNOVA_l; k1++)
									gf_set_add(
									    &whipped_F12[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1],
									    gf_mult(F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1],
									            whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr + k1 * SNOVA_r + j1]));
		}

		// compute the coefficients of Xo and put into gauss matrix and compute
		// the coefficients of Xo^t and add into gauss matrix
		for (int mi = 0; mi < SNOVA_o; mi++) {
			for (int idx = 0; idx < SNOVA_o; idx++) {
				for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
					gf_t gf16m_temp0[SNOVA_l2] = {0};
					gf_t gf16m_temp1[SNOVA_lr] = {0};

					gf_t X_tmp[SNOVA_l2] = {0};
					int mi_prime = i_prime(mi, alpha);

					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_lr; i1++)
							gf_set_add(&gf16m_temp1[i1],
							           gf_mult(whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1],
							                   q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));
					gf_mat_mul_add_lr(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r,
					                  SNOVA_l);
					gf_mat_mul_add(X_tmp, &Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2], gf16m_temp0);

					for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
						for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
							for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
								for (int tj2 = 0; tj2 < SNOVA_r; tj2++) {
									int ti = ti1 * SNOVA_l + ti2;
									int tj = tj1 * SNOVA_r + tj2;
									gf_set_add(&gauss[mi * SNOVA_lr + ti][idx * SNOVA_lr + tj],
									           gf_mult(X_tmp[tj1 * SNOVA_l + ti2],
									                   Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r + tj2]));
								}
				}
			}
		}

		for (int mi = 0; mi < SNOVA_o; mi++) {
			for (int idx = 0; idx < SNOVA_o; idx++) {
				for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
					gf_t gf16m_temp0[SNOVA_lr] = {0};
					gf_t gf16m_temp1[SNOVA_lr] = {0};

					gf_t X_tmp[SNOVA_lr] = {0};
					int mi_prime = i_prime(mi, alpha);

					// Transpose
					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								gf_set_add(
								    &gf16m_temp1[j1 * SNOVA_l + i1],
								    gf_mult(
								        whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1],
								        q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));

					gf_mat_mul_add_lr(gf16m_temp0, &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp1, SNOVA_r, SNOVA_r,
					                  SNOVA_l);
					gf_mat_mul_add_lr(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l,
					                  SNOVA_l);

					for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
						for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
							for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
								for (int tj2 = 0; tj2 < SNOVA_r; tj2++) {
									int ti = ti1 * SNOVA_l + ti2;
									int tj = tj1 * SNOVA_r + tj2;
									gf_set_add(&gauss[mi * SNOVA_lr + ti][idx * SNOVA_lr + tj],
									           gf_mult(X_tmp[ti1 * SNOVA_l + tj1],
									                   Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2]));
								}
				}
			}
		}

		// Gaussian elimination
		for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) {
			if (gauss[i][i] == 0) {
				for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
					if (gauss[j][i] != 0) {
						for (int k = i; k < SNOVA_o * SNOVA_lr + 1; ++k) {
							gf_t t_GF16 = gauss[i][k];
							gauss[i][k] = gauss[j][k];
							gauss[j][k] = t_GF16;
						}
						break;
					}
				}
			}
			if (gauss[i][i] == 0) {
				flag_redo = 1;
				break;
			}

			gf_t t_GF16 = gf_inv(gauss[i][i]);
			for (int k = i; k < SNOVA_o * SNOVA_lr + 1; ++k) {
				gauss[i][k] = gf_mult(gauss[i][k], t_GF16);
			}

			for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
				if (gauss[j][i] != 0) {
					gf_t gji = gauss[j][i];
					for (int k = i; k < SNOVA_o * SNOVA_lr + 1; ++k) {
						gauss[j][k] = gf_sub(gauss[j][k], gf_mult(gauss[i][k], gji));
					}
				}
			}
		}

		if (!flag_redo) {
			// Last step of Gaussian elimination
			memset(solution, 0, sizeof(solution));

			for (int i = SNOVA_o * SNOVA_lr - 1; i >= 0; --i) {
				gf_t sum = 0;
				for (int k = i + 1; k < SNOVA_o * SNOVA_lr; ++k) {
					gf_set_add(&sum, gf_mult(gauss[i][k], solution[k]));
				}
				solution[i] = gf_sub(gauss[i][SNOVA_o * SNOVA_lr], sum);
			}
			memcpy(signature_in_GF + SNOVA_v * SNOVA_lr, solution, SNOVA_o * SNOVA_lr);

			// Establish signature using T12
			for (int index = 0; index < SNOVA_v; ++index) {
				for (int mi = 0; mi < SNOVA_o; ++mi) {
					gf_mat_mul_add_lr(&signature_in_GF[index * SNOVA_lr], &T12[(index * SNOVA_o + mi) * SNOVA_l2],
					                  &solution[mi * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);
				}
			}
		}
	} while (flag_redo);

	compress_gf(sig, signature_in_GF, SNOVA_n * SNOVA_lr);
	memcpy(sig + BYTES_SIGNATURE - BYTES_SALT, salt, BYTES_SALT);

	return 0;
}

/**
 * PK expansion.
 */
int SNOVA_NAMESPACE(pk_expand)(expanded_PK* pkx, const uint8_t *pk) {
	if (first_time) {
		snova_init();
	}

	memset(pkx, 0, sizeof(expanded_PK));
	memcpy(pkx->pk_seed, pk, SEED_LENGTH_PUBLIC);
#if HASH_PK
	shake256(pkx->pk_hash, BYTES_PK_HASH, pk, BYTES_PK);
#endif

	/**
	 * Create P matrix
	 */
	gf_t P_matrix[NUM_PUB_GF];
	gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];

	gf_t* P11 = P_matrix;
	gf_t* P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t* P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;

	if (expand_pk(P22, pk + SEED_LENGTH_PUBLIC)) {
		return -1;
	}
	expand_public(P_matrix, pk);

	for (int mi = 0; mi < SNOVA_m1; ++mi) {
		for (int ni = 0; ni < SNOVA_v; ++ni) {
			for (int nj = 0; nj < SNOVA_v; ++nj) {
				for (int idx = 0; idx < SNOVA_l2; idx++)
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + idx] =
					    P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + idx];
			}

			for (int nj = SNOVA_v; nj < SNOVA_n; ++nj) {
				for (int idx = 0; idx < SNOVA_l2; idx++)
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + idx] =
					    P12[((mi * SNOVA_v + ni) * SNOVA_o + (nj - SNOVA_v)) * SNOVA_l2 + idx];
			}
		}
		for (int ni = SNOVA_v; ni < SNOVA_n; ++ni) {
			for (int nj = 0; nj < SNOVA_v; ++nj) {
				for (int idx = 0; idx < SNOVA_l2; idx++)
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + idx] =
					    P21[((mi * SNOVA_o + (ni - SNOVA_v)) * SNOVA_v + nj) * SNOVA_l2 + idx];
			}

			for (int nj = SNOVA_v; nj < SNOVA_n; ++nj) {
				for (int idx = 0; idx < SNOVA_l2; idx++)
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + idx] =
					    P22[((mi * SNOVA_o + (ni - SNOVA_v)) * SNOVA_o + nj - SNOVA_v) * SNOVA_l2 + idx];
			}
		}
	}

	/**
	 * Create AB matrices, improve q
	 */
	gf_t* A = P_matrix + (SNOVA_m1 * (SNOVA_n * SNOVA_n - SNOVA_o * SNOVA_o)) * SNOVA_l2;
	gf_t* B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
	gf_t* q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
	gf_t* q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

#if FIXED_ABQ
	memcpy(A, fixed_abq, sizeof(fixed_abq));
#endif

	for (size_t idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&(pkx->Am[idx * SNOVA_r2]), &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
		be_invertible_by_add_aS(&(pkx->Bm[idx * SNOVA_lr]), &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);

#if ROUND2_T12
		if (!q1[idx * SNOVA_l + SNOVA_l - 1]) {
			q1[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q1[idx * SNOVA_l] + (q1[idx * SNOVA_l] == 0));
		}
		if (!q2[idx * SNOVA_l + SNOVA_l - 1]) {
			q2[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q2[idx * SNOVA_l] + (q2[idx * SNOVA_l] == 0));
		}
#endif
	}

	memcpy(pkx->q1, q1, SNOVA_o * SNOVA_alpha * SNOVA_l);
	memcpy(pkx->q2, q2, SNOVA_o * SNOVA_alpha * SNOVA_l);

	return 0;
}

/**
 * Reference version of verify.
 */
int SNOVA_NAMESPACE(verify)(const expanded_PK* pkx, const uint8_t *sig, const uint8_t *digest, const size_t len_digest) {
	if (first_time) {
		snova_init();
	}

	// Actual verify
	gf_t signature_in_GF[NUMGF_SIGNATURE];
	if (expand_gf(signature_in_GF, sig, NUMGF_SIGNATURE)) {
		return -1;
	}

	/**
	 * Evaluate whipped central map
	 */
	gf_t hash_in_GF[SNOVA_o * SNOVA_lr] = {0};
	gf_t sum_t0[SNOVA_m1 * SNOVA_l * SNOVA_n * SNOVA_lr] = {0};
	gf_t sum_t1[SNOVA_m1 * SNOVA_l2 * SNOVA_r2] = {0};

	/**
	 * Whip signature
	 */
	gf_t whipped_sig[SNOVA_l * SNOVA_n * SNOVA_lr] = {0};

	for (int ab = 0; ab < SNOVA_l; ++ab)
		for (int idx = 0; idx < SNOVA_n; ++idx)
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = 0; j1 < SNOVA_r; j1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						gf_set_add(&whipped_sig[(ab * SNOVA_n + idx) * SNOVA_lr + i1 * SNOVA_r + j1],
						           gf_mult(gf_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1],
						                   signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r + j1]));

	for (int mi = 0; mi < SNOVA_m1; ++mi) {
		// Right
		for (int ni = 0; ni < SNOVA_n; ++ni)
			for (int b1 = 0; b1 < SNOVA_l; ++b1)
				for (int nj = 0; nj < SNOVA_n; ++nj)
					gf_mat_mul_add16(&sum_t0[((mi * SNOVA_l + b1) * SNOVA_n + ni) * SNOVA_lr],
					                 &pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2],
					                 &whipped_sig[(b1 * SNOVA_n + nj) * SNOVA_lr]);

		// Left, transposed whipped_sig
		for (int a1 = 0; a1 < SNOVA_l; ++a1)
			for (int b1 = 0; b1 < SNOVA_l; ++b1)
				for (int ni = 0; ni < SNOVA_n; ++ni)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int j1 = 0; j1 < SNOVA_r; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								gf_set_add(
								    &sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
								    gf_mult(whipped_sig[(a1 * SNOVA_n + ni) * SNOVA_lr + k1 * SNOVA_r + i1],
								            sum_t0[((mi * SNOVA_l + b1) * SNOVA_n + ni) * SNOVA_lr + k1 * SNOVA_r + j1]));
	}

	/**
	 * Apply A, B, q1 and q2, aka E matrix
	 */
	for (int mi = 0; mi < SNOVA_o; ++mi)
		for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
			int mi_prime = i_prime(mi, alpha);

			gf_t gf16m_temp1[SNOVA_r2] = {0};
			gf_t gf16m_temp2[SNOVA_lr] = {0};

			// apply q1 and q2
			for (int a1 = 0; a1 < SNOVA_l; ++a1) {
				gf_t sumb[SNOVA_r2] = {0};

				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int j1 = 0; j1 < SNOVA_r; j1++)
							gf_set_add(&sumb[i1 * SNOVA_r + j1],
							           gf_mult(sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
							                   pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));

				for (int i1 = 0; i1 < SNOVA_r; i1++)
					for (int j1 = 0; j1 < SNOVA_r; j1++)
						gf_set_add(&gf16m_temp1[i1 * SNOVA_r + j1],
						           gf_mult(sumb[i1 * SNOVA_r + j1], pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]));
			}

			// A and B
			gf_mat_mul_add_lr(gf16m_temp2, gf16m_temp1, &(pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr]), SNOVA_r, SNOVA_r,
			                  SNOVA_l);
			gf_mat_mul_add_lr(&hash_in_GF[mi * SNOVA_lr], &(pkx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2]), gf16m_temp2,
			                  SNOVA_r, SNOVA_r, SNOVA_l);
		}

	/**
	 * Check hashes
	 */
	uint8_t signed_bytes[BYTES_HASH];
	uint8_t signed_gf[GF16_HASH] = {0};
	const uint8_t *salt = sig + BYTES_SIGNATURE - BYTES_SALT;
#if HASH_PK
	hash_combined(signed_bytes, digest, len_digest, pkx->pk_hash, salt);
#else
	hash_combined(signed_bytes, digest, len_digest, pkx->pk_seed, salt);
#endif
	expand_gf(signed_gf, signed_bytes, GF16_HASH);

	int result = 0;
	for (int i = 0; i < GF16_HASH; ++i) {
		if (hash_in_GF[i] != signed_gf[i]) {
			result = -1;
			break;
		}
	}

	return result;
}
