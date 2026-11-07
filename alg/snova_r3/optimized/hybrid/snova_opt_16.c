// SPDX-License-Identifier: MIT

/**
 * Optimized q=16 implementation. Uses GFNI if available
 *
 * Copyright (c) 2026 SNOVA TEAM
 */

#include <stdalign.h>
#include <stdint.h>
#include <string.h>

#include "snova.h"
#include "symmetric.h"
#ifndef R3_FULL_KEY
#define R3_FULL_KEY 1
#endif
extern int r3_port_genkeys(uint8_t *,uint8_t *,const uint8_t *);
extern int r3_port_sk_expand(void *,const uint8_t *);

#if SNOVA_q != 16
#error "SNOVA_q != 16"
#include "stop"
#endif

#ifndef USE_AVX2
#define USE_AVX2 __AVX2__
#endif
#define USE_AVX2_VERIFY (USE_AVX2 && (SNOVA_l != SNOVA_r))

#ifndef USE_GFNI
#define USE_GFNI (USE_AVX2 && __GFNI__ && ((SNOVA_l != SNOVA_r) || (SNOVA_l == 5)))
#endif

#if USE_AVX2
#include <immintrin.h>
#endif

#define SNOVA_olr16 (SNOVA_o * SNOVA_lr / 16 + 1)
#define SNOVA_olr (SNOVA_olr16 * 16)

#if USE_GFNI
#define SNOVA_lr16 ((SNOVA_lr + 31) / 32)
#define SNOVA_lr32 (SNOVA_lr16 * 32)
#else
#define SNOVA_lr16 ((SNOVA_lr + 15) / 16)
#define SNOVA_lr32 (SNOVA_lr16 * 16)
#endif

#define SNOVA_vl16 ((SNOVA_v * SNOVA_l + 15) / 16)
#define SNOVA_vl (SNOVA_vl16 * 16)

#define SNOVA_mvl16 ((SNOVA_m1 * SNOVA_v * SNOVA_l + 15) / 16)
#define SNOVA_mvl (SNOVA_mvl16 * 16)

#if 0
static inline uint16_t gf16_expand(const uint8_t a) {
	uint16_t val = a | (a << 3) | (a << 6) | (a << 9);
	return val & 0x1111;
}

static inline uint16_t gf16_compress(const uint16_t a) {
	uint16_t val = (a & 0xf) ^ ((a & 0xf0) >> 3) ^ ((a & 0xf00) >> 6) ^ ((a & 0xf000) >> 9);
	return (val ^ ((val & 0xf0) >> 3) ^ (val >> 4)) & 0xf;
}

#if USE_GFNI
static inline uint8_t gfni_cleanup(const uint8_t val) {
	return (val ^ ((val & 0xf0) >> 3) ^ (val >> 4)) & 0xf;
}
#endif

/**
 * Constant time function. CT is according to valgrind
 */
static inline uint32_t ct_is_not_zero(uint8_t val) {
	// return (val | (val >> 1) | (val >> 2) | (val >> 3)) & 1;
	return val != 0;
}
#else

#define gf16_expand(a) (((a) | ((a) << 3) | ((a) << 6) | ((a) << 9)) & 0x1111)
#define gf16_compress(a) (((a) & 0xf) ^ (((a) >> 3) & 0xe) ^ (((a) >> 6) & 0xe) ^ (((a) >> 7) & 0x1) ^ (((a) >> 9) & 0xe) ^ (((a) >> 10) & 0x3) ^ (((a) >> 12) & 0xe) ^ (((a) >> 13) & 0x7))
#define gfni_cleanup(val) (((val) ^ (((val) & 0xf0) >> 3) ^ ((val) >> 4)) & 0xf)

#define ct_is_not_zero(val) ((val) != 0)

#endif

/**
 * Initialization
 */
static int first_time = 1;

#if FIXED_ABQ
#include "round3_neon.h"

static inline void gen_fixed_ABQ(const char *abq_seed);
#endif

#if SNOVA_l == 4
static uint8_t gf_multtab[SNOVA_q * SNOVA_q] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
	13, 14, 15, 0, 2, 4, 6, 8, 10, 12, 14, 3, 1, 7, 5, 11, 9, 15, 13, 0, 3, 6, 5, 12, 15, 10, 9, 11, 8,
	13, 14, 7, 4, 1, 2, 0, 4, 8, 12, 3, 7, 11, 15, 6, 2, 14, 10, 5, 1, 13, 9, 0, 5, 10, 15, 7, 2, 13,
	8, 14, 11, 4, 1, 9, 12, 3, 6, 0, 6, 12, 10, 11, 13, 7, 1, 5, 3, 9, 15, 14, 8, 2, 4, 0, 7, 14, 9,
	15, 8, 1, 6, 13, 10, 3, 4, 2, 5, 12, 11, 0, 8, 3, 11, 6, 14, 5, 13, 12, 4, 15, 7, 10, 2, 9, 1, 0,
	9, 1, 8, 2, 11, 3, 10, 4, 13, 5, 12, 6, 15, 7, 14, 0, 10, 7, 13, 14, 4, 9, 3, 15, 5, 8, 2, 1, 11,
	6, 12, 0, 11, 5, 14, 10, 1, 15, 4, 7, 12, 2, 9, 13, 6, 8, 3, 0, 12, 11, 7, 5, 9, 14, 2, 10, 6, 1,
	13, 15, 3, 4, 8, 0, 13, 9, 4, 1, 12, 8, 5, 2, 15, 11, 6, 3, 14, 10, 7, 0, 14, 15, 1, 13, 3, 2, 12,
	9, 7, 6, 8, 4, 10, 11, 5, 0, 15, 13, 2, 9, 6, 4, 11, 1, 14, 12, 3, 8, 7, 5, 10
};

static uint16_t gf_Sx[SNOVA_l * SNOVA_l2] = {
	0x1, 0x0, 0x0, 0x0, 0x0, 0x1, 0x0, 0x0, 0x0, 0x0, 0x1, 0x0, 0x0, 0x0, 0x0, 0x1,
	0x1000, 0x111, 0x110, 0x101, 0x111, 0x110, 0x101, 0x100, 0x110, 0x101, 0x100, 0x11, 0x101, 0x100, 0x11, 0x10,
	0x1111, 0x110, 0x1001, 0x1, 0x110, 0x0, 0x111, 0x100, 0x1001, 0x111, 0x11, 0x0, 0x1, 0x100, 0x0, 0x0,
	0x110, 0x11, 0x1000, 0x111, 0x11, 0x1010, 0x100, 0x1100, 0x1000, 0x100, 0x111, 0x1, 0x111, 0x1100, 0x1, 0x110
};

static inline void snova_init(void) {
#if FIXED_ABQ
	gen_fixed_ABQ("SNOVA_ABQ");
#endif
}

#else

static uint8_t gf_multtab[SNOVA_q * SNOVA_q] = {0};
static uint8_t gf_invtab[SNOVA_q] = {0};
static uint8_t gf_addtab[SNOVA_q * SNOVA_q] = {0};
static uint8_t gf_S[SNOVA_l * SNOVA_l2] = {0};
static uint16_t gf_Sx[SNOVA_l * SNOVA_l2] = {0};

static inline uint8_t gf_mult(const uint8_t a, const uint8_t b) {
	return gf_multtab[a * SNOVA_q + b];
}

#if SNOVA_l == 5
static inline uint8_t gf_add(const uint8_t a, const uint8_t b) {
	return gf_addtab[a * SNOVA_q + b];
}
#endif

static inline void gf_set_add(uint8_t *a, const uint8_t b) {
	*a = gf_addtab[*a * SNOVA_q + b];
}

static inline uint8_t gf_sub(const uint8_t a, const uint8_t b) {
#if SNOVA_q != 16
	return gf_addtab[a * SNOVA_q + (SNOVA_q - b) % SNOVA_q];
#else
	return gf_addtab[a * SNOVA_q + b];
#endif
}

static inline void init_uint8_tables(void) {
	// GF(16)
	uint8_t F_star[15] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9}; // Z2[x]/(x^4+x+1)
	for (int i1 = 0; i1 < 16; i1++) {
		gf_multtab[i1] = 0;
		gf_multtab[i1 * SNOVA_q] = 0;
	}
	for (int i1 = 0; i1 < SNOVA_q - 1; i1++)
		for (int j1 = 0; j1 < SNOVA_q - 1; j1++) {
			gf_multtab[F_star[i1] * SNOVA_q + F_star[j1]] = F_star[(i1 + j1) % (SNOVA_q - 1)];
		}

	for (int i1 = 0; i1 < SNOVA_q; i1++)
		for (int j1 = 0; j1 < SNOVA_q; j1++) {
			gf_addtab[i1 * SNOVA_q + j1] = (i1 ^ j1);
		}
	// Use that x^q = x and therefore x^(q-2) = x^-1
	for (int i1 = 0; i1 < SNOVA_q; i1++) {
		uint8_t val = i1;
		for (int j1 = 3; j1 < SNOVA_q; j1++) {
			val = gf_mult(val, i1);
		}
		gf_invtab[i1] = val;
	}
}

static inline void gf_mat_mul(uint8_t *a, const uint8_t *b, const uint8_t *c) {
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			uint8_t sum = 0;
			for (int k1 = 0; k1 < SNOVA_l; k1++) {
				gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
			}
			a[i1 * SNOVA_l + j1] = sum;
		}
}

// Set the irreducible S matrix
static inline void set_S(uint8_t *gf_S1) {
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_S1[i1 * SNOVA_l + j1] = 8 - (i1 + j1);
		}
#if SNOVA_l == 5
	gf_S1[SNOVA_l2 - 1] = 9;
#endif
}

static inline void gen_S_array(void) {
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

	for (int i1 = 0; i1 < SNOVA_l * SNOVA_l2; i1++) {
		gf_Sx[i1] = gf16_expand(gf_S[i1]);
	}
}

static inline void snova_init(void) {
	init_uint8_tables();
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

#endif

#define SNOVA_INIT      \
    if (first_time)     \
    {                   \
        first_time = 0; \
        snova_init();   \
    }

/**
 * Utilities
 */

#if SNOVA_l == 4

static inline uint8_t gf_mat_det(uint8_t *a) {
#define DET_SUB(a, b) (a ^ b)
#define DET_MULT(a, b) gf_multtab[a * SNOVA_q + b]
	uint8_t det = 0;
	uint8_t det_l;
	uint8_t det_r;
#define DET_L(x, y) det_l = DET_SUB(DET_MULT(a[x], a[4 + y]), DET_MULT(a[y], a[4 + x]))
#define DET_R(x, y) det_r = DET_SUB(DET_MULT(a[8 + x], a[12 + y]), DET_MULT(a[8 + y], a[12 + x]))
#define DET22(x1, y1, x2, y2) \
    DET_L(x1, y1);            \
    DET_R(x2, y2);            \
    det ^= DET_MULT(det_l, det_r)
	DET22(0, 1, 2, 3);
	DET22(0, 2, 3, 1);
	DET22(0, 3, 1, 2);
	DET22(1, 2, 0, 3);
	DET22(1, 3, 2, 0);
	DET22(2, 3, 0, 1);

	return det;
}

#else
static inline uint8_t gf_mat_det(uint8_t *a) {
	uint8_t det = 0;
#if SNOVA_l == 1
	det = a[0];
#elif SNOVA_l == 2
	det = gf_sub(gf_mult(a[0], a[3]), gf_mult(a[1], a[2]));
#elif SNOVA_l == 3
	det = gf_mult(a[0], gf_sub(gf_mult(a[4], a[8]), gf_mult(a[5], a[7])));
	gf_set_add(&det, gf_mult(a[1], gf_sub(gf_mult(a[5], a[6]), gf_mult(a[3], a[8]))));
	gf_set_add(&det, gf_mult(a[2], gf_sub(gf_mult(a[3], a[7]), gf_mult(a[4], a[6]))));
#elif SNOVA_l == 5
	uint8_t det_l;
	uint8_t det_r;
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
}
#endif

static inline void convert_bytes_to_GF(uint8_t *gf_array, const uint8_t *byte_array, size_t num) {
	for (size_t idx = 0; idx < num / 2; idx++) {
		gf_array[2 * idx] = byte_array[idx] & 0xf;
		gf_array[2 * idx + 1] = byte_array[idx] >> 4;
	}
	if (num & 1) {
		gf_array[num - 1] = byte_array[num / 2] & 0xf;
	}
}

// Used to compress PK (genkey) and SIG(sign)
static inline void compress_gf(uint8_t *byte_array, const uint8_t *gf_array, size_t num) {
	for (size_t idx = 0; idx < num / 2; idx++) {
		byte_array[idx] = gf_array[2 * idx] ^ (gf_array[2 * idx + 1] << 4);
	}
	if (num & 1) {
		byte_array[num / 2] = gf_array[num - 1];
	}
}

/**
 * Expand the public key from a seed.
 */
static inline void expand_public(uint8_t *P_matrix, const uint8_t *seed) {
	uint8_t pk_bytes[NUM_GEN_PUB_BYTES];

	snova_pk_expand(pk_bytes, NUM_GEN_PUB_BYTES, seed, SEED_LENGTH_PUBLIC);

#if PUB_REORDER
	uint8_t pk_gf[NUM_GEN_PUB_GF];
	convert_bytes_to_GF(pk_gf, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);

	uint8_t *P11 = P_matrix;
	uint8_t *P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	uint8_t *P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
	uint8_t *curval = &pk_gf[0];

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

static inline void hash_combined(uint8_t *hash_out, const uint8_t *m, size_t mlen, const uint8_t *pk_seed,
                                 const uint8_t *salt) {
	shake_t state;
	shake256_init(&state);
#if HASH_PK
	shake_absorb(&state, pk_seed, BYTES_PK_HASH);
#else
	shake_absorb(&state, pk_seed, SEED_LENGTH_PUBLIC);
#endif
	shake_absorb(&state, m, mlen);
	shake_absorb(&state, salt, BYTES_SALT);
	shake_finalize(&state);
	shake_squeeze(hash_out, BYTES_HASH, &state);
}

/**
 * Improve q and calculate Q matrix
 */
static inline void gen_a_FqS(uint8_t *Qm, uint8_t *q) {
#if ROUND2_T12
	int16_t not_zero = -ct_is_not_zero(q[SNOVA_l - 1]);
	q[SNOVA_l - 1] = (not_zero & q[SNOVA_l - 1]) | ((not_zero ^ -1) & (SNOVA_q - (q[0] + 1 - ct_is_not_zero(q[0]))));
#endif

	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		uint16_t sum = 0;
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			sum ^= q[j1] * gf_Sx[j1 * SNOVA_l2 + i1];
		}
		Qm[i1] = gf16_compress(sum);
	}
}

static inline void gen_a_FqS_sec(uint8_t *Qm, uint8_t *q) {
#if ROUND2_T12
	int16_t not_zero = -ct_is_not_zero(q[SNOVA_l - 1]);
	q[SNOVA_l - 1] = (not_zero & q[SNOVA_l - 1]) | ((not_zero ^ -1) & (SNOVA_q - (q[0] + 1 - ct_is_not_zero(q[0]))));
#endif

	uint16_t sum;
	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		sum = 0;
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			sum ^= q[j1] * gf_Sx[j1 * SNOVA_l2 + i1];
		}
		Qm[i1] = gf16_compress(sum);
	}
	SNOVA_SECURE_CLEARI(sum);
}

/**
 * Expand T12 matrix and coefficients. Shared by genkey and sign
 */
#define SK_BLOCK_SIZE ((SNOVA_o * SNOVA_v * SNOVA_l + 1) / 2)
static inline void expand_T12(uint8_t *T12, const uint8_t *seed) {
	uint8_t T12coef[SNOVA_o * SNOVA_v * SNOVA_l];
	uint8_t sk_data[(SNOVA_o * SNOVA_v * SNOVA_l + 1) / 2];

	shake256(sk_data, (SNOVA_o * SNOVA_v * SNOVA_l + 1) / 2, seed, SEED_LENGTH_PRIVATE);
	convert_bytes_to_GF(T12coef, sk_data, SNOVA_o * SNOVA_v * SNOVA_l);

	for (size_t i1 = 0; i1 < SNOVA_o * SNOVA_v; i1++) {
		gen_a_FqS_sec(&T12[i1 * SNOVA_l2], &T12coef[i1 * SNOVA_l]);
	}

	SNOVA_SECURE_CLEAR(sk_data);
	SNOVA_SECURE_CLEAR(T12coef);
}

/**
 * Ensure that a matrix is invertible by adding multiples of S
 */
static inline void be_invertible_by_add_aS(uint8_t *mat, const uint8_t *orig, const int l1, const int l2) {
	memcpy(mat, orig, l1 * l2);
#if ABQ_ALG2
#if SNOVA_l > 1
	if ((l1 == SNOVA_l) && (l2 == SNOVA_l)) {
		for (uint8_t f1 = 1; gf_mat_det(mat) == 0; f1++)
			for (int i1 = 0; i1 < SNOVA_l2; i1++) {
				mat[i1] ^= gf16_compress(f1 * gf_Sx[SNOVA_l2 + i1]);
			}
	}
#endif
#endif
}

/**
 * Fix the ABQ to constants
 */
#if FIXED_ABQ
static uint8_t fixed_abq[SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l)] = {0};

static uint16_t fixedAm[SNOVA_o * SNOVA_alpha * SNOVA_r2];
static uint16_t fixedBm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
static uint16_t fixedQ1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
static uint16_t fixedQ2[SNOVA_o * SNOVA_alpha * SNOVA_l2];
static uint16_t fixedq1[SNOVA_o * SNOVA_alpha * SNOVA_l];
static uint16_t fixedq2[SNOVA_o * SNOVA_alpha * SNOVA_l];

static inline void gen_fixed_ABQ(const char *abq_seed) {
	uint8_t rng_out[SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l)] = {0};

	shake256(rng_out, SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l), (uint8_t *)abq_seed, strlen(abq_seed));
	convert_bytes_to_GF(fixed_abq, rng_out, SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l));

	uint8_t *A = fixed_abq;
	uint8_t *B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
	uint8_t *q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
	uint8_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	uint8_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
	uint8_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
	uint8_t Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
	uint8_t Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2];

	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&(Am[idx * SNOVA_r2]), &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
		be_invertible_by_add_aS(&(Bm[idx * SNOVA_lr]), &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);
		gen_a_FqS(&(Q1[idx * SNOVA_l2]), &q1[idx * SNOVA_l]);
		gen_a_FqS(&(Q2[idx * SNOVA_l2]), &q2[idx * SNOVA_l]);

#if ROUND2_T12
		if (!q1[idx * SNOVA_l + SNOVA_l - 1]) {
			q1[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q1[idx * SNOVA_l] + (q1[idx * SNOVA_l] == 0));
		}
		if (!q2[idx * SNOVA_l + SNOVA_l - 1]) {
			q2[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q2[idx * SNOVA_l] + (q2[idx * SNOVA_l] == 0));
		}
#endif
	}

	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_r2; ++idx) {
		fixedAm[idx] = gf16_expand(Am[idx]);
	}
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_lr; ++idx) {
		fixedBm[idx] = gf16_expand(Bm[idx]);
	}
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_l2; ++idx) {
		fixedQ1[idx] = gf16_expand(Q1[idx]);
		fixedQ2[idx] = gf16_expand(Q2[idx]);
	}
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_l; ++idx) {
		fixedq1[idx] = gf16_expand(q1[idx]);
		fixedq2[idx] = gf16_expand(q2[idx]);
	}
}
#else

/**
 * Use last part of the P matrix to establish ABQ
 */
static inline void gen_ABQ(uint8_t *A, uint8_t *Am, uint8_t *Bm, uint8_t *Q1m, uint8_t *Q2m) {
	uint8_t *B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
	uint8_t *q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
	uint8_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (size_t idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&Am[idx * SNOVA_r2], &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
		be_invertible_by_add_aS(&Bm[idx * SNOVA_lr], &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);
		gen_a_FqS(&Q1m[idx * SNOVA_l2], &q1[idx * SNOVA_l]);
		gen_a_FqS(&Q2m[idx * SNOVA_l2], &q2[idx * SNOVA_l]);
	}
}
#endif

/**
 * Optimized version of genkey.
 */
int SNOVA_NAMESPACE(genkeys)(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
#if defined(__aarch64__) && R3_FULL_KEY
return r3_port_genkeys(pk,sk,seed);
#else

	SNOVA_INIT

	uint16_t T12x[SNOVA_o * SNOVA_v * SNOVA_l2];
	{
		uint8_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];

		expand_T12(T12, seed + SEED_LENGTH_PUBLIC);
		for (int i1 = 0; i1 < SNOVA_o * SNOVA_v * SNOVA_l2; i1++) {
			T12x[i1] = gf16_expand(T12[i1]);
		}

		SNOVA_SECURE_CLEAR(T12);
	}

	/**
	 * Gen Public matrix but not ABQ
	 */
	uint8_t P_matrix[NUM_PUB_GF];

	expand_public(P_matrix, seed);

	/**
	 * Calculate F12 matrix, use P11
	 */
	uint8_t *P11 = P_matrix;
	uint8_t *P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	uint8_t *P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
	uint16_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int ni = 0; ni < SNOVA_o; ni++)
				for (int nk = 0; nk < SNOVA_v; nk++) {
					uint16_t P11x[SNOVA_l2];
					for (int i1 = 0; i1 < SNOVA_l2; i1++) {
						P11x[i1] = P11[((mi * SNOVA_v + nj) * SNOVA_v + nk) * SNOVA_l2 + i1];
					}
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								F12[((mi * SNOVA_v + nj) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    T12x[(nk * SNOVA_o + ni) * SNOVA_l2 + k1 * SNOVA_l + j1] * P11x[i1 * SNOVA_l + k1];
				}

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		F12[i1] = gf16_compress(F12[i1]) ^ P12[i1];
	}

	uint16_t P22x[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int ni = 0; ni < SNOVA_o; ni++)
				for (int nk = 0; nk < SNOVA_o; nk++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								P22x[((mi * SNOVA_o + nk) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    T12x[(nj * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    F12[((mi * SNOVA_v + nj) * SNOVA_o + ni) * SNOVA_l2 + k1 * SNOVA_l + j1];

	/**
	 * Calculate P22. Uses P21
	 */
	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nk = 0; nk < SNOVA_o; nk++)
			for (int nj = 0; nj < SNOVA_v; nj++) {
				uint16_t P21x[SNOVA_l2];
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					P21x[i1] = P21[((mi * SNOVA_o + nk) * SNOVA_v + nj) * SNOVA_l2 + i1];
				}
				for (int ni = 0; ni < SNOVA_o; ni++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								P22x[((mi * SNOVA_o + nk) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    P21x[i1 * SNOVA_l + k1] * T12x[(nj * SNOVA_o + ni) * SNOVA_l2 + k1 * SNOVA_l + j1];
			}

	uint8_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2; i1++) {
		P22[i1] = gf16_compress(P22x[i1]);
	}

	/**
	 * Output public and secret keys
	 */
	memcpy(pk, seed, SEED_LENGTH_PUBLIC);
	compress_gf(pk + SEED_LENGTH_PUBLIC, P22, NUMGF_PK);
	memcpy(sk, seed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);

#if HASH_PK
	shake256(sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, BYTES_PK_HASH, pk, BYTES_PK);
#endif

	SNOVA_SECURE_CLEAR(F12);
	SNOVA_SECURE_CLEAR(T12x);

	return 0;

#endif
}

/**
 * SK expansion.
 */
int SNOVA_NAMESPACE(sk_expand)(expanded_SK *skx, const uint8_t *sk) {
#if defined(__aarch64__) && R3_FULL_KEY
return r3_port_sk_expand(skx,sk);
#else

	SNOVA_INIT

	memcpy(skx->pk_seed, sk, SEED_LENGTH_PUBLIC);
	memcpy(skx->sk_seed, sk + SEED_LENGTH_PUBLIC, SEED_LENGTH_PRIVATE);
#if HASH_PK
	memcpy(skx->pk_hash, sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, BYTES_PK_HASH);
#endif

	{
		uint8_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];

		expand_T12(T12, skx->sk_seed);
		for (int i1 = 0; i1 < SNOVA_o * SNOVA_v * SNOVA_l2; i1++) {
			skx->T12[i1] = gf16_expand(T12[i1]);
		}

		SNOVA_SECURE_CLEAR(T12);
	}

	uint8_t P_matrix[NUM_PUB_GF];

	expand_public(P_matrix, skx->pk_seed);

	/**
	 * Calculate F12, F21
	 */
	uint8_t *P11 = P_matrix;
	uint8_t *P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	uint8_t *P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2; i1++) {
		skx->P11[i1] = P11[i1];
	}

#if SNOVA_l == 5
	uint16_t F21[SNOVA_mvl * SNOVA_o * SNOVA_l] = {0};

	uint16_t P11a[SNOVA_mvl * SNOVA_v * SNOVA_l] = {0};
	for (int nk = 0; nk < SNOVA_v; nk++)
		for (int k1 = 0; k1 < SNOVA_l; k1++)
			for (int mi = 0; mi < SNOVA_m1; mi++)
				for (int nj = 0; nj < SNOVA_v; nj++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						P11a[(nk * SNOVA_l + k1) * SNOVA_mvl + (mi * SNOVA_v + nj) * SNOVA_l + j1] =
						    skx->P11[((mi * SNOVA_v + nk) * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];

	for (int ni = 0; ni < SNOVA_o; ni++)
		for (int i1 = 0; i1 < SNOVA_l; i1++)
			for (int nk = 0; nk < SNOVA_v; nk++)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int mi = 0; mi < SNOVA_mvl; mi++)
						F21[(ni * SNOVA_l + i1) * SNOVA_mvl + mi] ^=
						    skx->T12[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + k1] *
						    P11a[(nk * SNOVA_l + k1) * SNOVA_mvl + mi];

	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int ni = 0; ni < SNOVA_o; ni++)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						skx->F21[((mi * SNOVA_o + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    gf16_compress(F21[(ni * SNOVA_l + i1) * SNOVA_mvl + (mi * SNOVA_v + nj) * SNOVA_l + j1]) ^
						    P21[((mi * SNOVA_o + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];

	for (int ni = 0; ni < SNOVA_v; ni++)
		for (int k1 = 0; k1 < SNOVA_l; k1++)
			for (int mi = 0; mi < SNOVA_m1; mi++)
				for (int nj = 0; nj < SNOVA_v; nj++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						P11a[(ni * SNOVA_l + k1) * SNOVA_mvl + (mi * SNOVA_v + nj) * SNOVA_l + i1] =
						    skx->P11[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + k1];

	uint16_t F12[SNOVA_mvl * SNOVA_o * SNOVA_l] = {0};
	for (int nk = 0; nk < SNOVA_o; nk++)
		for (int j1 = 0; j1 < SNOVA_l; j1++)
			for (int ni = 0; ni < SNOVA_v; ni++)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int mi = 0; mi < SNOVA_mvl; mi++)
						F12[(nk * SNOVA_l + j1) * SNOVA_mvl + mi] ^=
						    P11a[(ni * SNOVA_l + k1) * SNOVA_mvl + mi] *
						    skx->T12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1];

	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int nk = 0; nk < SNOVA_o; nk++)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						skx->F12[((mi * SNOVA_v + nj) * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    gf16_compress(F12[(nk * SNOVA_l + j1) * SNOVA_mvl + (mi * SNOVA_v + nj) * SNOVA_l + i1]) ^
						    P12[((mi * SNOVA_v + nj) * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1];
#else
	uint16_t F21[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int nk = 0; nk < SNOVA_v; nk++)
				for (int ni = 0; ni < SNOVA_o; ni++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								F21[((mi * SNOVA_o + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    skx->T12[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    skx->P11[((mi * SNOVA_v + nk) * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		skx->F21[i1] = gf16_compress(F21[i1]) ^ P21[i1];
	}

	uint16_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m1; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int ni = 0; ni < SNOVA_v; ni++)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int nk = 0; nk < SNOVA_o; nk++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								F12[((mi * SNOVA_v + nj) * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    skx->P11[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    skx->T12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1];

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		skx->F12[i1] = gf16_compress(F12[i1]) ^ P12[i1];
	}
#endif

#if !FIXED_ABQ
	// Generate ABQ
	uint8_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
	uint8_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
	uint8_t Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
	uint8_t Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2];

	uint8_t *aptr = P_matrix + SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2;
	uint8_t *q1 = aptr + SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr);
	uint8_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	uint8_t *A = P_matrix + (SNOVA_m1 * (SNOVA_n * SNOVA_n - SNOVA_o * SNOVA_o)) * SNOVA_l2;

	gen_ABQ(A, Am, Bm, Q1, Q2);

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_r2; i1++) {
		skx->Am[i1] = gf16_expand(Am[i1]);
	}
	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_lr; i1++) {
		skx->Bm[i1] = gf16_expand(Bm[i1]);
	}
	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_l2; i1++) {
		skx->Q1[i1] = gf16_expand(Q1[i1]);
		skx->Q2[i1] = gf16_expand(Q2[i1]);
	}

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_l; i1++) {
		skx->q1[i1] = gf16_expand(q1[i1]);
		skx->q2[i1] = gf16_expand(q2[i1]);
	}
#endif

	SNOVA_SECURE_CLEAR(F12);
	SNOVA_SECURE_CLEAR(F21);

	return 0;

#endif
}

/**
 * Optimized version of Sign. Deterministic using the salt provided
 */
int SNOVA_NAMESPACE(sign)(const expanded_SK *skx, uint8_t *sig, const uint8_t *digest, const size_t len_digest,
                          const uint8_t *salt) {
	SNOVA_INIT

#if FIXED_ABQ
#define Am fixedAm
#define Bm fixedBm
#define Q1 fixedQ1
#define Q2 fixedQ2
#define q1 fixedq1
#define q2 fixedq2
#else
#define Am skx->Am
#define Bm skx->Bm
#define Q1 skx->Q1
#define Q2 skx->Q2
#define q1 skx->q1
#define q2 skx->q2
#endif

	// Calculate message has of size l^2o
	uint8_t hash_in_GF16[GF16_HASH];

	uint8_t sign_hashb[BYTES_HASH];
#if HASH_PK
	hash_combined(sign_hashb, digest, len_digest, skx->pk_hash, salt);
#else
	hash_combined(sign_hashb, digest, len_digest, skx->pk_seed, salt);
#endif
	convert_bytes_to_GF(hash_in_GF16, sign_hashb, GF16_HASH);

	// Find a solution for T.X
	uint8_t signature_in_GF[SNOVA_n * SNOVA_lr];
	uint8_t num_sign = 0;

	while (1) {
		num_sign++;
		if (num_sign == 255) {
			// Probability of getting here is about q^{-255}
			memset(sig, 0, BYTES_SIGNATURE);
			return -1;
		}

		int flag_redo = 0;

		/**
		 * All variables that may contain secrets
		 */
		uint8_t solution[SNOVA_olr] = {0};
		alignas(32) uint16_t gauss16[SNOVA_olr][SNOVA_olr] = {0};

		uint8_t vinegar_in_byte[NUM_GEN_SEC_BYTES];
		uint16_t Fvv_in_GF16Matrix[SNOVA_olr] = {0};

		alignas(32) uint16_t whipped_sig[SNOVA_l * SNOVA_v * SNOVA_lr32] = {0};
#if USE_GFNI
		alignas(32) uint8_t whipped_sig2[SNOVA_l * SNOVA_v * SNOVA_lr32] = {0};
		alignas(32) uint8_t sum_t0[SNOVA_v * SNOVA_l * SNOVA_lr32];
		alignas(32) uint8_t sum_t1[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32] = {0};
#else
		
#if defined(__aarch64__) && !defined(ROUND3_DISABLE_NEON)
		alignas(32) uint8_t whipped_sig2[SNOVA_l * SNOVA_v * SNOVA_lr32] = {0};
#else
		alignas(32) uint16_t whipped_sig2[SNOVA_l * SNOVA_v * SNOVA_lr32] = {0};
#endif
		alignas(32) uint16_t sum_t0[SNOVA_v * SNOVA_l * SNOVA_lr32];
		alignas(32) uint16_t sum_t1[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32] = {0};
#endif
		alignas(32) uint16_t sum_t1s[SNOVA_m1 * SNOVA_l * SNOVA_l * SNOVA_r2] = {0};
		uint16_t gausstmp1[SNOVA_m1 * SNOVA_r * SNOVA_r * SNOVA_olr] = {0};
		uint16_t gausstmp2[SNOVA_m1 * SNOVA_r * SNOVA_r * SNOVA_olr] = {0};

		uint16_t gfm_temp0[SNOVA_r2];
		uint16_t gfm_temp1[SNOVA_r2];
		uint16_t gfm_temp2[SNOVA_lr];

		alignas(32) uint16_t whipped_F21[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr32] = {0};
		uint16_t whipped_F12[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr32] = {0};
		uint16_t whipped_F12s[SNOVA_m1 * SNOVA_l * SNOVA_olr] = {0};

		alignas(32) uint16_t coef_temp0[SNOVA_olr];
		alignas(32) uint16_t coef_temp1[SNOVA_olr];
		alignas(32) uint16_t coef_temp2[SNOVA_olr];

		uint16_t t_GF16;
		uint16_t fact;
		uint16_t solution16[SNOVA_o * SNOVA_lr] = {0};

		shake_t v_instance;

		/*********************************************************************** */

		// generate the vinegar value
		shake256_init(&v_instance);
		shake_absorb(&v_instance, skx->sk_seed, SEED_LENGTH_PRIVATE);
		shake_absorb(&v_instance, digest, BYTES_DIGEST);
		shake_absorb(&v_instance, salt, BYTES_SALT);
		shake_absorb(&v_instance, &num_sign, 1);
		shake_finalize(&v_instance);
		shake_squeeze(vinegar_in_byte, NUM_GEN_SEC_BYTES, &v_instance);

		convert_bytes_to_GF(signature_in_GF, vinegar_in_byte, SNOVA_v * SNOVA_lr);
		for (int ab = 0; ab < SNOVA_l; ++ab)
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_r; j1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							whipped_sig[(i1 * SNOVA_v + ni) * SNOVA_lr32 + ab * SNOVA_r + j1] ^=
							    gf_Sx[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r + j1];

		for (int i1 = 0; i1 < SNOVA_v * SNOVA_l * SNOVA_lr32; i1++) {
			whipped_sig2[i1] = gf16_compress(whipped_sig[i1]);
			whipped_sig[i1] = gf16_expand(whipped_sig2[i1]);
		}

		/**
		 * Evaluate whipped central map
		 */
#if USE_GFNI
		for (int mi = 0; mi < SNOVA_m1; ++mi) {
			// Right
			memset(sum_t0, 0, sizeof(sum_t0));
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int nj = 0; nj < SNOVA_v; ++nj)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
								__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32];
								__m256i q1_256 = _mm256_set1_epi8(
								                     skx->P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
								__m256i *wip_256 = (__m256i *)&whipped_sig2[(k1 * SNOVA_v + nj) * SNOVA_lr32];

								sum0_256[b1] ^= _mm256_gf2p8mul_epi8(q1_256, wip_256[b1]);
							}

			for (int i1 = 0; i1 < SNOVA_v * SNOVA_l * SNOVA_lr32; ++i1) {
				sum_t0[i1] = gfni_cleanup(sum_t0[i1]);
			}

			// Left, transposed whipped_sig
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
								__m256i *sum1_256 =
								    (__m256i *)&sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
								__m256i wip_256 =
								    _mm256_set1_epi8(whipped_sig2[(k1 * SNOVA_v + ni) * SNOVA_lr32 + a1 * SNOVA_r + i1]);
								__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32];

								sum1_256[b1] ^= _mm256_gf2p8mul_epi8(wip_256, sum0_256[b1]);
							}
		}

		for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; ++i1) {
			sum_t1[i1] = gfni_cleanup(sum_t1[i1]);
		}
#elif USE_AVX2
		// Right
		for (int mi = 0; mi < SNOVA_m1; ++mi) {
			memset(sum_t0, 0, sizeof(sum_t0));
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int nj = 0; nj < SNOVA_v; ++nj)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
								__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32];
								__m256i q1_256 = _mm256_set1_epi16(
								                     skx->P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
								__m256i *wip_256 = (__m256i *)&whipped_sig[(k1 * SNOVA_v + nj) * SNOVA_lr32];

								sum0_256[b1] ^= _mm256_mullo_epi16(q1_256, wip_256[b1]);
							}

			for (int i1 = 0; i1 < SNOVA_v * SNOVA_l * SNOVA_lr32; ++i1) {
				sum_t0[i1] = gf16_compress(sum_t0[i1]);
			}

			// Left, transposed whipped_sig
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
								__m256i *sum1_256 =
								    (__m256i *)&sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
								__m256i wip_256 =
								    _mm256_set1_epi16(whipped_sig[(k1 * SNOVA_v + ni) * SNOVA_lr32 + a1 * SNOVA_r + i1]);
								__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32];

								sum1_256[b1] ^= _mm256_mullo_epi16(wip_256, sum0_256[b1]);
							}
		}

		for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; ++i1) {
			sum_t1[i1] = gf16_compress(sum_t1[i1]);
		}
#elif defined(__aarch64__) && !defined(ROUND3_DISABLE_NEON)
		round3_neon_contract(skx->P11, whipped_sig2, sum_t1, SNOVA_v, 1);
#else
		// Right
		for (int mi = 0; mi < SNOVA_m1; ++mi)
			for (int ni = 0; ni < SNOVA_v; ++ni) {
				memset(sum_t0, 0, sizeof(sum_t0));
				for (int nj = 0; nj < SNOVA_v; ++nj)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int b1 = 0; b1 < SNOVA_lr32; ++b1)
								sum_t0[i1 * SNOVA_lr32 + b1] ^=
								    skx->P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    whipped_sig[(k1 * SNOVA_v + nj) * SNOVA_lr32 + b1];

				for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr32; ++i1) {
					sum_t0[i1] = gf16_compress(sum_t0[i1]);
				}

				// Left, transposed whipped_sig
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int i1 = 0; i1 < SNOVA_r; i1++)
							for (int b1 = 0; b1 < SNOVA_lr32; ++b1)
								sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32 + b1] ^=
								    whipped_sig[(k1 * SNOVA_v + ni) * SNOVA_lr32 + a1 * SNOVA_r + i1] *
								    sum_t0[k1 * SNOVA_lr32 + b1];
			}

		for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; i1++) {
			sum_t1[i1] = gf16_compress(sum_t1[i1]);
		}
#endif

		for (int mi = 0; mi < SNOVA_m1; ++mi)
			for (int a1 = 0; a1 < SNOVA_l; ++a1)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int j1 = 0; j1 < SNOVA_r; j1++)
							sum_t1s[(mi * SNOVA_l + a1) * SNOVA_l * SNOVA_r2 + b1 * SNOVA_r2 + i1 * SNOVA_r + j1] =
							    sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32 + b1 * SNOVA_r + j1];

		/**
		 * Apply A, B, q1 and q2, aka E matrix
		 */
		for (int mi = 0; mi < SNOVA_o; ++mi)
			for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
				int mi_prime = i_prime(mi, alpha);

				memset(gfm_temp1, 0, sizeof(gfm_temp1));
				memset(gfm_temp2, 0, sizeof(gfm_temp2));

				// apply q1 and q2
				for (int a1 = 0; a1 < SNOVA_l; ++a1) {
					memset(gfm_temp0, 0, sizeof(gfm_temp0));

					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_r2; i1++)
							gfm_temp0[i1] ^= sum_t1s[(mi_prime * SNOVA_l + a1) * SNOVA_l * SNOVA_r2 + b1 * SNOVA_r2 + i1] *
							                 q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

					for (int i1 = 0; i1 < SNOVA_r2; i1++) {
						gfm_temp0[i1] = gf16_compress(gfm_temp0[i1]);
					}

					for (int i1 = 0; i1 < SNOVA_r2; i1++) {
						gfm_temp1[i1] ^= gfm_temp0[i1] * q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
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
							    gfm_temp1[i1 * SNOVA_r + k1] * Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];
						}

				for (int i1 = 0; i1 < SNOVA_lr; i1++) {
					gfm_temp2[i1] = gf16_compress(gfm_temp2[i1]);
				}

				for (int i1 = 0; i1 < SNOVA_r; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int k1 = 0; k1 < SNOVA_r; k1++) {
							Fvv_in_GF16Matrix[mi * SNOVA_lr + i1 * SNOVA_l + j1] ^=
							    Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] * gfm_temp2[k1 * SNOVA_l + j1];
						}

				// Set the last column of gauss matrix
				for (int i1 = 0; i1 < SNOVA_lr; i1++) {
					gauss16[mi * SNOVA_lr + i1][SNOVA_o * SNOVA_lr] =
					    hash_in_GF16[mi * SNOVA_lr + i1] ^ Fvv_in_GF16Matrix[mi * SNOVA_lr + i1];
				}
			}

		// Whipped F21
		for (int mi = 0; mi < SNOVA_m1; mi++)
			for (int idx = 0; idx < SNOVA_o; idx++)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int nj = 0; nj < SNOVA_v; ++nj)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								for (int k1 = 0; k1 < SNOVA_l; k1++)
									whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1] ^=
									    skx->F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1] *
									    whipped_sig[(k1 * SNOVA_v + nj) * SNOVA_lr32 + b1 * SNOVA_r + j1];

		for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr32; i1++) {
			whipped_F21[i1] = gf16_compress(whipped_F21[i1]);
		}

		// Whipped F12
		for (int mi = 0; mi < SNOVA_m1; mi++)
			for (int nj = 0; nj < SNOVA_v; ++nj)
				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int b1 = 0; b1 < SNOVA_lr32; ++b1)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								whipped_F12[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32 + b1] ^=
								    skx->F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1] *
								    whipped_sig[(k1 * SNOVA_v + nj) * SNOVA_lr32 + b1];

		for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr32; i1++) {
			whipped_F12[i1] = gf16_compress(whipped_F12[i1]);
		}

		for (int mi = 0; mi < SNOVA_m1; mi++)
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int j1 = 0; j1 < SNOVA_r; j1++)
							whipped_F12s[(mi * SNOVA_l + b1) * SNOVA_olr + idx * SNOVA_lr + i1 * SNOVA_r + j1] =
							    whipped_F12[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32 + b1 * SNOVA_r + j1];

		// compute the coefficients of Xo and put into gauss matrix and compute
		// the coefficients of Xo^t and add into gauss matrix
		for (int mi = 0; mi < SNOVA_o; mi++)
			for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
				memset(coef_temp0, 0, sizeof(coef_temp0));
				memset(coef_temp1, 0, sizeof(coef_temp1));
				memset(coef_temp2, 0, sizeof(coef_temp2));

				int mi_prime = i_prime(mi, alpha);

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_r; j1++)
								coef_temp0[idx * SNOVA_lr + i1 * SNOVA_r + j1] ^=
								    whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1] *
								    q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
					coef_temp0[i1] = gf16_compress(coef_temp0[i1]);
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_r; k1++) {
								coef_temp1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    coef_temp0[idx * SNOVA_lr + i1 * SNOVA_r + k1] *
								    Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];
							}

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					coef_temp1[i1] = gf16_compress(coef_temp1[i1]);
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								coef_temp2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
								    Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    coef_temp1[idx * SNOVA_l2 + k1 * SNOVA_l + j1];

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					coef_temp2[i1] = gf16_compress(coef_temp2[i1]);
				}

				for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
					for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
						for (int idx = 0; idx < SNOVA_olr; idx++)
							gausstmp1[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * SNOVA_olr + idx] ^=
							    coef_temp2[idx] * Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r + tj2];
			}

		for (int mi = 0; mi < SNOVA_o; mi++)
			for (int idx = 0; idx < SNOVA_o; idx++)
				for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
					for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
						for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
							for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
								gauss16[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr + tj1 * SNOVA_r + tj2] =
								    gausstmp1[(mi * SNOVA_r * SNOVA_r + ti1 * SNOVA_r + tj2) * SNOVA_olr + idx * SNOVA_l2 +
								              tj1 * SNOVA_l + ti2];

		for (int mi = 0; mi < SNOVA_o; mi++) {
			for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
				memset(coef_temp0, 0, sizeof(coef_temp0));
				memset(coef_temp1, 0, sizeof(coef_temp1));
				memset(coef_temp2, 0, sizeof(coef_temp2));

				int mi_prime = i_prime(mi, alpha);

				// Transpose
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int idx = 0; idx < SNOVA_olr; idx++)
						coef_temp0[idx] ^= whipped_F12s[(mi_prime * SNOVA_l + b1) * SNOVA_olr + idx] *
						                   q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

				for (int i1 = 0; i1 < SNOVA_olr; i1++) {
					coef_temp0[i1] = gf16_compress(coef_temp0[i1]);
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_r; k1++)
								coef_temp1[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
								    Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] *
								    coef_temp0[idx * SNOVA_lr + j1 * SNOVA_r + k1];

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
					coef_temp1[i1] = gf16_compress(coef_temp1[i1]);
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								coef_temp2[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
								    coef_temp1[idx * SNOVA_lr + i1 * SNOVA_l + k1] *
								    Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + k1 * SNOVA_l + j1];

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
					coef_temp2[i1] = gf16_compress(coef_temp2[i1]);
				}

				for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
					for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
						for (int idx = 0; idx < SNOVA_olr; idx++)
							gausstmp2[(mi * SNOVA_r * SNOVA_r + ti2 * SNOVA_r + tj2) * SNOVA_olr + idx] ^=
							    coef_temp2[idx] * Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2];
			}
		}

		for (int mi = 0; mi < SNOVA_o; mi++)
			for (int idx = 0; idx < SNOVA_o; idx++)
				for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
					for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
						for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
							for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
								gauss16[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr + tj1 * SNOVA_r + tj2] ^=
								    gausstmp2[(mi * SNOVA_r * SNOVA_r + ti2 * SNOVA_r + tj2) * SNOVA_olr + idx * SNOVA_lr +
								              ti1 * SNOVA_l + tj1];

		for (int ti = 0; ti < SNOVA_o * SNOVA_lr; ti++)
			for (int tj = 0; tj < SNOVA_o * SNOVA_lr; tj++) {
				gauss16[ti][tj] = gf16_compress(gauss16[ti][tj]);
			}

		// Gaussian elimination in constant time
		for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) {
			for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
				int16_t mask = ct_is_not_zero(gauss16[i][i]) - 1;
				for (int k = 0; k < SNOVA_olr; ++k) {
					gauss16[i][k] ^= mask & gauss16[j][k];
				}
			}

			for (int k = 0; k < SNOVA_olr; ++k) {
				gauss16[i][k] = gf16_compress(gauss16[i][k]);
			}

			flag_redo |= 1 - ct_is_not_zero(gauss16[i][i]);

			/**
			 * Constant time GF(16) inverse
			 *
			 * Use that x^q = x and therefore x^(q-2) = x^-1
			 *
			 * t_GF16 = ct_gf_inverse(gauss16[i][i]);
			 */
			fact = gauss16[i][i] * gf16_expand(gauss16[i][i]);
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
				gauss16[i][k] = gauss16[i][k] * t_GF16;
			}

			for (int k = 0; k < SNOVA_olr; ++k) {
				gauss16[i][k] = gf16_compress(gauss16[i][k]);
			}

			for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
				t_GF16 = gf16_expand(gauss16[j][i]);
				for (int k = 0; k < SNOVA_olr; ++k) {
					gauss16[j][k] ^= gauss16[i][k] * t_GF16;
				}
			}

			for (int j = i + 1; j < SNOVA_o * SNOVA_lr; ++j) {
				gauss16[j][i + 1] = gf16_compress(gauss16[j][i + 1]);
			}
		}

		if (flag_redo) {
			continue;
		}

		// Last step of Gaussian elimination
		for (int i = SNOVA_o * SNOVA_lr - 1; i >= 0; --i) {
			t_GF16 = 0;
			for (int k = i + 1; k < SNOVA_o * SNOVA_lr; ++k) {
				t_GF16 ^= gauss16[i][k] * solution16[k];
			}
			t_GF16 ^= gauss16[i][SNOVA_o * SNOVA_lr];
			t_GF16 = gf16_compress(t_GF16);
			solution16[i] = gf16_expand(t_GF16);
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
						t_GF16 = 0;
						for (int k1 = 0; k1 < SNOVA_l; k1++) {
							t_GF16 ^= skx->T12[(index * SNOVA_o + mi) * SNOVA_l2 + i1 * SNOVA_l + k1] *
							          solution[mi * SNOVA_lr + k1 * SNOVA_r + j1];
						}
						signature_in_GF[index * SNOVA_lr + i1 * SNOVA_r + j1] ^= gf16_compress(t_GF16);
					}

		memcpy(signature_in_GF + SNOVA_v * SNOVA_lr, solution, SNOVA_o * SNOVA_lr);

		SNOVA_SECURE_CLEAR(solution);
		SNOVA_SECURE_CLEAR(gauss16);
		SNOVA_SECURE_CLEAR(vinegar_in_byte);
		SNOVA_SECURE_CLEAR(Fvv_in_GF16Matrix);
		SNOVA_SECURE_CLEAR(whipped_sig);
		SNOVA_SECURE_CLEAR(whipped_sig2);
		SNOVA_SECURE_CLEAR(sum_t0);
		SNOVA_SECURE_CLEAR(sum_t1);
		SNOVA_SECURE_CLEAR(sum_t1s);
		SNOVA_SECURE_CLEAR(gausstmp1);
		SNOVA_SECURE_CLEAR(gausstmp2);
		SNOVA_SECURE_CLEAR(gfm_temp0);
		SNOVA_SECURE_CLEAR(gfm_temp1);
		SNOVA_SECURE_CLEAR(gfm_temp2);
		SNOVA_SECURE_CLEAR(whipped_F21);
		SNOVA_SECURE_CLEAR(whipped_F12);
		SNOVA_SECURE_CLEAR(whipped_F12s);
		SNOVA_SECURE_CLEAR(coef_temp0);
		SNOVA_SECURE_CLEAR(coef_temp1);
		SNOVA_SECURE_CLEAR(coef_temp2);
		SNOVA_SECURE_CLEARI(t_GF16);
		SNOVA_SECURE_CLEARI(fact);
		SNOVA_SECURE_CLEAR(solution16);
		SNOVA_SECURE_CLEARI(v_instance);

		break;
	}

	compress_gf(sig, signature_in_GF, SNOVA_n * SNOVA_lr);
	memcpy(sig + BYTES_SIGNATURE - BYTES_SALT, salt, BYTES_SALT);

	return 0;
}

#undef Am
#undef Bm
#undef Q1
#undef Q2
#undef q1
#undef q2

/**
 * PK expansion.
 */
int SNOVA_NAMESPACE(pk_expand)(expanded_PK *pkx, const uint8_t *pk) {
	SNOVA_INIT

	memcpy(pkx->pk_seed, pk, SEED_LENGTH_PUBLIC);
#if HASH_PK
	shake256(pkx->pk_hash, BYTES_PK_HASH, pk, BYTES_PK);
#endif

	/**
	 * Create P matrix
	 */
	uint8_t P_matrix[NUM_PUB_GF];
	uint8_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];

	uint8_t *P11 = P_matrix;
	uint8_t *P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
	uint8_t *P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;

	convert_bytes_to_GF(P22, pk + SEED_LENGTH_PUBLIC, NUMGF_PK);
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

#if !FIXED_ABQ
	/**
	 * Create AB matrices, improve q
	 */
	uint8_t *A = P_matrix + SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2;
	uint8_t *B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
	uint8_t *q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
	uint8_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (size_t idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&(pkx->Am[idx * SNOVA_r2]), &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
		be_invertible_by_add_aS(&(pkx->Bm[idx * SNOVA_lr]), &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);

		if (!q1[idx * SNOVA_l + SNOVA_l - 1]) {
			q1[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q1[idx * SNOVA_l] + (q1[idx * SNOVA_l] == 0));
		}
		if (!q2[idx * SNOVA_l + SNOVA_l - 1]) {
			q2[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q2[idx * SNOVA_l] + (q2[idx * SNOVA_l] == 0));
		}
	}

	memcpy(pkx->q1, q1, SNOVA_o * SNOVA_alpha * SNOVA_l);
	memcpy(pkx->q2, q2, SNOVA_o * SNOVA_alpha * SNOVA_l);
#endif

	return 0;
}

/**
 * Optimized version of verify.
 */
int SNOVA_NAMESPACE(verify)(const expanded_PK *pkx, const uint8_t *sig, const uint8_t *digest, const size_t len_digest) {
	SNOVA_INIT

	uint8_t signature_in_GF[NUMGF_SIGNATURE];
	convert_bytes_to_GF(signature_in_GF, sig, NUMGF_SIGNATURE);

	/**
	 * Whip signature
	 */
	alignas(32) uint16_t whipped_sig[SNOVA_l * SNOVA_n * SNOVA_lr32] = {0};
	alignas(32) uint8_t whipped_sig2[SNOVA_l * SNOVA_n * SNOVA_lr32] = {0};

	for (int ab = 0; ab < SNOVA_l; ++ab)
		for (int idx = 0; idx < SNOVA_n; ++idx)
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = 0; j1 < SNOVA_r; j1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						whipped_sig[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32 + ab * SNOVA_r + j1] ^=
						    gf_Sx[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r + j1];

	for (int i1 = 0; i1 < SNOVA_l * SNOVA_n * SNOVA_lr32; i1++) {
		whipped_sig2[i1] = gf16_compress(whipped_sig[i1]);
		whipped_sig[i1] = gf16_expand(whipped_sig2[i1]);
	}

	/**
	 * Evaluate whipped central map
	 */
	uint16_t hash_in_GF[SNOVA_o * SNOVA_lr] = {0};

#if USE_GFNI
	alignas(32) uint8_t sum_t1[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32] = {0};
	for (int mi = 0; mi < SNOVA_m1; ++mi) {
		// Right
		alignas(32) uint8_t sum_t0[SNOVA_n * SNOVA_l * SNOVA_lr32] = {0};
		for (int ni = 0; ni < SNOVA_n; ++ni)
			for (int nj = 0; nj < SNOVA_n; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
							__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32];
							__m256i q1_256 =
							    _mm256_set1_epi8(pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
							__m256i *wip_256 = (__m256i *)&whipped_sig2[nj * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32];

							sum0_256[b1] ^= _mm256_gf2p8mul_epi8(q1_256, wip_256[b1]);
						}

		for (int i1 = 0; i1 < SNOVA_n * SNOVA_l * SNOVA_lr32; ++i1) {
			sum_t0[i1] = gfni_cleanup(sum_t0[i1]);
		}

		// Left, transposed whipped_sig
		for (int ni = 0; ni < SNOVA_n; ++ni)
			for (int a1 = 0; a1 < SNOVA_l; ++a1)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
							__m256i *sum1_256 = (__m256i *)&sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
							__m256i wip_256 =
							    _mm256_set1_epi8(whipped_sig2[ni * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32 + a1 * SNOVA_r + i1]);
							__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32];

							sum1_256[b1] ^= _mm256_gf2p8mul_epi8(wip_256, sum0_256[b1]);
						}
	}

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; ++i1) {
		sum_t1[i1] = gfni_cleanup(sum_t1[i1]);
	}
#elif USE_AVX2_VERIFY
	alignas(32) uint16_t sum_t1[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32] = {0};
	for (int mi = 0; mi < SNOVA_m1; ++mi) {
		// Right
		alignas(32) uint16_t sum_t0[SNOVA_n * SNOVA_l * SNOVA_lr32] = {0};
		for (int ni = 0; ni < SNOVA_n; ++ni)
			for (int nj = 0; nj < SNOVA_n; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
							__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32];
							__m256i q1_256 =
							    _mm256_set1_epi16(pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
							__m256i *wip_256 = (__m256i *)&whipped_sig[nj * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32];

							sum0_256[b1] ^= _mm256_mullo_epi16(q1_256, wip_256[b1]);
						}

		for (int i1 = 0; i1 < SNOVA_n * SNOVA_l * SNOVA_lr32; ++i1) {
			sum_t0[i1] = gf16_compress(sum_t0[i1]);
		}

		// Left, transposed whipped_sig
		for (int ni = 0; ni < SNOVA_n; ++ni)
			for (int a1 = 0; a1 < SNOVA_l; ++a1)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
							__m256i *sum1_256 = (__m256i *)&sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
							__m256i wip_256 =
							    _mm256_set1_epi16(whipped_sig[ni * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32 + a1 * SNOVA_r + i1]);
							__m256i *sum0_256 = (__m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32];

							sum1_256[b1] ^= _mm256_mullo_epi16(wip_256, sum0_256[b1]);
						}
	}

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; ++i1) {
		sum_t1[i1] = gf16_compress(sum_t1[i1]);
	}
#elif defined(__aarch64__) && !defined(ROUND3_DISABLE_NEON)
	alignas(32) uint16_t sum_t1[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32] = {0};
	round3_neon_contract(pkx->P, whipped_sig2, sum_t1, SNOVA_n, 0);
#else
	alignas(32) uint16_t sum_t1[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32] = {0};
	for (int mi = 0; mi < SNOVA_m1; ++mi) {
		// Right
		for (int ni = 0; ni < SNOVA_n; ++ni) {
			alignas(32) uint16_t sum_t0[SNOVA_n * SNOVA_l * SNOVA_lr32] = {0};
			for (int nj = 0; nj < SNOVA_n; ++nj)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int b1 = 0; b1 < SNOVA_lr32; ++b1)
							sum_t0[i1 * SNOVA_lr32 + b1] ^=
							    pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + k1] *
							    whipped_sig[nj * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32 + b1];

			for (int i1 = 0; i1 < SNOVA_l * SNOVA_lr32; ++i1) {
				sum_t0[i1] = gf16_compress(sum_t0[i1]);
			}

			// Left, transposed whipped_sig
			for (int a1 = 0; a1 < SNOVA_l; ++a1)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int i1 = 0; i1 < SNOVA_r; i1++)
						for (int b1 = 0; b1 < SNOVA_lr32; ++b1)
							sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32 + b1] ^=
							    whipped_sig[ni * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32 + a1 * SNOVA_r + i1] *
							    sum_t0[k1 * SNOVA_lr32 + b1];
		}
	}

	for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; i1++) {
		sum_t1[i1] = gf16_compress(sum_t1[i1]);
	}
#endif

	alignas(32) uint16_t sum_t1s[SNOVA_m1 * SNOVA_l * SNOVA_l * SNOVA_r2] = {0};

	for (int mi = 0; mi < SNOVA_m1; ++mi)
		for (int a1 = 0; a1 < SNOVA_l; ++a1)
			for (int i1 = 0; i1 < SNOVA_r; i1++)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int j1 = 0; j1 < SNOVA_r; j1++)
						sum_t1s[((mi * SNOVA_l + a1) * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1] =
						    sum_t1[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32 + b1 * SNOVA_r + j1];

	/**
	 * Prepare
	 */
#if FIXED_ABQ
#define Amx fixedAm
#define Bmx fixedBm
#define q1x fixedq1
#define q2x fixedq2
#else

	uint16_t Amx[SNOVA_o * SNOVA_alpha * SNOVA_r2];
	uint16_t Bmx[SNOVA_o * SNOVA_alpha * SNOVA_lr];
	uint16_t q1x[SNOVA_o * SNOVA_alpha * SNOVA_l];
	uint16_t q2x[SNOVA_o * SNOVA_alpha * SNOVA_l];

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_r2; i1++) {
		Amx[i1] = gf16_expand(pkx->Am[i1]);
	}

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_lr; i1++) {
		Bmx[i1] = gf16_expand(pkx->Bm[i1]);
	}

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_l; i1++) {
		q1x[i1] = gf16_expand(pkx->q1[i1]);
	}

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_alpha * SNOVA_l; i1++) {
		q2x[i1] = gf16_expand(pkx->q2[i1]);
	}
#endif

	/**
	 * Apply A, B, q1 and q2, aka E matrix
	 */
	for (int mi = 0; mi < SNOVA_o; ++mi)
		for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
			int mi_prime = i_prime(mi, alpha);

			uint16_t gfm_temp1[SNOVA_r2] = {0};
			uint16_t gfm_temp2[SNOVA_lr] = {0};

			// apply q1 and q2
			for (int a1 = 0; a1 < SNOVA_l; ++a1) {
				uint16_t gfm_temp0[SNOVA_r2] = {0};

				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int i1 = 0; i1 < SNOVA_r2; i1++)
						gfm_temp0[i1] ^= q2x[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1] *
						                 sum_t1s[(mi_prime * SNOVA_l + a1) * SNOVA_l * SNOVA_r2 + b1 * SNOVA_r2 + i1];

				for (int i1 = 0; i1 < SNOVA_r2; i1++) {
					gfm_temp0[i1] = gf16_compress(gfm_temp0[i1]);
				}

				for (int i1 = 0; i1 < SNOVA_r2; i1++) {
					gfm_temp1[i1] ^= gfm_temp0[i1] * q1x[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
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
						    gfm_temp1[i1 * SNOVA_r + k1] * Bmx[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];

			for (int i1 = 0; i1 < SNOVA_lr; i1++) {
				gfm_temp2[i1] = gf16_compress(gfm_temp2[i1]);
			}

			for (int i1 = 0; i1 < SNOVA_r; i1++)
				for (int j1 = 0; j1 < SNOVA_l; j1++)
					for (int k1 = 0; k1 < SNOVA_r; k1++)
						hash_in_GF[mi * SNOVA_lr + i1 * SNOVA_l + j1] ^=
						    Amx[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] * gfm_temp2[k1 * SNOVA_l + j1];
		}

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
		hash_in_GF[i1] = gf16_compress(hash_in_GF[i1]);
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
	convert_bytes_to_GF(signed_gf, signed_bytes, GF16_HASH);

	for (int i = 0; i < GF16_HASH; ++i) {
		if (hash_in_GF[i] != signed_gf[i]) {
			return -1;
		}
	}

	return 0;
}
