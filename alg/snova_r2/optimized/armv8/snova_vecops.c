#include "snova_vecops.h"

alignas(32) uint8_t mt4b[256] = {0};
alignas(32) uint8_t inv4b[16] = {0};

alignas(16 * VEC_LENGTH) static uint8_t mt4b2_16[256][VEC_LENGTH] = {0};
VECTOR *mtk2_16 = (VECTOR *)mt4b2_16;

VECTOR vector_inv_table = {0};
VECTOR l_mask = {0};
VECTOR vtl_multmask1;
VECTOR vtl_multmask2;
VECTOR vtl_multmask4;
VECTOR vtl_multmask8;
VECTOR vtl_mult_table1;
VECTOR vtl_mult_table2;
VECTOR vtl_mult_table4;
VECTOR vtl_mult_table8;
VECTOR zero256 = {0};

gf16m_t snova_S[l_SNOVA] = {0};
uint32_t snova_xS[l_SNOVA][lsq_SNOVA] = {0};

void snova_init_gf16_tables(void) {
	static int gf16_tables_is_init = 0;
	if (gf16_tables_is_init) {
		return;
	}
	gf16_tables_is_init = 1;

	uint8_t F_star[15] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9};
	for (int i = 0; i < 16; ++i) {
		mt(0, i) = mt(i, 0) = 0;
	}

	for (int i = 0; i < 15; ++i) {
		for (int j = 0; j < 15; ++j) {
			mt(F_star[i], F_star[j]) = F_star[(i + j) % 15];
		}
	}

	{
		int g = F_star[1], g_inv = F_star[14], gn = 1, gn_inv = 1;
		inv4b[0] = 0;
		inv4b[1] = 1;
		for (int index = 0; index < 14; index++) {
			inv4b[(gn = mt(gn, g))] = (gn_inv = mt(gn_inv, g_inv));
		}
	}
}

void snova_init_vector_table(void) {
	static int vector_table_init_flag = 0;
	if (vector_table_init_flag) {
		return;
	}
	vector_table_init_flag = 1;

	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 16; ++j) {
			for (int k = 0; k < 16; ++k) {
				uint8_t temp = (mt(i, k) << 4) ^ mt(j, k);
				mt4b2_16[i * 16 + j][k] = temp;
#if VEC_LENGTH > 16
				mt4b2_16[i * 16 + j][k + 16] = temp;
#endif
			}
		}
	}

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	uint8_t inv_table[16] = {0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8};
	vector_inv_table = vld1q_u8(inv_table);

	uint8_t numf = 0x0f;
	l_mask = vld1q_dup_u8(&numf);

	uint8_t num1 = 1;
	uint8_t num2 = 2;
	uint8_t num4 = 4;
	uint8_t num8 = 8;
	vtl_multmask1 = vld1q_dup_u8(&num1);
	vtl_multmask2 = vld1q_dup_u8(&num2);
	vtl_multmask4 = vld1q_dup_u8(&num4);
	vtl_multmask8 = vld1q_dup_u8(&num8);
#else
#if VEC_LENGTH > 16
	vector_inv_table = _mm256_setr_epi8(0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8, 0, 1, 9, 14, 13, 11, 7, 6, 15, 2,
	                                    12, 5, 10, 4, 3, 8);
#else
	vector_inv_table = _mm_setr_epi8(0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8);
#endif
	l_mask = _mm256_set1_epi8(0x0f);

	vtl_multmask1 = _mm256_set1_epi8(1);
	vtl_multmask2 = _mm256_set1_epi8(2);
	vtl_multmask4 = _mm256_set1_epi8(4);
	vtl_multmask8 = _mm256_set1_epi8(8);
#endif

	vtl_mult_table1 = mtk2_16[1];
	vtl_mult_table2 = mtk2_16[2];
	vtl_mult_table4 = mtk2_16[4];
	vtl_mult_table8 = mtk2_16[8];
}

void snova_gen_S_array(void) {
	static int S_is_init = 0;
	if (S_is_init) {
		return;
	}
	S_is_init = 1;

	be_aI(snova_S[0], 1);
	be_the_S(snova_S[1]);
	for (int index = 2; index < l_SNOVA; index++) {
		gf16m_mul(snova_S[index - 1], snova_S[1], snova_S[index]);
	}

	for (int index = 0; index < l_SNOVA; index++) {
		for (int ij = 0; ij < lsq_SNOVA; ij++) {
			snova_xS[index][ij] = gf16_from_nibble(snova_S[index][ij]);
		}
	}
}

#define POD(entry, a, b, c, d, e, f, g, h, i, j)                                                                    \
    gf16_get_mul(get_gf16m(entry, a, b), gf16_get_add(gf16_get_mul(get_gf16m(entry, c, d), get_gf16m(entry, e, f)), \
                                                      gf16_get_mul(get_gf16m(entry, g, h), get_gf16m(entry, i, j))))

#if rank == 5
static gf16_t gf16m_det3(const gf16m_t entry, int i0, int i1, int i2, int j0, int j1, int j2) {
	return gf16_get_add(
	           gf16_get_add(gf16_get_mul(get_gf16m(entry, j0, i0),
	                                     gf16_get_add(gf16_get_mul(get_gf16m(entry, j1, i1), get_gf16m(entry, j2, i2)),
	                                             gf16_get_mul(get_gf16m(entry, j1, i2), get_gf16m(entry, j2, i1)))),
	                        gf16_get_mul(get_gf16m(entry, j0, i1),
	                                     gf16_get_add(gf16_get_mul(get_gf16m(entry, j1, i0), get_gf16m(entry, j2, i2)),
	                                             gf16_get_mul(get_gf16m(entry, j1, i2), get_gf16m(entry, j2, i0))))),
	           gf16_get_mul(get_gf16m(entry, j0, i2), gf16_get_add(gf16_get_mul(get_gf16m(entry, j1, i0), get_gf16m(entry, j2, i1)),
	                        gf16_get_mul(get_gf16m(entry, j1, i1), get_gf16m(entry, j2, i0)))));
}

static gf16_t gf16m_det5(const gf16m_t entry) {
	uint8_t acc = 0;
	static const int trip[10][3] = {
		{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {0, 2, 3}, {0, 2, 4},
		{0, 3, 4}, {1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}
	};

	for (int t = 0; t < 10; ++t) {

		int r0 = -1, r1 = -1;
		for (int r = 0; r < 5; ++r) {
			if (r != trip[t][0] && r != trip[t][1] && r != trip[t][2]) {
				if (r0 < 0) {
					r0 = r;
				} else {
					r1 = r;
				}
			}
		}

		uint8_t a = gf16m_det3(entry, trip[t][0], trip[t][1], trip[t][2], 0, 1, 2);
		uint8_t b = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, r0), get_gf16m(entry, 4, r1)),
		                         gf16_get_mul(get_gf16m(entry, 3, r1), get_gf16m(entry, 4, r0)));
		acc ^= gf16_get_mul(a, b);
	}

	return acc;
}
#endif

gf16_t snova_gf16m_det(const gf16m_t entry) {
#if rank == 2
	return gf16_get_add(gf16_get_mul(get_gf16m(entry, 0, 0), get_gf16m(entry, 1, 1)),
	                    gf16_get_mul(get_gf16m(entry, 0, 1), get_gf16m(entry, 1, 0)));
#elif rank == 3
	return gf16_get_add(
	           gf16_get_add(
	               gf16_get_mul(get_gf16m(entry, 0, 0), gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 2)),
	                            gf16_get_mul(get_gf16m(entry, 1, 2), get_gf16m(entry, 2, 1)))),
	               gf16_get_mul(get_gf16m(entry, 0, 1), gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2, 2)),
	                            gf16_get_mul(get_gf16m(entry, 1, 2), get_gf16m(entry, 2, 0))))),
	           gf16_get_mul(get_gf16m(entry, 0, 2), gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2, 1)),
	                        gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 0)))));
#elif rank == 4
	gf16_t d0 = gf16_get_mul(get_gf16m(entry, 0, 0), gf16_get_add(gf16_get_add(POD(entry, 1, 1, 2, 2, 3, 3, 2, 3, 3, 2),
	                         POD(entry, 1, 2, 2, 1, 3, 3, 2, 3, 3, 1)),
	                         POD(entry, 1, 3, 2, 1, 3, 2, 2, 2, 3, 1)));

	gf16_t d1 = gf16_get_mul(get_gf16m(entry, 0, 1), gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 2, 3, 3, 2, 3, 3, 2),
	                         POD(entry, 1, 2, 2, 0, 3, 3, 2, 3, 3, 0)),
	                         POD(entry, 1, 3, 2, 0, 3, 2, 2, 2, 3, 0)));

	gf16_t d2 = gf16_get_mul(get_gf16m(entry, 0, 2), gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 1, 3, 3, 2, 3, 3, 1),
	                         POD(entry, 1, 1, 2, 0, 3, 3, 2, 3, 3, 0)),
	                         POD(entry, 1, 3, 2, 0, 3, 1, 2, 1, 3, 0)));

	gf16_t d3 = gf16_get_mul(get_gf16m(entry, 0, 3), gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 1, 3, 2, 2, 2, 3, 1),
	                         POD(entry, 1, 1, 2, 0, 3, 2, 2, 2, 3, 0)),
	                         POD(entry, 1, 2, 2, 0, 3, 1, 2, 1, 3, 0)));

	return gf16_get_add(gf16_get_add(gf16_get_add(d0, d1), d2), d3);
#elif rank == 5
	return gf16m_det5(entry);
#else
#error "rank > 5 is not supported"
#endif
}

void snova_be_invertible_by_add_aS(gf16m_t source) {
	gf16m_t temp;

	if (snova_gf16m_det(source) == 0) {
		for (uint8_t a = 1; a < 16; ++a) {
			be_aS(temp, a);
			gf16m_add(temp, source, source);
			if (snova_gf16m_det(source) != 0) {
				return;
			}
		}
	}
}

void snova_gen_a_FqS(const gf16_t* c, gf16m_t pt_matrix) {
	gf16m_t temp;

	be_aI(pt_matrix, c[0]);
	for (int i = 1; i < rank - 1; ++i) {
		gf16m_scale(snova_S[i], c[i], temp);
		gf16m_add(pt_matrix, temp, pt_matrix);
	}
	gf16m_scale(snova_S[rank - 1], (c[rank - 1] != 0) ? c[rank - 1] : 16 - (c[0] + (c[0] == 0)), temp);
	gf16m_add(pt_matrix, temp, pt_matrix);
}

void snova_gen_a_FqS_ct(const gf16_t* c, gf16m_t pt_matrix) {
	uint32_t xTemp[lsq_SNOVA] = {0};
	uint32_t cX = gf16_from_nibble(c[0]);

	for (int ij = 0; ij < l_SNOVA; ij++) {
		xTemp[ij * l_SNOVA + ij] = cX;
	}

	for (int i1 = 1; i1 < l_SNOVA - 1; i1++) {
		cX = gf16_from_nibble(c[i1]);
		for (int ij = 0; ij < lsq_SNOVA; ij++) {
			xTemp[ij] ^= cX * snova_xS[i1][ij];
		}
	}

	uint8_t zero = ct_gf16_is_not_zero(c[rank - 1]);
	uint8_t val = zero * c[rank - 1] + (1 - zero) * (15 + ct_gf16_is_not_zero(c[0]) - c[0]);

	cX = gf16_from_nibble(val);
	for (int ij = 0; ij < lsq_SNOVA; ij++) {
		xTemp[ij] ^= cX * snova_xS[l_SNOVA - 1][ij];
	}

	for (int ij = 0; ij < lsq_SNOVA; ij++) {
		pt_matrix[ij] = gf16_to_nibble(xTemp[ij]);
	}
}
