#include "snova_public_layout.h"

void snova_pk_pack(uint8_t* pk, const snova_key_elems* key_elems) {
	memcpy(pk, &key_elems->pk, bytes_pk);
}

void snova_sk_pack(uint8_t* esk, const snova_key_elems* key_elems, const uint8_t* pt_private_key_seed) {
	const uint8_t* sk_gf16_ptr = (const uint8_t*)(key_elems->map1.Aalpha);

	snova_gf16_to_bytes_merge(esk, sk_gf16_ptr, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
	memcpy(esk + (bytes_sk - (seed_length_public + seed_length_private)), key_elems->pk.pt_public_key_seed,
	       seed_length_public);
	memcpy(esk + (bytes_sk - seed_length_private), pt_private_key_seed, seed_length_private);
}

void snova_sk_unpack(sk_gf16* skupk, const uint8_t* esk) {
	snova_bytes_to_gf16_split(esk, (uint8_t*)skupk, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
	memcpy(skupk->pt_public_key_seed, esk + (bytes_sk - (seed_length_public + seed_length_private)),
	       seed_length_public + seed_length_private);
}

void snova_input_P22(uint8_t* P22_gf16s, const uint8_t* P22_bytes) {
	snova_bytes_to_gf16(P22_bytes, P22_gf16s, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

#define PUBLIC_ROW_STRIDE vtl_mainRow_x_rank(n_SNOVA)

static inline void clear_public_layout_padding(uint8_t* pj) {
	const int logical_cols = n_SNOVA * rank;
	const int padding = PUBLIC_ROW_STRIDE - logical_cols;

	if (padding == 0) {
		return;
	}

	for (int row = 0; row < logical_cols; ++row) {
		memset(pj + row * PUBLIC_ROW_STRIDE + logical_cols, 0, (size_t)padding);
	}
}

void snova_fill_public_layout_mi(VECTOR (*layout)[vtl_len], const public_key_expand* pkx, int mi) {
	uint8_t* pj = (uint8_t*)layout;

	clear_public_layout_padding(pj);

	for (int ni = 0; ni < v_SNOVA; ++ni) {
		uint8_t* row0 = pj + (ni * rank) * PUBLIC_ROW_STRIDE;

		snova_pack_block_run(row0, PUBLIC_ROW_STRIDE, (const uint8_t*)pkx->map1.P11[mi][ni], v_SNOVA);
		snova_pack_block_run(row0 + v_SNOVA * rank, PUBLIC_ROW_STRIDE, (const uint8_t*)pkx->map1.P12[mi][ni], o_SNOVA);
	}

	for (int ni = 0; ni < o_SNOVA; ++ni) {
		uint8_t* row0 = pj + ((v_SNOVA + ni) * rank) * PUBLIC_ROW_STRIDE;

		snova_pack_block_run(row0, PUBLIC_ROW_STRIDE, (const uint8_t*)pkx->map1.P21[mi][ni], v_SNOVA);
		snova_pack_block_run(row0 + v_SNOVA * rank, PUBLIC_ROW_STRIDE, (const uint8_t*)pkx->P22[mi][ni], o_SNOVA);
	}
}
