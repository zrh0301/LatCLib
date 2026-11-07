#ifndef SNOVA_PUBLIC_LAYOUT_H
#define SNOVA_PUBLIC_LAYOUT_H

#include "snova_transform.h"

#define GF16s_hash (o_SNOVA * lsq_SNOVA)
#define GF16s_signature (n_SNOVA * lsq_SNOVA)
#define bytes_hash ((GF16s_hash + 1) >> 1)
#define bytes_signature ((GF16s_signature + 1) >> 1)
#define bytes_sig_with_salt (bytes_signature + bytes_salt)

#define GF16s_prng_public                                                                          \
    (sq_rank * (2 * (m_SNOVA * alpha_SNOVA) + m_SNOVA * (n_SNOVA * n_SNOVA - m_SNOVA * m_SNOVA)) + \
     rank * 2 * m_SNOVA * alpha_SNOVA)
#define bytes_prng_public ((GF16s_prng_public + 1) >> 1)

#define GF16s_prng_private (v_SNOVA * o_SNOVA * rank)
#define bytes_prng_private ((GF16s_prng_private + 1) >> 1)

#define bytes_pk (seed_length_public + ((m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA + 1) >> 1))

#define bytes_sk                                                                                                      \
    (((sq_rank * (4 * m_SNOVA * alpha_SNOVA + m_SNOVA * (v_SNOVA * v_SNOVA + v_SNOVA * o_SNOVA + o_SNOVA * v_SNOVA) + \
                  v_SNOVA * o_SNOVA) +                                                                                \
       1) >>                                                                                                          \
      1) +                                                                                                            \
     seed_length_public + seed_length_private)

typedef gf16m_t P11_t[SNOVA_m1][SNOVA_v][SNOVA_v];
typedef gf16m_t P12_t[SNOVA_m1][SNOVA_v][SNOVA_o];
typedef gf16m_t P21_t[SNOVA_m1][SNOVA_o][SNOVA_v];
typedef gf16m_t Aalpha_t[SNOVA_m1][SNOVA_alpha];
typedef gf16m_t Balpha_t[SNOVA_m1][SNOVA_alpha];
typedef gf16m_t Qalpha1_t[SNOVA_m1][SNOVA_alpha];
typedef gf16m_t Qalpha2_t[SNOVA_m1][SNOVA_alpha];
typedef gf16_t Qalpha_coeff_t[SNOVA_m1][SNOVA_alpha][SNOVA_l];

typedef struct {
	P11_t P11;
	P12_t P12;
	P21_t P21;
	Aalpha_t Aalpha;
	Balpha_t Balpha;
	Qalpha1_t Qalpha1;
	Qalpha2_t Qalpha2;
} map_group1;

typedef gf16m_t T12_t[SNOVA_v][SNOVA_o];
typedef gf16m_t F11_t[SNOVA_m1][SNOVA_v][SNOVA_v];
typedef gf16m_t F12_t[SNOVA_m1][SNOVA_v][SNOVA_o];
typedef gf16m_t F21_t[SNOVA_m1][SNOVA_o][SNOVA_v];

typedef struct {
	F11_t F11;
	F12_t F12;
	F21_t F21;
} map_group2;

typedef struct {
	Aalpha_t Aalpha;
	Balpha_t Balpha;
	Qalpha1_t Qalpha1;
	Qalpha2_t Qalpha2;
	T12_t T12;
	F11_t F11;
	F12_t F12;
	F21_t F21;
	uint8_t pt_public_key_seed[SEED_LENGTH_PUBLIC];
	uint8_t pt_private_key_seed[SEED_LENGTH_PRIVATE];
} sk_gf16;

typedef gf16m_t P22_t[SNOVA_m1][SNOVA_o][SNOVA_o];
typedef uint8_t P22_byte_t[(SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2 + 1) >> 1];

typedef struct {
	uint8_t pt_public_key_seed[SEED_LENGTH_PUBLIC];
	P22_byte_t P22;
} public_key;

typedef struct {
	uint8_t pt_public_key_seed[SEED_LENGTH_PUBLIC];
	P22_t P22;
	map_group1 map1;
#if !FIXED_ABQ
	Qalpha_coeff_t q1coeff;
	Qalpha_coeff_t q2coeff;
#endif
} public_key_expand;

typedef struct {
	map_group1 map1;
	T12_t T12;
	map_group2 map2;
	public_key pk;
} snova_key_elems;

void snova_pk_pack(uint8_t* pk, const snova_key_elems* key_elems);

void snova_sk_pack(uint8_t* esk, const snova_key_elems* key_elems, const uint8_t* pt_private_key_seed);

void snova_sk_unpack(sk_gf16* skupk, const uint8_t* esk);

void snova_input_P22(uint8_t* P22_gf16s, const uint8_t* P22_bytes);

void snova_fill_public_layout_mi(VECTOR (*layout)[vtl_len], const public_key_expand* pkx, int mi);

#endif
