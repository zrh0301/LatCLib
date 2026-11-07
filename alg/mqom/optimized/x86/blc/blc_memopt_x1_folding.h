#ifndef __BLC_MEMOPT_X1_FOLDING_H__
#define __BLC_MEMOPT_X1_FOLDING_H__

#include "fields.h"
#include "enc.h"
#include "blc_memopt_x1.h"

#ifdef BLC_SEEDEXPAND_CACHE
#define BLC_NB_SEEDEXPAND_ENC_CTX ((PRG_BLC_SIZE+MQOM2_PARAM_SEED_SIZE-1) / MQOM2_PARAM_SEED_SIZE)
#endif

typedef struct {
	uint8_t acc[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N) + BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_ETA)];
	uint8_t data[MQOM2_PARAM_NB_EVALS_LOG][BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N) + BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_ETA)];
	uint32_t e;
	const uint8_t* salt;
#ifdef BLC_SEEDEXPAND_CACHE
	enc_ctx_ecb DECL_VAR(enc_ctx[BLC_NB_SEEDEXPAND_ENC_CTX]);
#endif
} folding_sign_t;

typedef struct {
	uint8_t acc[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N) + BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_ETA)];
	uint8_t data[MQOM2_PARAM_NB_EVALS_LOG][BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N) + BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_ETA)];
	uint32_t e;
	const uint8_t* salt;
#ifdef BLC_SEEDEXPAND_CACHE
	enc_ctx_pub_ecb DECL_VAR(enc_ctx[BLC_NB_SEEDEXPAND_ENC_CTX]);
#endif
} folding_verify_t;

#if PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE != 0
static inline void xor_blocks_partial(const uint8_t seed_in[PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE], const uint8_t delta[PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE], uint8_t seed_out[PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE]) {
	unsigned int i;

	for (i = 0; i < PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE; i++) {
		seed_out[i] = seed_in[i] ^ delta[i];
	}

	return;
}
#endif

static inline int InitializeFolding_sign(folding_sign_t* folding, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], uint32_t e) {
	int ret = -1;

	memset(folding->acc, 0, sizeof(folding->acc));
	memset(folding->data, 0, sizeof(folding->data));
	folding->e = e;
	folding->salt = salt;

#ifdef BLC_SEEDEXPAND_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	for (uint32_t j = 0; j < BLC_NB_SEEDEXPAND_ENC_CTX; j++) {
		TweakSalt(salt, tweaked_salt, 3, e, j);
		ret = enc_key_sched_ecb(&folding->enc_ctx[j], tweaked_salt);
		ERR(ret, err);
	}
#endif

	ret = 0;
#ifdef BLC_SEEDEXPAND_CACHE
err:
#endif
	return ret;
}

static inline int InitializeFolding_verify(folding_verify_t* folding, const uint8_t salt[MQOM2_PARAM_SALT_SIZE], uint32_t e) {
	int ret = -1;

	memset(folding->acc, 0, sizeof(folding->acc));
	memset(folding->data, 0, sizeof(folding->data));
	folding->e = e;
	folding->salt = salt;

#ifdef BLC_SEEDEXPAND_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	for (uint32_t j = 0; j < BLC_NB_SEEDEXPAND_ENC_CTX; j++) {
		TweakSalt(salt, tweaked_salt, 3, e, j);
		ret = enc_key_sched_pub_ecb(&folding->enc_ctx[j], tweaked_salt);
		ERR(ret, err);
	}
#endif

	ret = 0;
#ifdef BLC_SEEDEXPAND_CACHE
err:
#endif
	return ret;
}

static inline void folding_sign_clean_ctx(folding_sign_t* ctx) {
	(void) ctx;
#ifdef BLC_SEEDEXPAND_CACHE
	for (uint32_t j = 0; j < BLC_NB_SEEDEXPAND_ENC_CTX; j++) {
		enc_clean_ctx_ecb(&ctx->enc_ctx[j]);
	}
#endif
}

static inline void folding_verify_clean_ctx(folding_verify_t* ctx) {
	(void) ctx;
#ifdef BLC_SEEDEXPAND_CACHE
	for (uint32_t j = 0; j < BLC_NB_SEEDEXPAND_ENC_CTX; j++) {
		enc_clean_ctx_pub_ecb(&ctx->enc_ctx[j]);
	}
#endif
}

static inline int SeedExpandThenAccumulate_sign(folding_sign_t* folding, uint32_t i, const uint8_t lseed[][MQOM2_PARAM_SEED_SIZE]) {
	int ret = -1;
	uint32_t j, i_;

#ifndef BLC_SEEDEXPAND_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	enc_ctx_ecb DECL_VAR(enc_ctx_);
	enc_ctx_ecb* enc_ctx_ptr = &enc_ctx_;
#else
	enc_ctx_ecb* enc_ctx_ptr;
#endif

	uint8_t linortho_seed[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	uint8_t out_data[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	uint8_t pos[BLC_NB_LEAF_SEEDS_IN_PARALLEL];

	/* Compute Psi(seed) once and for all */
	for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		LinOrtho(lseed[i_], linortho_seed[i_]);
		pos[i_] = get_gray_code_bit_position(i+i_);
	}

	__BENCHMARK_START__(BS_BLC_ARITH);
	for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		xor_blocks(&folding->acc[0], lseed[i_], &folding->acc[0]);
		xor_blocks(&folding->data[pos[i_]][0], &folding->acc[0], &folding->data[pos[i_]][0]);
	}
	__BENCHMARK_STOP__(BS_BLC_ARITH);

	for (j = 0; j < (PRG_BLC_SIZE / MQOM2_PARAM_SEED_SIZE); j++) {
		/* Key schedule */
		__BENCHMARK_START__(BS_BLC_PRG);
	#ifndef BLC_SEEDEXPAND_CACHE
		TweakSalt(folding->salt, tweaked_salt, 3, folding->e, j);
		ret = enc_key_sched_ecb(enc_ctx_ptr, tweaked_salt);
		ERR(ret, err);
	#else
		enc_ctx_ptr = &folding->enc_ctx[j];
	#endif
		
		/* Encryption */
		ret = enc_encrypt_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) out_data);
		ERR(ret, err);
		__BENCHMARK_STOP__(BS_BLC_PRG);

		for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
			/* Xor with LinOrtho seed */
			xor_blocks(out_data[i_], linortho_seed[i_], out_data[i_]);

			__BENCHMARK_START__(BS_BLC_ARITH);
			xor_blocks(&folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], out_data[i_], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE]);
			xor_blocks(&folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE]);
			__BENCHMARK_STOP__(BS_BLC_ARITH);
		}
	}
	/* Deal with the possible leftover incomplete block */
#if PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE != 0
	/* Key schedule */
	__BENCHMARK_START__(BS_BLC_PRG);
#ifndef BLC_SEEDEXPAND_CACHE
	TweakSalt(folding->salt, tweaked_salt, 3, folding->e, j);
	ret = enc_key_sched_ecb(enc_ctx_ptr, tweaked_salt);
	ERR(ret, err);
#else
	enc_ctx_ptr = &folding->enc_ctx[BLC_NB_SEEDEXPAND_ENC_CTX-1];
#endif

	/* Encryption */
	ret = enc_encrypt_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) out_data);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_PRG);

	for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		/* Xor with LinOrtho seed */
		xor_blocks(out_data[i_], linortho_seed[i_], out_data[i_]);

		__BENCHMARK_START__(BS_BLC_ARITH);
		xor_blocks_partial(&folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], out_data[i_], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE]);
		xor_blocks_partial(&folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE]);
		__BENCHMARK_STOP__(BS_BLC_ARITH);
	}
#endif

	ret = 0;
err:
#ifndef BLC_SEEDEXPAND_CACHE
	enc_clean_ctx_ecb(&enc_ctx_);
#endif
	return ret;
}

static inline int SeedExpandThenAccumulate_verify(folding_verify_t* folding, uint32_t i, const uint8_t lseed[][MQOM2_PARAM_SEED_SIZE], uint32_t hidden_index) {
	int ret = -1;
	uint32_t j, i_;

#ifndef BLC_SEEDEXPAND_CACHE
	uint8_t tweaked_salt[MQOM2_PARAM_SALT_SIZE];
	enc_ctx_pub_ecb DECL_VAR(enc_ctx_);
	enc_ctx_pub_ecb* enc_ctx_ptr = &enc_ctx_;
#else
	enc_ctx_pub_ecb* enc_ctx_ptr;
#endif

	uint8_t linortho_seed[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	uint8_t out_data[BLC_NB_LEAF_SEEDS_IN_PARALLEL][MQOM2_PARAM_SEED_SIZE];
	uint8_t pos[BLC_NB_LEAF_SEEDS_IN_PARALLEL];

	/* Compute Psi(seed) once and for all */
	for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		LinOrtho(lseed[i_], linortho_seed[i_]);
		pos[i_] = get_gray_code_bit_position(i+i_);
	}

	__BENCHMARK_START__(BS_BLC_ARITH);
	for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		if (hidden_index != i + i_) {
			xor_blocks(&folding->acc[0], lseed[i_], &folding->acc[0]);
		}
		xor_blocks(&folding->data[pos[i_]][0], &folding->acc[0], &folding->data[pos[i_]][0]);
	}
	__BENCHMARK_STOP__(BS_BLC_ARITH);

	for (j = 0; j < (PRG_BLC_SIZE / MQOM2_PARAM_SEED_SIZE); j++) {
		/* Key schedule */
		__BENCHMARK_START__(BS_BLC_PRG);
	#ifndef BLC_SEEDEXPAND_CACHE
		TweakSalt(folding->salt, tweaked_salt, 3, folding->e, j);
		ret = enc_key_sched_pub_ecb(enc_ctx_ptr, tweaked_salt);
		ERR(ret, err);
	#else
		enc_ctx_ptr = &folding->enc_ctx[j];
	#endif
	
		/* Encryption */
		ret = enc_encrypt_pub_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) out_data);
		ERR(ret, err);
		__BENCHMARK_STOP__(BS_BLC_PRG);

		for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
			/* Xor with LinOrtho seed */
			xor_blocks(out_data[i_], linortho_seed[i_], out_data[i_]);

			if (hidden_index == i + i_) {
				memset(out_data[i_], 0, MQOM2_PARAM_SEED_SIZE);
			}

			__BENCHMARK_START__(BS_BLC_ARITH);
			xor_blocks(&folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], out_data[i_], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE]);
			xor_blocks(&folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE]);
			__BENCHMARK_STOP__(BS_BLC_ARITH);
		}
	}
	/* Deal with the possible leftover incomplete block */
#if PRG_BLC_SIZE % MQOM2_PARAM_SEED_SIZE != 0
	/* Key schedule */
	__BENCHMARK_START__(BS_BLC_PRG);
#ifndef BLC_SEEDEXPAND_CACHE
	TweakSalt(folding->salt, tweaked_salt, 3, folding->e, j);
	ret = enc_key_sched_pub_ecb(enc_ctx_ptr, tweaked_salt);
	ERR(ret, err);
#else
	enc_ctx_ptr = &folding->enc_ctx[BLC_NB_SEEDEXPAND_ENC_CTX-1];
#endif

	/* Encryption */
	ret = enc_encrypt_pub_ecb(enc_ctx_ptr, BLC_NB_LEAF_SEEDS_IN_PARALLEL, (uint8_t*) lseed, (uint8_t*) out_data);
	ERR(ret, err);
	__BENCHMARK_STOP__(BS_BLC_PRG);

	for (i_ = 0; i_ < BLC_NB_LEAF_SEEDS_IN_PARALLEL; i_++) {
		/* Xor with LinOrtho seed */
		xor_blocks(out_data[i_], linortho_seed[i_], out_data[i_]);

		if (hidden_index == i + i_) {
			memset(out_data[i_], 0, MQOM2_PARAM_SEED_SIZE);
		}

		__BENCHMARK_START__(BS_BLC_ARITH);
		xor_blocks_partial(&folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], out_data[i_], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE]);
		xor_blocks_partial(&folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->acc[(j+1)*MQOM2_PARAM_SEED_SIZE], &folding->data[pos[i_]][(j+1)*MQOM2_PARAM_SEED_SIZE]);
		__BENCHMARK_STOP__(BS_BLC_ARITH);
	}
#endif

	ret = 0;
err:
#ifndef BLC_SEEDEXPAND_CACHE
	enc_clean_ctx_pub_ecb(&enc_ctx_);
#endif
	return ret;
}

static inline void FinalizeFolding_sign(const folding_sign_t* folding, const field_base_elt x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)], uint8_t partial_delta_x[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N)-MQOM2_PARAM_SEED_SIZE], field_ext_elt x0[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], field_ext_elt u0[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], field_ext_elt u1[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)]) {
	field_base_elt bar_x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)];
	field_ext_elt bar_u[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)];
	field_ext_elt tmp_n[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
	/* Alias tmp_eta to save stack space */
	field_ext_elt *tmp_eta = bar_u;
	/* Alias acc_x to save stack space */
	field_base_elt *acc_x = bar_x;

	memset(x0, 0, BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_MQ_N));
	for (uint32_t j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		field_base_parse(folding->data[j], MQOM2_PARAM_MQ_N, bar_x);
		field_ext_base_constant_vect_mult((1 << j), bar_x, tmp_n, MQOM2_PARAM_MQ_N);
		field_ext_vect_add(x0, tmp_n, x0, MQOM2_PARAM_MQ_N);
	}

	memset(u0, 0, BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_ETA));
	for (uint32_t j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		field_ext_parse(folding->data[j] + BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N), MQOM2_PARAM_ETA, bar_u);
		field_ext_constant_vect_mult((1 << j), bar_u, tmp_eta, MQOM2_PARAM_ETA);
		field_ext_vect_add(u0, tmp_eta, u0, MQOM2_PARAM_ETA);
	}

	field_base_parse(folding->acc, MQOM2_PARAM_MQ_N, acc_x);
	field_ext_parse(folding->acc + BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N), MQOM2_PARAM_ETA, u1);

	field_base_elt delta_x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)];
	uint8_t serialized_delta_x[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N)];
	field_base_vect_add(x, acc_x, delta_x, MQOM2_PARAM_MQ_N);
	field_base_serialize(delta_x, MQOM2_PARAM_MQ_N, serialized_delta_x);
	memcpy(partial_delta_x, serialized_delta_x + MQOM2_PARAM_SEED_SIZE, BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N) - MQOM2_PARAM_SEED_SIZE);
}

static inline void FinalizeFolding_verify(const folding_verify_t* folding, uint16_t i_star, const uint8_t partial_delta_x[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N)-MQOM2_PARAM_SEED_SIZE], field_ext_elt x_eval[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], field_ext_elt u_eval[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)]) {
	field_base_elt bar_x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)];
	field_ext_elt bar_u[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)];
	field_ext_elt tmp_n[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
	/* Alias tmp_eta to save stack space */
	field_ext_elt *tmp_eta = bar_u;
	/* Alias acc_x to save stack space */
	field_base_elt *acc_x = bar_x;

	field_ext_elt r = get_evaluation_point(i_star);

	field_base_elt delta_x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)];
	uint8_t serialized_delta_x[BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N)];
	memset(serialized_delta_x, 0, MQOM2_PARAM_SEED_SIZE);
	memcpy(serialized_delta_x + MQOM2_PARAM_SEED_SIZE, partial_delta_x, BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N) - MQOM2_PARAM_SEED_SIZE);
	field_base_parse(serialized_delta_x, MQOM2_PARAM_MQ_N, delta_x);

	memset(x_eval, 0, BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_MQ_N));
	for (uint32_t j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		field_base_parse(folding->data[j], MQOM2_PARAM_MQ_N, bar_x);
		field_ext_base_constant_vect_mult((1 << j), bar_x, tmp_n, MQOM2_PARAM_MQ_N);
		field_ext_vect_add(x_eval, tmp_n, x_eval, MQOM2_PARAM_MQ_N);
	}
	field_base_parse(folding->acc, MQOM2_PARAM_MQ_N, acc_x);
	field_base_vect_add(acc_x, delta_x, acc_x, MQOM2_PARAM_MQ_N);
	field_ext_base_constant_vect_mult(r, acc_x, tmp_n, MQOM2_PARAM_MQ_N);
	field_ext_vect_add(x_eval, tmp_n, x_eval, MQOM2_PARAM_MQ_N);

	memset(u_eval, 0, BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_ETA));
	for (uint32_t j = 0; j < MQOM2_PARAM_NB_EVALS_LOG; j++) {
		field_ext_parse(folding->data[j] + BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N), MQOM2_PARAM_ETA, bar_u);
		field_ext_constant_vect_mult((1 << j), bar_u, tmp_eta, MQOM2_PARAM_ETA);
		field_ext_vect_add(u_eval, tmp_eta, u_eval, MQOM2_PARAM_ETA);
	}
	field_ext_parse(folding->acc + BYTE_SIZE_FIELD_BASE(MQOM2_PARAM_MQ_N), MQOM2_PARAM_ETA, tmp_eta);
	field_ext_constant_vect_mult(r, tmp_eta, tmp_eta, MQOM2_PARAM_ETA);
	field_ext_vect_add(u_eval, tmp_eta, u_eval, MQOM2_PARAM_ETA);
}


#endif /* __BLC_MEMOPT_X1_FOLDING_H__ */
