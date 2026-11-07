#include "piop.h"
#if MQOM2_PARAM_WITH_STATISTICAL_BATCHING == 1
#include "xof.h"
#endif
#include "piop_cache.h"
#include "benchmark.h"
#include "expand_mq.h"

#include "fields_bitsliced.h"

#if MQOM2_PARAM_WITH_STATISTICAL_BATCHING == 1
static int ExpandBatchingChallenge(const uint8_t com[MQOM2_PARAM_DIGEST_SIZE], field_ext_elt Gamma[MQOM2_PARAM_ETA][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)]) {
	int ret = -1;
	uint32_t i;
	uint8_t stream[MQOM2_PARAM_ETA * BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)];
	xof_context DECL_VAR(xof_ctx);

	ret = xof_init(&xof_ctx);
	ERR(ret, err);
	ret = xof_update(&xof_ctx, (const uint8_t*) "\x08", 1);
	ERR(ret, err);
	ret = xof_update(&xof_ctx, com, MQOM2_PARAM_DIGEST_SIZE);
	ERR(ret, err);
	ret = xof_squeeze(&xof_ctx, stream, MQOM2_PARAM_ETA * BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU));
	ERR(ret, err);
	for (i = 0; i < MQOM2_PARAM_ETA; i++) {
		field_ext_parse(&stream[i * BYTE_SIZE_FIELD_EXT(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)], MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU, Gamma[i]);
	}

	ret = 0;
err:
	xof_clean_ctx(&xof_ctx);
	return ret;
}
#endif

/* Here, we perform as many 32 bitslice packets as we can, and perform the leftover in a non-bitsliced manner */
#ifndef FLOOR
#define FLOOR(x,y) ((y) * ((x) / (y)))
#endif
/* Compute what we have to compute in a non-bitslice manner */
#ifndef BITSLICE_HYBDID_LEFTOVER_LIMIT
#define BITSLICE_HYBDID_LEFTOVER_LIMIT 10 /* If we have less than BITSLICE_HYBDID_LEFTOVER_LIMIT elements, we perform regular non-bitslice operations */
#endif
#if BITSLICE_HYBDID_LEFTOVER_LIMIT > 32
#error "Error: BITSLICE_HYBDID_LEFTOVER_LIMIT > 32"
#endif

/* Some macros for ComputePz_xTau_bitslice */
#define BS_LEFTOVER ((MQOM2_PARAM_TAU + 1) - FLOOR((MQOM2_PARAM_TAU + 1), 32))
#if BS_LEFTOVER > BITSLICE_HYBDID_LEFTOVER_LIMIT
#define BITSLICE_LANES (MQOM2_PARAM_TAU + 1)
#define NON_BITSLICE_LEFTOVER 0
#define NON_BITSLICE_OFFSET 0
#else
#define BITSLICE_LANES FLOOR((MQOM2_PARAM_TAU + 1), 32)
#if BITSLICE_LANES > 0
#define NON_BITSLICE_LEFTOVER BS_LEFTOVER
#define NON_BITSLICE_OFFSET (BITSLICE_LANES - 1)
#else
/* In this specific case, no bitslice at all so t1 is computed in the regular way */
#define NON_BITSLICE_LEFTOVER MQOM2_PARAM_TAU
#define NON_BITSLICE_OFFSET 0
#endif
#endif

#if (BITSLICE_LANES == 0) && (NON_BITSLICE_LEFTOVER == 0)
#error "Error: (BITSLICE_LANES == 0) && (NON_BITSLICE_LEFTOVER == 0)"
#endif
static int ComputePz_xTau_bitslice(const field_ext_elt x0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], const field_base_elt x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)], const uint8_t mseed_eq[2 * MQOM2_PARAM_SEED_SIZE], field_ext_elt z0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)], field_ext_elt z1[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)]) {
	int ret = -1;
	uint32_t i, j, e;
	ExpandEquations_ctx EEctx = { 0 };

	/* Only use rows for A_hat and b_hat to save memory */
	field_ext_elt A_hat_row[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
	/* NOTE: we reuse the A_hat_row memory slot to save memory */
	field_ext_elt *b_hat_row = A_hat_row;

	field_ext_elt t1[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
#if (BITSLICE_LANES > 0)
	felt_ext_elt_bitsliced_t x0_bitsliced[BITSLICED_PACKING(MQOM2_PARAM_MQ_N, BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t t0_bitsliced[BITSLICED_PACKING(MQOM2_PARAM_MQ_N, BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t t1_x0_bitslice[BITSLICED_PACKING(1, BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t t0_x_bitslice[BITSLICED_PACKING(1, BITSLICE_LANES)];

	/* Aliasing to save stack space */
	felt_ext_elt_bitsliced_t *z0_bitsliced_i = t1_x0_bitslice;
	felt_ext_elt_bitsliced_t *z1_bitsliced_i = t0_x_bitslice;
	/* Transformation to bitslice */
	field_ext_bitslice_vect_pack_pre(x0_bitsliced, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
	for (e = 0; e < (uint32_t)(BITSLICE_LANES - 1); e++) {
		field_ext_bitslice_vect_pack(x0[e], x0_bitsliced, e, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
	}
	field_ext_bitslice_vect_pack_base(x, x0_bitsliced, BITSLICE_LANES-1, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
	field_ext_bitslice_vect_pack_post(x0_bitsliced, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
#endif

#if (NON_BITSLICE_LEFTOVER > 0)
	/* Leftover allocations */
        field_ext_elt t0[NON_BITSLICE_LEFTOVER][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
        field_ext_elt t1_x0[NON_BITSLICE_LEFTOVER];
        const field_ext_elt *x0_ptr[NON_BITSLICE_LEFTOVER];
        field_ext_elt *t0_ptr[NON_BITSLICE_LEFTOVER];
	field_ext_elt z_0i, z_1i;
        for (e = 0; e < NON_BITSLICE_LEFTOVER; e++) {
                x0_ptr[e] = x0[NON_BITSLICE_OFFSET + e];
        }
#endif
	/* Compute the equations expansion in a streaming way to save memory */
	ret = ExpandEquations_memopt_init(mseed_eq, &EEctx);
	ERR(ret, err);

	for (i = 0; i < MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU; i++) {
		for (j = 0; j < MQOM2_PARAM_MQ_N; j++) {
			__BENCHMARK_START__(BS_PIOP_EXPAND_MQ);
			/* Extract row from A_hat */
			ret = ExpandEquations_memopt_update(&EEctx, A_hat_row);
			ERR(ret, err);
			__BENCHMARK_STOP__(BS_PIOP_EXPAND_MQ);
#if (BITSLICE_LANES > 0)
			/* Simultaneously compute t0 for each tau repetition, and t1 */
			__BENCHMARK_START__(BS_PIOP_MAT_MUL_EXT);
			field_ext_bitslice_vect_mult_hybrid_public(
			    field_ext_bitslice_vect_get(t0_bitsliced, j, BITSLICE_LANES),
			    A_hat_row, x0_bitsliced, j + 1, BITSLICE_LANES 
			);
			__BENCHMARK_STOP__(BS_PIOP_MAT_MUL_EXT);
#else
                        /* Compute t1, common to all tau repetitions */
                        __BENCHMARK_START__(BS_PIOP_COMPUTE_T1);
                        t1[j] = field_ext_base_vect_mult(A_hat_row, x, j + 1);
                        __BENCHMARK_STOP__(BS_PIOP_COMPUTE_T1);
#endif
#if (NON_BITSLICE_LEFTOVER > 0)
			__BENCHMARK_START__(BS_PIOP_MAT_MUL_EXT);
			/* Compute the leftover in a non-bitsliced fashion */
			for(e = 0; e < NON_BITSLICE_LEFTOVER; e++){
				t0_ptr[e] = &t0[e][j];
			}
			field_ext_vect_mult_multiple_public(t0_ptr, A_hat_row, x0_ptr, j + 1, NON_BITSLICE_LEFTOVER);
			__BENCHMARK_STOP__(BS_PIOP_MAT_MUL_EXT);
#endif
		}
		/* Finish t1 computation with b_hat_row */
		__BENCHMARK_START__(BS_PIOP_EXPAND_MQ);
		ret = ExpandEquations_memopt_update(&EEctx, b_hat_row);
		ERR(ret, err);
		__BENCHMARK_STOP__(BS_PIOP_EXPAND_MQ);
		__BENCHMARK_START__(BS_PIOP_COMPUTE_T1);
#if (BITSLICE_LANES > 0)
		field_ext_bitslice_vect_unpack_pre(t0_bitsliced, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
		field_ext_bitslice_vect_unpack(t0_bitsliced, t1, BITSLICE_LANES - 1, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
		field_ext_bitslice_vect_pack_post(t0_bitsliced, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
		field_ext_vect_add(t1, b_hat_row, t1, MQOM2_PARAM_MQ_N);
#else
                field_ext_vect_add(t1, b_hat_row, t1, MQOM2_PARAM_MQ_N);
#endif
		__BENCHMARK_STOP__(BS_PIOP_COMPUTE_T1);
		__BENCHMARK_START__(BS_PIOP_COMPUTE_PZI);
#if (BITSLICE_LANES > 0)
		/* Compute the rest with bitslice */
		field_ext_bitslice_vect_mult_hybrid(t1_x0_bitslice, t1, x0_bitsliced, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
		field_ext_bitslice_vect_mult_hybrid_base(t0_x_bitslice, x, t0_bitsliced, MQOM2_PARAM_MQ_N, BITSLICE_LANES);
		field_ext_bitslice_add(
		    field_ext_bitslice_vect_get(z1_bitsliced_i, 0, BITSLICE_LANES),
		    t0_x_bitslice, t1_x0_bitslice,
		    BITSLICE_LANES
		);
		field_ext_bitslice_vect_mult(
		    field_ext_bitslice_vect_get(z0_bitsliced_i, 0, BITSLICE_LANES),
		    x0_bitsliced, t0_bitsliced, MQOM2_PARAM_MQ_N,
		    BITSLICE_LANES
		);
		field_ext_bitslice_vect_unpack_pre(z0_bitsliced_i, 1, BITSLICE_LANES);
		field_ext_bitslice_vect_unpack_pre(z1_bitsliced_i, 1, BITSLICE_LANES);
		for (e = 0; e < (uint32_t)(BITSLICE_LANES - 1); e++) {
			field_ext_bitslice_vect_unpack(z0_bitsliced_i, &z0[e][i], e, 1, BITSLICE_LANES);
			field_ext_bitslice_vect_unpack(z1_bitsliced_i, &z1[e][i], e, 1, BITSLICE_LANES);
		}
#endif
#if (NON_BITSLICE_LEFTOVER > 0)
		/* Compute the rest with non-bitslice (leftover) */
                for (e = 0; e < NON_BITSLICE_LEFTOVER; e++) {
			t1_x0[e] = field_ext_vect_mult(t1, x0[NON_BITSLICE_OFFSET + e], MQOM2_PARAM_MQ_N); /* t1^T x0[e] */
                }
                for (e = 0; e < NON_BITSLICE_LEFTOVER; e++) {
                        field_ext_elt t0_x = field_ext_base_vect_mult(t0[e], x, MQOM2_PARAM_MQ_N);   /* t0^T x[e] */
                        field_ext_vect_add(&t0_x, &t1_x0[e], &z_1i, 1);
                        field_ext_vect_pack(z_1i, z1[NON_BITSLICE_OFFSET + e], i);
                }
                for (e = 0; e < NON_BITSLICE_LEFTOVER; e++) {
                        z_0i = field_ext_vect_mult(t0[e], x0[NON_BITSLICE_OFFSET + e], MQOM2_PARAM_MQ_N);
                        field_ext_vect_pack(z_0i, z0[NON_BITSLICE_OFFSET + e], i);
                }
#endif
		__BENCHMARK_STOP__(BS_PIOP_COMPUTE_PZI);
	}

	ret = 0;
err:
	ExpandEquations_memopt_final(&EEctx);
	return ret;
}

int ComputePAlpha_bitslice(const uint8_t com[MQOM2_PARAM_DIGEST_SIZE], const field_ext_elt x0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], const field_ext_elt u0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], const field_ext_elt u1[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], const field_base_elt x[FIELD_BASE_PACKING(MQOM2_PARAM_MQ_N)], const uint8_t mseed_eq[2 * MQOM2_PARAM_SEED_SIZE], field_ext_elt alpha0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], field_ext_elt alpha1[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)]) {
	int ret = -1;
	uint32_t e;

	/* Initialize the PIOP cache for t1 */
	__BENCHMARK_START__(BS_PIOP_EXPAND_BATCHING_MAT);
#if MQOM2_PARAM_WITH_STATISTICAL_BATCHING == 1
	uint32_t i;
	field_ext_elt Gamma[MQOM2_PARAM_ETA][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)];
	ret = ExpandBatchingChallenge(com, Gamma);
	ERR(ret, err);
#else
	(void) com;
#endif
	__BENCHMARK_STOP__(BS_PIOP_EXPAND_BATCHING_MAT);

	field_ext_elt z0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)], z1[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)];
	ret = ComputePz_xTau_bitslice(x0, x, mseed_eq, z0, z1);
	ERR(ret, err);
	for (e = 0; e < MQOM2_PARAM_TAU; e++) {
		__BENCHMARK_START__(BS_PIOP_BATCH_AND_MASK);
#if MQOM2_PARAM_WITH_STATISTICAL_BATCHING == 1
		field_ext_elt tmp[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)];
		for (i = 0; i < MQOM2_PARAM_ETA; i++) {
			field_ext_vect_pack(
			    field_ext_vect_mult(Gamma[i], z0[e], MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU),
			    tmp, i
			);
		}
		field_ext_vect_add(tmp, u0[e], alpha0[e], MQOM2_PARAM_ETA);
		for (i = 0; i < MQOM2_PARAM_ETA; i++) {
			field_ext_vect_pack(
			    field_ext_vect_mult(Gamma[i], z1[e], MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU),
			    tmp, i
			);
		}
		field_ext_vect_add(tmp, u1[e], alpha1[e], MQOM2_PARAM_ETA);
#else
		field_ext_vect_add(z0[e], u0[e], alpha0[e], MQOM2_PARAM_ETA);
		field_ext_vect_add(z1[e], u1[e], alpha1[e], MQOM2_PARAM_ETA);
#endif
		__BENCHMARK_STOP__(BS_PIOP_BATCH_AND_MASK);
	}

	ret = 0;
err:
	return ret;
}

/***************************************************************/
/***************************************************************/
/* Some macros for ComputePzEval_xTau_bitslice */
#define EVAL_BS_LEFTOVER (MQOM2_PARAM_TAU - FLOOR(MQOM2_PARAM_TAU, 32))
#if EVAL_BS_LEFTOVER > BITSLICE_HYBDID_LEFTOVER_LIMIT
#define EVAL_BITSLICE_LANES MQOM2_PARAM_TAU
#define EVAL_NON_BITSLICE_LEFTOVER 0
#define EVAL_NON_BITSLICE_OFFSET 0
#else
#define EVAL_BITSLICE_LANES FLOOR(MQOM2_PARAM_TAU, 32)
#define EVAL_NON_BITSLICE_LEFTOVER EVAL_BS_LEFTOVER
#define EVAL_NON_BITSLICE_OFFSET EVAL_BITSLICE_LANES
#endif

static int ComputePzEval_xTau_bitslice(const field_ext_elt r[MQOM2_PARAM_TAU], const field_ext_elt v_x[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], const uint8_t mseed_eq[2 * MQOM2_PARAM_SEED_SIZE], const field_ext_elt y[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)], field_ext_elt v_z[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)]) {
	int ret = -1;
	uint32_t i, j, e;
	ExpandEquations_ctx EEctx = { 0 };

	/* Only use rows for A_hat and b_hat to save memory */
	field_ext_elt A_hat_row[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
	/* NOTE: we reuse the A_hat_row memory slot to save memory */
	field_ext_elt *b_hat_row = A_hat_row;

#if (EVAL_BITSLICE_LANES > 0)
	felt_ext_elt_bitsliced_t r_bitsliced[BITSLICED_PACKING(1, EVAL_BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t r2_bitsliced[BITSLICED_PACKING(1, EVAL_BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t v_x_bitsliced[BITSLICED_PACKING(MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t v_t_bitsliced[BITSLICED_PACKING(MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES)];
	felt_ext_elt_bitsliced_t tmp_bitsliced[BITSLICED_PACKING(MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES)];
	/* Aliasing to save memory slot */
	felt_ext_elt_bitsliced_t *y_r2_bitsliced_i = tmp_bitsliced; // [BITSLICED_PACKING(1, EVAL_BITSLICE_LANES)]
	felt_ext_elt_bitsliced_t v_z_bitsliced_i[BITSLICED_PACKING(1, EVAL_BITSLICE_LANES)];

	field_ext_bitslice_vect_pack_pre(r_bitsliced, 1, EVAL_BITSLICE_LANES);
	field_ext_bitslice_vect_pack_pre(v_x_bitsliced, MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES);
	for (e = 0; e < EVAL_BITSLICE_LANES; e++) {
		field_ext_bitslice_vect_pack(&r[e], r_bitsliced, e, 1, EVAL_BITSLICE_LANES);
		field_ext_bitslice_vect_pack(v_x[e], v_x_bitsliced, e, MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES);
	}
	field_ext_bitslice_vect_pack_post(r_bitsliced, 1, EVAL_BITSLICE_LANES);
	field_ext_bitslice_vect_pack_post(v_x_bitsliced, MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES);
	field_ext_bitslice_mult(r2_bitsliced, r_bitsliced, r_bitsliced, EVAL_BITSLICE_LANES);
#endif
#if (EVAL_NON_BITSLICE_LEFTOVER > 0)
        field_ext_elt v_t[EVAL_BS_LEFTOVER][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
        field_ext_elt tmp[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)];
        field_ext_elt y_r2[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)];
        field_ext_elt v_zi;

        const field_ext_elt* vx_ptr[EVAL_BS_LEFTOVER];
        field_ext_elt* vt_ptr[EVAL_BS_LEFTOVER];

        for (e = 0; e < EVAL_BS_LEFTOVER; e++) {
                vx_ptr[e] = v_x[EVAL_NON_BITSLICE_OFFSET + e];
        }
#endif
	/* Compute the equations expansion in a streaming way to save memory */
	ret = ExpandEquations_memopt_init(mseed_eq, &EEctx);
	ERR(ret, err);

	for (i = 0; i < MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU; i++) {
		/* Perform operations row by row for A_hat */
		for (j = 0; j < MQOM2_PARAM_MQ_N; j++) {
			ret = ExpandEquations_memopt_update(&EEctx, A_hat_row);
			ERR(ret, err);
#if (EVAL_BITSLICE_LANES > 0)
			field_ext_bitslice_vect_mult_hybrid_public(
			    field_ext_bitslice_vect_get(v_t_bitsliced, j, EVAL_BITSLICE_LANES),
			    A_hat_row, v_x_bitsliced, j + 1, EVAL_BITSLICE_LANES
			);
#endif
#if (EVAL_NON_BITSLICE_LEFTOVER > 0)
                        for (e = 0; e < EVAL_BS_LEFTOVER; e++) {
                                vt_ptr[e] = &v_t[e][j];
                        }
                        field_ext_vect_mult_multiple_public(vt_ptr, A_hat_row, vx_ptr, j + 1, EVAL_BS_LEFTOVER);
#endif
		}
		/* Generate and add b_hat row */
		ret = ExpandEquations_memopt_update(&EEctx, b_hat_row);
		ERR(ret, err);
#if (EVAL_BITSLICE_LANES > 0)
		field_ext_bitslice_const_vect_mult_hybrid(tmp_bitsliced, r_bitsliced, b_hat_row, MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES);
		field_ext_bitslice_vect_add(v_t_bitsliced, v_t_bitsliced, tmp_bitsliced, MQOM2_PARAM_MQ_N, EVAL_BITSLICE_LANES);
		/* Compute v_{z,i} = P_{z,i}(r) = v_t^T v_x - y_i r^2 */
		field_ext_bitslice_vect_mult(
		    field_ext_bitslice_vect_get(v_z_bitsliced_i, 0, EVAL_BITSLICE_LANES),
		    v_t_bitsliced, v_x_bitsliced, MQOM2_PARAM_MQ_N,
		    EVAL_BITSLICE_LANES
		);
		/**/
		field_ext_bitslice_const_vect_mult_hybrid(y_r2_bitsliced_i, r2_bitsliced, &y[i], 1, EVAL_BITSLICE_LANES);
		field_ext_bitslice_vect_add(v_z_bitsliced_i, v_z_bitsliced_i, y_r2_bitsliced_i, 1, EVAL_BITSLICE_LANES);
		field_ext_bitslice_vect_unpack_pre(v_z_bitsliced_i, 1, EVAL_BITSLICE_LANES);
		for (e = 0; e < EVAL_BITSLICE_LANES; e++) {
			field_ext_bitslice_vect_unpack(v_z_bitsliced_i, &v_z[e][i], e, 1, EVAL_BITSLICE_LANES);
		}
#endif
#if (EVAL_NON_BITSLICE_LEFTOVER > 0)
                for (e = 0; e < EVAL_BS_LEFTOVER; e++) {
                        field_ext_constant_vect_mult(r[EVAL_NON_BITSLICE_OFFSET + e], b_hat_row, tmp, MQOM2_PARAM_MQ_N);
                        field_ext_vect_add(v_t[e], tmp, v_t[e], MQOM2_PARAM_MQ_N);
                        /* Compute v_{z,i} = P_{z,i}(r) = v_t^T v_x - y_i r^2 */
                        v_zi = field_ext_vect_mult(v_t[e], v_x[EVAL_NON_BITSLICE_OFFSET + e], MQOM2_PARAM_MQ_N);
                        field_ext_vect_pack(v_zi, v_z[EVAL_NON_BITSLICE_OFFSET + e], i);
                }
#endif
	}

#if (EVAL_NON_BITSLICE_LEFTOVER > 0)
        for (e = 0; e < EVAL_BS_LEFTOVER; e++) {
                field_ext_elt r2 = field_ext_mult(r[EVAL_NON_BITSLICE_OFFSET + e], r[EVAL_NON_BITSLICE_OFFSET + e]);
                field_ext_constant_vect_mult(r2, y, y_r2, MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU);
                field_ext_vect_add(v_z[EVAL_NON_BITSLICE_OFFSET + e], y_r2, v_z[EVAL_NON_BITSLICE_OFFSET + e], MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU);
        }
#endif
	ret = 0;
err:
	ExpandEquations_memopt_final(&EEctx);
	return ret;
}

int RecomputePAlpha_bitslice(const uint8_t com[MQOM2_PARAM_DIGEST_SIZE], const field_ext_elt alpha1[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], const uint16_t i_star[MQOM2_PARAM_TAU], const field_ext_elt x_eval[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_N)], const field_ext_elt u_eval[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)], const uint8_t mseed_eq[2 * MQOM2_PARAM_SEED_SIZE], const field_ext_elt y[FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)], field_ext_elt alpha0[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_ETA)]) {
	int ret = -1;
	uint32_t e;

#if MQOM2_PARAM_WITH_STATISTICAL_BATCHING == 1
	uint32_t i;
	field_ext_elt Gamma[MQOM2_PARAM_ETA][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)];
	ret = ExpandBatchingChallenge(com, Gamma);
	ERR(ret, err);
#else
	(void) com;
#endif

	field_ext_elt v_z[MQOM2_PARAM_TAU][FIELD_EXT_PACKING(MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU)];
	field_ext_elt v_alpha[FIELD_EXT_PACKING(MQOM2_PARAM_ETA)];
	field_ext_elt r[MQOM2_PARAM_TAU];
	for (e = 0; e < MQOM2_PARAM_TAU; e++) {
		r[e] = get_evaluation_point(i_star[e]);
	}
	ret = ComputePzEval_xTau_bitslice(r, x_eval, mseed_eq, y, v_z);
	ERR(ret, err);
	for (e = 0; e < MQOM2_PARAM_TAU; e++) {
#if MQOM2_PARAM_WITH_STATISTICAL_BATCHING == 1
		for (i = 0; i < MQOM2_PARAM_ETA; i++) {
			field_ext_vect_pack(
			    field_ext_vect_mult(Gamma[i], v_z[e], MQOM2_PARAM_MQ_M / MQOM2_PARAM_MU),
			    v_alpha, i
			);
		}
		field_ext_vect_add(v_alpha, u_eval[e], v_alpha, MQOM2_PARAM_ETA);
		field_ext_constant_vect_mult(r[e], alpha1[e], alpha0[e], MQOM2_PARAM_ETA);
		field_ext_vect_add(v_alpha, alpha0[e], alpha0[e], MQOM2_PARAM_ETA);
#else
		field_ext_vect_add(v_z[e], u_eval[e], v_alpha, MQOM2_PARAM_ETA);
		field_ext_constant_vect_mult(r[e], alpha1[e], alpha0[e], MQOM2_PARAM_ETA);
		field_ext_vect_add(v_alpha, alpha0[e], alpha0[e], MQOM2_PARAM_ETA);
#endif
	}

	ret = 0;
err:
	return ret;
}
