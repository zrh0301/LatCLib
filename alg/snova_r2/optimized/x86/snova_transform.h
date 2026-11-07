#ifndef SNOVA_TRANSFORM_H
#define SNOVA_TRANSFORM_H

#include "snova_vecops.h"

int snova_bytes_to_gf16(const uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s);

void snova_gf16_to_bytes(uint8_t* byte_array, const uint8_t* gf16_array, int num_of_GF16s);

void snova_bytes_to_gf16_split(const uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s);

void snova_gf16_to_bytes_merge(uint8_t* byte_array, const uint8_t* gf16_array, int num_of_GF16s);

#if SNOVA_RANK4_FAST

static const uint8_t snova_run_idx[16] = {
	0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51
};

static const uint8_t snova_run_tr_idx[16] = {
	0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60
};

static const uint8_t snova_tr_block_idx[16] = {
	0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15
};

static const uint8_t snova_vtl_tr_idx[16] = {
	0, 16, 32, 48, 4, 20, 36, 52, 8, 24, 40, 56, 12, 28, 44, 60
};
#endif

#if SNOVA_AVX2_RANK4
static inline __m128i snova_avx2_pack_row(__m128i b0, __m128i b1, __m128i b2, __m128i b3, int row) {
	const int base = row * 4;
	const __m128i m0 = _mm_setr_epi8((char)(base + 0), (char)(base + 1), (char)(base + 2), (char)(base + 3),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80);
	const __m128i m1 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(base + 0), (char)(base + 1), (char)(base + 2), (char)(base + 3),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80);
	const __m128i m2 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(base + 0), (char)(base + 1), (char)(base + 2), (char)(base + 3),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80);
	const __m128i m3 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(base + 0), (char)(base + 1), (char)(base + 2), (char)(base + 3));
	return _mm_or_si128(_mm_or_si128(_mm_shuffle_epi8(b0, m0), _mm_shuffle_epi8(b1, m1)),
	                    _mm_or_si128(_mm_shuffle_epi8(b2, m2), _mm_shuffle_epi8(b3, m3)));
}

static inline __m128i snova_avx2_pack_row_tr(__m128i b0, __m128i b1, __m128i b2, __m128i b3, int row) {
	const int base = row;
	const __m128i m0 = _mm_setr_epi8((char)(base + 0), (char)(base + 4), (char)(base + 8), (char)(base + 12),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80);
	const __m128i m1 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(base + 0), (char)(base + 4), (char)(base + 8), (char)(base + 12),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80);
	const __m128i m2 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(base + 0), (char)(base + 4), (char)(base + 8), (char)(base + 12),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80);
	const __m128i m3 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)0x80, (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(base + 0), (char)(base + 4), (char)(base + 8), (char)(base + 12));
	return _mm_or_si128(_mm_or_si128(_mm_shuffle_epi8(b0, m0), _mm_shuffle_epi8(b1, m1)),
	                    _mm_or_si128(_mm_shuffle_epi8(b2, m2), _mm_shuffle_epi8(b3, m3)));
}

static inline __m128i snova_avx2_transpose_row(__m128i r0, __m128i r1, __m128i r2, __m128i r3, int col) {
	const __m128i m0 = _mm_setr_epi8((char)(col + 0), (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(col + 4), (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(col + 8), (char)0x80, (char)0x80, (char)0x80,
	                                 (char)(col + 12), (char)0x80, (char)0x80, (char)0x80);
	const __m128i m1 = _mm_setr_epi8((char)0x80, (char)(col + 0), (char)0x80, (char)0x80,
	                                 (char)0x80, (char)(col + 4), (char)0x80, (char)0x80,
	                                 (char)0x80, (char)(col + 8), (char)0x80, (char)0x80,
	                                 (char)0x80, (char)(col + 12), (char)0x80, (char)0x80);
	const __m128i m2 = _mm_setr_epi8((char)0x80, (char)0x80, (char)(col + 0), (char)0x80,
	                                 (char)0x80, (char)0x80, (char)(col + 4), (char)0x80,
	                                 (char)0x80, (char)0x80, (char)(col + 8), (char)0x80,
	                                 (char)0x80, (char)0x80, (char)(col + 12), (char)0x80);
	const __m128i m3 = _mm_setr_epi8((char)0x80, (char)0x80, (char)0x80, (char)(col + 0),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)(col + 4),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)(col + 8),
	                                 (char)0x80, (char)0x80, (char)0x80, (char)(col + 12));
	return _mm_or_si128(_mm_or_si128(_mm_shuffle_epi8(r0, m0), _mm_shuffle_epi8(r1, m1)),
	                    _mm_or_si128(_mm_shuffle_epi8(r2, m2), _mm_shuffle_epi8(r3, m3)));
}

static inline __m128i snova_avx2_gather_pair(__m128i p0, __m128i p1, int nibble) {
	const int base = nibble;
	const __m128i m0 = _mm_setr_epi8((char)(base + 0), (char)(base + 0), (char)0x80, (char)0x80,
	                                 (char)(base + 4), (char)(base + 4), (char)0x80, (char)0x80,
	                                 (char)(base + 8), (char)(base + 8), (char)0x80, (char)0x80,
	                                 (char)(base + 12), (char)(base + 12), (char)0x80, (char)0x80);
	const __m128i m1 = _mm_setr_epi8((char)0x80, (char)0x80, (char)(base + 0), (char)(base + 0),
	                                 (char)0x80, (char)0x80, (char)(base + 4), (char)(base + 4),
	                                 (char)0x80, (char)0x80, (char)(base + 8), (char)(base + 8),
	                                 (char)0x80, (char)0x80, (char)(base + 12), (char)(base + 12));
	return _mm_or_si128(_mm_shuffle_epi8(p0, m0), _mm_shuffle_epi8(p1, m1));
}
#endif

static inline void snova_pack_block_run_ex(uint8_t* dst, int row_stride, const uint8_t* blocks, int block_stride,
                                            int count) {
	int b = 0;

#if SNOVA_NEON_RANK4
	const uint8x16_t idx0 = vld1q_u8(snova_run_idx);
	const uint8x16_t idx1 = vaddq_u8(idx0, vdupq_n_u8(4));
	const uint8x16_t idx2 = vaddq_u8(idx0, vdupq_n_u8(8));
	const uint8x16_t idx3 = vaddq_u8(idx0, vdupq_n_u8(12));

	for (; b + 4 <= count; b += 4) {
		uint8x16x4_t t;
		t.val[0] = vld1q_u8(blocks + (size_t)(b + 0) * block_stride);
		t.val[1] = vld1q_u8(blocks + (size_t)(b + 1) * block_stride);
		t.val[2] = vld1q_u8(blocks + (size_t)(b + 2) * block_stride);
		t.val[3] = vld1q_u8(blocks + (size_t)(b + 3) * block_stride);

		vst1q_u8(dst + 0 * row_stride + b * rank, vqtbl4q_u8(t, idx0));
		vst1q_u8(dst + 1 * row_stride + b * rank, vqtbl4q_u8(t, idx1));
		vst1q_u8(dst + 2 * row_stride + b * rank, vqtbl4q_u8(t, idx2));
		vst1q_u8(dst + 3 * row_stride + b * rank, vqtbl4q_u8(t, idx3));
	}
#elif SNOVA_AVX2_RANK4
	for (; b + 4 <= count; b += 4) {
		const __m128i b0 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 0) * block_stride));
		const __m128i b1 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 1) * block_stride));
		const __m128i b2 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 2) * block_stride));
		const __m128i b3 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 3) * block_stride));
		for (int row = 0; row < 4; ++row) {
			_mm_storeu_si128((__m128i*)(dst + row * row_stride + b * rank),
			                 snova_avx2_pack_row(b0, b1, b2, b3, row));
		}
	}
#endif

#if rank == 4

	for (; b < count; ++b) {
		const uint8_t* blk = blocks + (size_t)b * block_stride;
		uint8_t* base = dst + b * rank;

		memcpy(base + 0 * row_stride, blk + 0, 4);
		memcpy(base + 1 * row_stride, blk + 4, 4);
		memcpy(base + 2 * row_stride, blk + 8, 4);
		memcpy(base + 3 * row_stride, blk + 12, 4);
	}
#else
	for (; b < count; ++b) {
		for (int ei = 0; ei < rank; ++ei) {
			for (int ej = 0; ej < rank; ++ej) {
				dst[ei * row_stride + b * rank + ej] = blocks[(size_t)b * block_stride + ei * rank + ej];
			}
		}
	}
#endif
}

static inline void snova_pack_block_run_tr_ex(uint8_t* dst, int row_stride, const uint8_t* blocks, int block_stride,
                                               int count) {
	int b = 0;

#if SNOVA_NEON_RANK4
	const uint8x16_t idx0 = vld1q_u8(snova_run_tr_idx);
	const uint8x16_t idx1 = vaddq_u8(idx0, vdupq_n_u8(1));
	const uint8x16_t idx2 = vaddq_u8(idx0, vdupq_n_u8(2));
	const uint8x16_t idx3 = vaddq_u8(idx0, vdupq_n_u8(3));

	for (; b + 4 <= count; b += 4) {
		uint8x16x4_t t;
		t.val[0] = vld1q_u8(blocks + (size_t)(b + 0) * block_stride);
		t.val[1] = vld1q_u8(blocks + (size_t)(b + 1) * block_stride);
		t.val[2] = vld1q_u8(blocks + (size_t)(b + 2) * block_stride);
		t.val[3] = vld1q_u8(blocks + (size_t)(b + 3) * block_stride);

		vst1q_u8(dst + 0 * row_stride + b * rank, vqtbl4q_u8(t, idx0));
		vst1q_u8(dst + 1 * row_stride + b * rank, vqtbl4q_u8(t, idx1));
		vst1q_u8(dst + 2 * row_stride + b * rank, vqtbl4q_u8(t, idx2));
		vst1q_u8(dst + 3 * row_stride + b * rank, vqtbl4q_u8(t, idx3));
	}

	{
		const uint8x16_t tr = vld1q_u8(snova_tr_block_idx);
		for (; b < count; ++b) {
			uint32x4_t w = vreinterpretq_u32_u8(vqtbl1q_u8(vld1q_u8(blocks + (size_t)b * block_stride), tr));
			uint8_t* base = dst + b * rank;
			uint32_t w0 = vgetq_lane_u32(w, 0);
			uint32_t w1 = vgetq_lane_u32(w, 1);
			uint32_t w2 = vgetq_lane_u32(w, 2);
			uint32_t w3 = vgetq_lane_u32(w, 3);

			memcpy(base + 0 * row_stride, &w0, 4);
			memcpy(base + 1 * row_stride, &w1, 4);
			memcpy(base + 2 * row_stride, &w2, 4);
			memcpy(base + 3 * row_stride, &w3, 4);
		}
	}
#elif SNOVA_AVX2_RANK4
	for (; b + 4 <= count; b += 4) {
		const __m128i b0 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 0) * block_stride));
		const __m128i b1 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 1) * block_stride));
		const __m128i b2 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 2) * block_stride));
		const __m128i b3 = _mm_loadu_si128((const __m128i*)(blocks + (size_t)(b + 3) * block_stride));
		for (int row = 0; row < 4; ++row) {
			_mm_storeu_si128((__m128i*)(dst + row * row_stride + b * rank),
			                 snova_avx2_pack_row_tr(b0, b1, b2, b3, row));
		}
	}
	for (; b < count; ++b) {
		const __m128i block = _mm_loadu_si128((const __m128i*)(blocks + (size_t)b * block_stride));
		const __m128i tr = _mm_shuffle_epi8(block, _mm_loadu_si128((const __m128i*)snova_tr_block_idx));
		uint8_t* base = dst + b * rank;
		_mm_storeu_si32(base + 0 * row_stride, tr);
		_mm_storeu_si32(base + 1 * row_stride, _mm_srli_si128(tr, 4));
		_mm_storeu_si32(base + 2 * row_stride, _mm_srli_si128(tr, 8));
		_mm_storeu_si32(base + 3 * row_stride, _mm_srli_si128(tr, 12));
	}
#else
	for (; b < count; ++b) {
		for (int ei = 0; ei < rank; ++ei) {
			for (int ej = 0; ej < rank; ++ej) {
				dst[ei * row_stride + b * rank + ej] = blocks[(size_t)b * block_stride + ej * rank + ei];
			}
		}
	}
#endif
}

static inline void snova_pack_block_run(uint8_t* dst, int row_stride, const uint8_t* blocks, int count) {
	snova_pack_block_run_ex(dst, row_stride, blocks, sq_rank, count);
}

static inline void snova_pack_block_run_tr(uint8_t* dst, int row_stride, const uint8_t* blocks, int count) {
	snova_pack_block_run_tr_ex(dst, row_stride, blocks, sq_rank, count);
}

static inline void snova_pack_blocks(VECTOR* AJ, const uint8_t* A, int mainCol, int mainRow) {
	const int row_stride = vtl_mainRow_x_rank(mainRow);

	for (int mi = 0; mi < mainCol; ++mi) {
		snova_pack_block_run((uint8_t*)AJ + (mi * rank) * row_stride, row_stride, A + (mi * mainRow) * sq_rank, mainRow);
	}
}

static inline void snova_pack_blocks_tr(VECTOR* AJ, const uint8_t* A, int mainCol, int mainRow) {
	const int row_stride = vtl_mainRow_x_rank(mainRow);

	for (int mi = 0; mi < mainCol; ++mi) {
		snova_pack_block_run_tr((uint8_t*)AJ + (mi * rank) * row_stride, row_stride, A + (mi * mainRow) * sq_rank, mainRow);
	}
}

#if rank % 2
#define MJ_MAX mainRow
#else
#define MJ_MAX (vtl_mainRow_x_rank(mainRow)) / rank
#endif

static inline void snova_vtl_transpose(VECTOR* dst, const VECTOR* src, int mainCol, int mainRow) {
#if SNOVA_NEON_RANK4
	const int row_stride = vtl_mainRow_x_rank(mainRow);
	const uint8x16_t idx0 = vld1q_u8(snova_vtl_tr_idx);
	const uint8x16_t idx1 = vaddq_u8(idx0, vdupq_n_u8(1));
	const uint8x16_t idx2 = vaddq_u8(idx0, vdupq_n_u8(2));
	const uint8x16_t idx3 = vaddq_u8(idx0, vdupq_n_u8(3));

	for (int mi = 0; mi < mainCol; ++mi) {
		const uint8_t* s = (const uint8_t*)src + (mi * rank) * row_stride;
		uint8_t* d = (uint8_t*)dst + (mi * rank) * row_stride;

		for (int col = 0; col < row_stride; col += 16) {
			uint8x16x4_t t;
			t.val[0] = vld1q_u8(s + 0 * row_stride + col);
			t.val[1] = vld1q_u8(s + 1 * row_stride + col);
			t.val[2] = vld1q_u8(s + 2 * row_stride + col);
			t.val[3] = vld1q_u8(s + 3 * row_stride + col);

			vst1q_u8(d + 0 * row_stride + col, vqtbl4q_u8(t, idx0));
			vst1q_u8(d + 1 * row_stride + col, vqtbl4q_u8(t, idx1));
			vst1q_u8(d + 2 * row_stride + col, vqtbl4q_u8(t, idx2));
			vst1q_u8(d + 3 * row_stride + col, vqtbl4q_u8(t, idx3));
		}
	}
#elif SNOVA_AVX2_RANK4
	const int row_stride = vtl_mainRow_x_rank(mainRow);
	for (int mi = 0; mi < mainCol; ++mi) {
		const uint8_t* s = (const uint8_t*)src + (mi * rank) * row_stride;
		uint8_t* d = (uint8_t*)dst + (mi * rank) * row_stride;
		for (int col = 0; col < row_stride; col += 16) {
			const __m128i r0 = _mm_loadu_si128((const __m128i*)(s + 0 * row_stride + col));
			const __m128i r1 = _mm_loadu_si128((const __m128i*)(s + 1 * row_stride + col));
			const __m128i r2 = _mm_loadu_si128((const __m128i*)(s + 2 * row_stride + col));
			const __m128i r3 = _mm_loadu_si128((const __m128i*)(s + 3 * row_stride + col));
			const __m128i d0 = snova_avx2_transpose_row(r0, r1, r2, r3, 0);
			const __m128i d1 = snova_avx2_transpose_row(r0, r1, r2, r3, 1);
			const __m128i d2 = snova_avx2_transpose_row(r0, r1, r2, r3, 2);
			const __m128i d3 = snova_avx2_transpose_row(r0, r1, r2, r3, 3);
			_mm_storeu_si128((__m128i*)(d + 0 * row_stride + col), d0);
			_mm_storeu_si128((__m128i*)(d + 1 * row_stride + col), d1);
			_mm_storeu_si128((__m128i*)(d + 2 * row_stride + col), d2);
			_mm_storeu_si128((__m128i*)(d + 3 * row_stride + col), d3);
		}
	}
#else
	for (int mi = 0; mi < mainCol; ++mi) {
		for (int mj = 0; mj < MJ_MAX; ++mj) {
			for (int ei = 0; ei < rank; ++ei) {
				for (int ej = 0; ej < rank; ++ej) {
					VEC_BYTE(dst, ((mi * rank + ei) * vtl_mainRow_x_rank(mainRow)) + (mj * rank + ej)) =
					    VEC_BYTE(src, ((mi * rank + ej) * vtl_mainRow_x_rank(mainRow)) + (mj * rank + ei));
				}
			}
		}
	}
#endif
}

#if rank == 4

static inline void snova_pack_rankpair_tr(VECTOR* restrict dst_pair0, VECTOR* restrict dst_pair1,
                                           const VECTOR* restrict src, int mainRow) {
	const int row_stride = vtl_mainRow_x_rank(mainRow);
	const uint8_t* restrict src01 = (const uint8_t*)src;
	const uint8_t* restrict src23 = (const uint8_t*)src + 2 * row_stride;
	uint8_t* restrict out0 = (uint8_t*)dst_pair0;
	uint8_t* restrict out1 = (uint8_t*)dst_pair1;

#if SNOVA_NEON_RANK4
	static const uint8_t gather_idx[16] = {
		0, 0, 16, 16, 4, 4, 20, 20, 8, 8, 24, 24, 12, 12, 28, 28
	};
	static const uint8_t even_bytes[16] = {
		0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00,
		0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00
	};

	const uint8x16_t idx_x0 = vld1q_u8(gather_idx);
	const uint8x16_t idx_y0 = vaddq_u8(idx_x0, vdupq_n_u8(1));
	const uint8x16_t idx_x1 = vaddq_u8(idx_x0, vdupq_n_u8(2));
	const uint8x16_t idx_y1 = vaddq_u8(idx_x0, vdupq_n_u8(3));
	const uint8x16_t even_mask = vld1q_u8(even_bytes);
	const uint8x16_t low_mask = vdupq_n_u8(0x0f);
	const uint8x16_t high_mask = vdupq_n_u8(0xf0);

	for (int col = 0; col < row_stride; col += 16) {
		uint8x16x2_t tbl;
		tbl.val[0] = vld1q_u8(src01 + col);
		tbl.val[1] = vld1q_u8(src23 + col);

		uint8x16_t x0 = vqtbl2q_u8(tbl, idx_x0);
		uint8x16_t y0 = vqtbl2q_u8(tbl, idx_y0);
		uint8x16_t x1 = vqtbl2q_u8(tbl, idx_x1);
		uint8x16_t y1 = vqtbl2q_u8(tbl, idx_y1);

		uint8x16_t lo0 = vorrq_u8(vandq_u8(x0, low_mask), vshlq_n_u8(vandq_u8(y0, low_mask), 4));
		uint8x16_t hi0 = vorrq_u8(vshrq_n_u8(x0, 4), vandq_u8(y0, high_mask));
		uint8x16_t lo1 = vorrq_u8(vandq_u8(x1, low_mask), vshlq_n_u8(vandq_u8(y1, low_mask), 4));
		uint8x16_t hi1 = vorrq_u8(vshrq_n_u8(x1, 4), vandq_u8(y1, high_mask));

		vst1q_u8(out0 + col, vbslq_u8(even_mask, lo0, hi0));
		vst1q_u8(out1 + col, vbslq_u8(even_mask, lo1, hi1));
	}
#elif SNOVA_AVX2_RANK4
	const __m128i low_mask = _mm_set1_epi8(0x0f);
	const __m128i high_mask = _mm_set1_epi8((char)0xf0);
	const __m128i even_mask = _mm_setr_epi8((char)0xff, 0, (char)0xff, 0, (char)0xff, 0, (char)0xff, 0,
	                                        (char)0xff, 0, (char)0xff, 0, (char)0xff, 0, (char)0xff, 0);
	for (int col = 0; col < row_stride; col += 16) {
		const __m128i p0 = _mm_loadu_si128((const __m128i*)(src01 + col));
		const __m128i p1 = _mm_loadu_si128((const __m128i*)(src23 + col));
		const __m128i x0 = snova_avx2_gather_pair(p0, p1, 0);
		const __m128i y0 = snova_avx2_gather_pair(p0, p1, 1);
		const __m128i x1 = snova_avx2_gather_pair(p0, p1, 2);
		const __m128i y1 = snova_avx2_gather_pair(p0, p1, 3);
		const __m128i lo0 = _mm_or_si128(_mm_and_si128(x0, low_mask), _mm_slli_epi16(_mm_and_si128(y0, low_mask), 4));
		const __m128i hi0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(x0, 4), low_mask), _mm_and_si128(y0, high_mask));
		const __m128i lo1 = _mm_or_si128(_mm_and_si128(x1, low_mask), _mm_slli_epi16(_mm_and_si128(y1, low_mask), 4));
		const __m128i hi1 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(x1, 4), low_mask), _mm_and_si128(y1, high_mask));
		_mm_storeu_si128((__m128i*)(out0 + col), _mm_or_si128(_mm_and_si128(even_mask, lo0), _mm_andnot_si128(even_mask, hi0)));
		_mm_storeu_si128((__m128i*)(out1 + col), _mm_or_si128(_mm_and_si128(even_mask, lo1), _mm_andnot_si128(even_mask, hi1)));
	}
#else
	for (int col = 0; col < row_stride; col += rank) {
		const uint8_t a0 = src01[col + 0];
		const uint8_t a1 = src01[col + 1];
		const uint8_t a2 = src01[col + 2];
		const uint8_t a3 = src01[col + 3];
		const uint8_t b0 = src23[col + 0];
		const uint8_t b1 = src23[col + 1];
		const uint8_t b2 = src23[col + 2];
		const uint8_t b3 = src23[col + 3];

		out0[col + 0] = (a0 & 0x0f) ^ ((a1 & 0x0f) << 4);
		out0[col + 1] = (a0 >> 4) ^ (a1 & 0xf0);
		out0[col + 2] = (b0 & 0x0f) ^ ((b1 & 0x0f) << 4);
		out0[col + 3] = (b0 >> 4) ^ (b1 & 0xf0);

		out1[col + 0] = (a2 & 0x0f) ^ ((a3 & 0x0f) << 4);
		out1[col + 1] = (a2 >> 4) ^ (a3 & 0xf0);
		out1[col + 2] = (b2 & 0x0f) ^ ((b3 & 0x0f) << 4);
		out1[col + 3] = (b2 >> 4) ^ (b3 & 0xf0);
	}
#endif
}
#endif

static inline void snova_pack_rankpair_tr_rows(uint8_t* dst, const VECTOR* src_pairs, int mainRow) {
	const int row_stride = vtl_mainRow_x_rank(mainRow);

#if rank == 4
	snova_pack_rankpair_tr((VECTOR*)dst, (VECTOR*)(dst + row_stride), src_pairs, mainRow);
#else

	const int max_stride = vtl_mainRow_x_rank(SNOVA_n);
	alignas(VEC_LENGTH) uint8_t split[rank_next2 * vtl_mainRow_x_rank(SNOVA_n)] = {0};
	alignas(VEC_LENGTH) uint8_t tr[rank_next2 * vtl_mainRow_x_rank(SNOVA_n)] = {0};
	const uint8_t* src = (const uint8_t*)src_pairs;

	(void)max_stride;

	for (int p = 0; p < rank_next32; ++p) {
		for (int c = 0; c < row_stride; ++c) {
			uint8_t packed = src[(2 * p) * row_stride + c];
			split[(2 * p) * row_stride + c] = packed & 0x0f;
			split[(2 * p + 1) * row_stride + c] = packed >> 4;
		}
	}

	for (int ei = 0; ei < rank; ++ei) {
		for (int mj = 0; mj < row_stride / rank; ++mj) {
			for (int ej = 0; ej < rank; ++ej) {
				tr[ei * row_stride + mj * rank + ej] = split[ej * row_stride + mj * rank + ei];
			}
		}
	}

	for (int p = 0; p < rank_next32; ++p) {
		for (int c = 0; c < row_stride; ++c) {
			uint8_t lo = tr[(2 * p) * row_stride + c];
			uint8_t hi = (2 * p + 1 < rank) ? tr[(2 * p + 1) * row_stride + c] : 0;
			dst[p * row_stride + c] = (uint8_t)(lo | (hi << 4));
		}
	}
#endif
}

#endif
