#include "snova_transform.h"

int snova_bytes_to_gf16(const uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s) {
	const int pairs = num_of_GF16s >> 1;
	int i = 0;

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	const uint8x16_t low_mask = vdupq_n_u8(0x0f);

	for (; i + 16 <= pairs; i += 16) {
		uint8x16_t packed = vld1q_u8(byte_array + i);
		uint8x16x2_t split;

		split.val[0] = vandq_u8(packed, low_mask);
		split.val[1] = vshrq_n_u8(packed, 4);
		vst2q_u8(gf16_array + 2 * i, split);
	}
#endif

	for (; i < pairs; ++i) {
		gf16_array[i * 2] = byte_array[i] & 0x0f;
		gf16_array[i * 2 + 1] = (byte_array[i] >> 4) & 0x0f;
	}

	if (num_of_GF16s & 1) {
		gf16_array[num_of_GF16s - 1] = byte_array[pairs] & 0x0f;
		return byte_array[pairs] & 0xf0;
	}

	return 0;
}

void snova_gf16_to_bytes(uint8_t* byte_array, const uint8_t* gf16_array, int num_of_GF16s) {
	const int pairs = num_of_GF16s >> 1;
	int i = 0;

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
	for (; i + 16 <= pairs; i += 16) {
		uint8x16x2_t split = vld2q_u8(gf16_array + 2 * i);
		vst1q_u8(byte_array + i, vorrq_u8(split.val[0], vshlq_n_u8(split.val[1], 4)));
	}
#endif

	for (; i < pairs; ++i) {
		byte_array[i] = gf16_array[i * 2] | (gf16_array[i * 2 + 1] << 4);
	}

	if (num_of_GF16s & 1) {
		byte_array[pairs] = gf16_array[num_of_GF16s - 1];
	}
}

void snova_bytes_to_gf16_split(const uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s) {
	const int half_GF16s = (num_of_GF16s + 1) >> 1;

	for (int i = 0; i < half_GF16s; ++i) {
		gf16_array[i] = byte_array[i] & 0x0f;
	}

	for (int i = 0; i < (num_of_GF16s >> 1); ++i) {
		gf16_array[i + half_GF16s] = byte_array[i] >> 4;
	}
}

void snova_gf16_to_bytes_merge(uint8_t* byte_array, const uint8_t* gf16_array, int num_of_GF16s) {
	const int half_GF16s = (num_of_GF16s + 1) >> 1;
	int i;

	for (i = 0; i < (num_of_GF16s >> 1); ++i) {
		byte_array[i] = gf16_array[i] | (gf16_array[i + half_GF16s] << 4);
	}

	if (num_of_GF16s & 1) {
		byte_array[i] = gf16_array[i];
	}
}
