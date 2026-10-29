/* SPDX-License-Identifier: MIT */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "snova.h"
#include "gf16_r3_rvv.h"

static uint32_t state = 0x6d2b79f5u;

static uint16_t next_value(void) {
	state = state * 1664525u + 1013904223u;
	return (uint16_t)(state & 0xffffu);
}

static void reference_axpy(uint16_t *dst, const uint16_t *src, uint16_t scalar, size_t n) {
	for (size_t i = 0; i < n; ++i)
		dst[i] ^= (uint16_t)((uint32_t)scalar * (uint32_t)src[i]);
}

static int check_axpy(void) {
	static const size_t lengths[] = {1, 4, 16, 31, 96};
	uint16_t expected[96], actual[96], src[96];
	for (size_t c = 0; c < sizeof(lengths) / sizeof(lengths[0]); ++c) {
		size_t n = lengths[c];
		for (size_t i = 0; i < n; ++i) {
			expected[i] = actual[i] = next_value();
			src[i] = next_value();
		}
		uint16_t scalar = next_value();
		reference_axpy(expected, src, scalar, n);
		snova_r3_gf16_axpy(actual, src, scalar, n);
		if (memcmp(expected, actual, n * sizeof(uint16_t)) != 0)
			return -1;
	}
	return 0;
}

static int check_axpy_reuse_vl(void) {
	static const size_t lengths[] = {1, 4, 16, 31, 96};
	uint16_t expected[96], actual[96], src[96];
	for (size_t c = 0; c < sizeof(lengths) / sizeof(lengths[0]); ++c) {
		size_t n = lengths[c];
		for (size_t i = 0; i < n; ++i) {
			expected[i] = actual[i] = next_value();
			src[i] = next_value();
		}
		uint16_t scalar = next_value();
		size_t vl = snova_r3_gf16_setvl(n);
		reference_axpy(expected, src, scalar, n);
		snova_r3_gf16_axpy_vl(actual, src, scalar, n, vl);
		if (memcmp(expected, actual, n * sizeof(uint16_t)) != 0)
			return -1;
	}
	return 0;
}

static int check_compress_and_normalize(void) {
	uint16_t expected[96], actual[96];
	for (size_t i = 0; i < 96; ++i)
		expected[i] = actual[i] = next_value();
	for (size_t i = 0; i < 96; ++i)
		expected[i] = snova_r3_gf16_compress_scalar(expected[i]);
	snova_r3_gf16_compress(actual, 96);
	if (memcmp(expected, actual, sizeof(expected)) != 0)
		return -1;

	for (size_t i = 0; i < 96; ++i)
		expected[i] = actual[i] = next_value();
	uint16_t scalar = (uint16_t)(next_value() & 0x1111u);
	for (size_t i = 0; i < 96; ++i)
		expected[i] = snova_r3_gf16_compress_scalar((uint16_t)(expected[i] * scalar));
	snova_r3_gf16_mul_compress(actual, scalar, 96);
	return memcmp(expected, actual, sizeof(expected)) == 0 ? 0 : -1;
}

static int check_accumulate_terms(size_t terms, size_t count) {
	static uint16_t sources[SNOVA_v * SNOVA_l][96];
	const uint16_t *source_ptrs[SNOVA_v * SNOVA_l];
	uint16_t scalars[SNOVA_v * SNOVA_l];
	uint16_t expected[96], axpy_expected[96], actual[96];
	for (size_t term = 0; term < terms; ++term) {
		source_ptrs[term] = sources[term];
		scalars[term] = next_value();
		for (size_t i = 0; i < count; ++i)
			sources[term][i] = next_value();
	}
	for (size_t i = 0; i < count; ++i)
		expected[i] = axpy_expected[i] = actual[i] = next_value();
	for (size_t term = 0; term < terms; ++term)
		for (size_t i = 0; i < count; ++i)
			expected[i] ^= (uint16_t)((uint32_t)scalars[term] *
				(uint32_t)sources[term][i]);
	for (size_t term = 0; term < terms; ++term)
		snova_r3_gf16_axpy(axpy_expected, source_ptrs[term], scalars[term], count);
	snova_r3_gf16_accumulate_terms(actual, source_ptrs, scalars, terms, count);
	for (size_t i = 0; i < count; ++i) {
		if (expected[i] != axpy_expected[i] || expected[i] != actual[i]) {
			printf("accum mismatch terms=%zu count=%zu index=%zu\n",
				terms, count, i);
			printf("  scalar=0x%04x repeated_axpy=0x%04x fused=0x%04x\n",
				expected[i], axpy_expected[i], actual[i]);
			return -1;
		}
	}
	return 0;
}

static int check_mask(void) {
	uint16_t expected[96], actual[96], src[96];
	for (size_t i = 0; i < 96; ++i) {
		expected[i] = actual[i] = next_value();
		src[i] = next_value();
	}
	uint16_t mask = 0xffffu;
	for (size_t i = 0; i < 96; ++i)
		expected[i] ^= (uint16_t)(mask & src[i]);
	snova_r3_xor_mask(actual, src, mask, 96);
	return memcmp(expected, actual, sizeof(expected)) == 0 ? 0 : -1;
}

int main(void) {
	int axpy = check_axpy();
	int axpy_vl = check_axpy_reuse_vl();
	int compress = check_compress_and_normalize();
	int accum_small = check_accumulate_terms(4, 96);
	int accum_sum_t0 = check_accumulate_terms(SNOVA_v * SNOVA_l, SNOVA_lr);
	int mask = check_mask();
	int result = axpy | axpy_vl | compress | accum_small | accum_sum_t0 | mask;
	printf("kernel-test axpy=%s axpy-reuse-vl=%s\n",
		axpy ? "FAIL" : "PASS", axpy_vl ? "FAIL" : "PASS");
	printf("kernel-test compress/mul=%s accum-4x96=%s\n",
		compress ? "FAIL" : "PASS",
		accum_small ? "FAIL" : "PASS");
	printf("kernel-test sum_t0_accum_108x16=%s mask=%s\n",
		accum_sum_t0 ? "FAIL" : "PASS", mask ? "FAIL" : "PASS");
	printf("Round3 RVV kernel self-test: %s\n", result ? "FAIL" : "PASS");
	return result ? 1 : 0;
}
