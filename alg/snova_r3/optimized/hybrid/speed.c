// SPDX-License-Identifier: MIT

/**
 * Test program to determine scheme cycle counts.
 * Also serves as a test to verify behavior for both legitimate and forged signatures.
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "api.h"
#include "rng.h"

// Max test duration in seconds
#ifndef DURATION
#define DURATION 3
#endif

// Max number of tests
#ifndef MAX_TESTS
#define MAX_TESTS 2048
#endif

#define TEXT_LEN 64
#define BYTES_DIGEST 64

/**
 * To get detailed API timing info use
 *
 * make clean speed P="-D DETAILS"
 *
 * make clean speed P="-D ANALYZE_TIMING -D AVERAGE"
 */

uint64_t cycles0 = 0, cycles1 = 0, cycles2 = 0, cycles3 = 0, cycles4 = 0, cycles5 = 0, cycles6 = 0, cycles7 = 0;

double get_cpu_f(void);

#ifdef __ARM_ARCH
// cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_cur_freq
// sudo cpufreq-set -d 2.4GHz
#define CYCLES_PER_NS 2.4
#else
// Uncomment for times in ns
// #define CYCLES_PER_NS 1
#endif

#ifdef CYCLES_PER_NS
uint64_t get_cycles(void) {
	struct timespec time = {0};
	timespec_get(&time, TIME_UTC);
	return (int64_t)((time.tv_sec * 1e9 + time.tv_nsec) * CYCLES_PER_NS);
}
#else
uint64_t get_cycles(void) {
	uint32_t lo = 0, hi = 0;
	uint64_t o = 0;
	__asm__ __volatile__("rdtscp" : "=a"(lo), "=d"(hi) : : "%ecx");
	o = hi;
	o <<= 32;
	return (o | lo);
}
#endif

void print_number(unsigned long n) {
	if (n < 1000) {
		printf("%lu", n);
		return;
	}
	print_number(n / 1000);
	printf(",%03lu", n % 1000);
}

int cmp_uint64(const void *a, const void *b) {
	if (*(uint64_t * )a < * (uint64_t * )b) {
		return -1;
	}
	if (*(uint64_t * )a > *(uint64_t * )b) {
		return 1;
	}
	return 0;
}

uint64_t median(uint64_t *l, size_t len) {
	qsort(l, len, sizeof(uint64_t), cmp_uint64);

	if (len % 2) {
		return l[len / 2];
	} else {
		return (l[len / 2 - 1] + l[len / 2]) / 2;
	}
}

uint64_t summary(uint64_t *t, int test_n) {
	uint64_t td[MAX_TESTS];

	for (int i = 0; i < test_n; ++i) {
		td[i] = t[i * 2 + 1] - t[i * 2];
	}

	uint64_t res = median(td, test_n);

	return res;
}

uint64_t average(uint64_t *t, size_t len) {
	uint64_t acc = 0;
	for (size_t i = 0; i < len; i++) {
		acc += t[i];
	}
	return acc / len;
}

uint64_t avg_summary(uint64_t *t, int test_n) {
	uint64_t td[MAX_TESTS];

	for (int i = 0; i < test_n; ++i) {
		td[i] = t[i * 2 + 1] - t[i * 2];
	}

	uint64_t res = average(td, test_n);

	return res;
}

#ifdef AVERAGE
#define SUMMARY avg_summary
#else
#define SUMMARY summary
#endif

int main(void) {
	uint8_t seed[48] = {0};
	uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
	uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
	uint8_t sm[CRYPTO_BYTES + TEXT_LEN] = {0};

#ifndef DETAILS
	uint8_t text[TEXT_LEN] = {0};
	uint8_t text1[TEXT_LEN] = {0};
	unsigned long long smlen = 0;
	unsigned long long len1 = 0;
	int res = 0;
	int fail = 0;
	uint64_t t0[MAX_TESTS * 2] = {0};
	uint64_t t1[MAX_TESTS * 2] = {0};
	uint64_t t2[MAX_TESTS * 2] = {0};
#else
	uint64_t t0[MAX_TESTS * 2] = {0};
	uint64_t t1[MAX_TESTS * 2] = {0};
	uint64_t t2[MAX_TESTS * 2] = {0};
	uint64_t t3[MAX_TESTS * 2] = {0};
	uint64_t t4[MAX_TESTS * 2] = {0};
#endif
	int r = 0;

	randombytes_init(seed, NULL, 256);

	struct timespec start_time = {0};
	timespec_get(&start_time, TIME_UTC);
	int i = 0;
	for (; i < MAX_TESTS; i++) {
#ifndef DETAILS
		t0[i * 2] = get_cycles();
		res += crypto_sign_keypair(pk, sk) & 1;
		t0[i * 2 + 1] = get_cycles();

		t1[i * 2] = get_cycles();
		res += crypto_sign(sm, &smlen, text, TEXT_LEN, sk) & 1;
		t1[i * 2 + 1] = get_cycles();

		// Also test failures
		fail = i % 2;
		if (fail) {
			sm[(i / 8) % CRYPTO_BYTES] ^= 1 << (i % 8);
		}
		t2[i * 2] = get_cycles();
		res = crypto_sign_open(text1, &len1, sm, smlen, pk);
		t2[i * 2 + 1] = get_cycles();
		r += (res ^ fail) & 1;
#else
		uint8_t bseed[SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE];
		uint8_t salt[BYTES_SALT];
		uint8_t digest[BYTES_DIGEST] = {0};
		expanded_SK skx_d = {0};
		expanded_PK pkx = {0};

		randombytes(bseed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);
		randombytes(salt, BYTES_SALT);

		digest[0] = i & 0xff;
		digest[1] = (i >> 8) & 0xff;
		digest[2] = (i >> 16) & 0xff;

		t0[i * 2] = get_cycles();
		r += SNOVA_NAMESPACE(genkeys)(pk, sk, bseed);
		t0[i * 2 + 1] = get_cycles();

		t1[i * 2] = get_cycles();
		r += SNOVA_NAMESPACE(sk_expand)(&skx_d, sk);
		t1[i * 2 + 1] = get_cycles();

		t2[i * 2] = get_cycles();
		r += SNOVA_NAMESPACE(sign)(&skx_d, sm, digest, BYTES_DIGEST, salt);
		t2[i * 2 + 1] = get_cycles();

		t3[i * 2] = get_cycles();
		r += SNOVA_NAMESPACE(pk_expand)(&pkx, pk);
		t3[i * 2 + 1] = get_cycles();

		t4[i * 2] = get_cycles();
		r += SNOVA_NAMESPACE(verify)(&pkx, sm, digest, BYTES_DIGEST);
		t4[i * 2 + 1] = get_cycles();
#endif

		struct timespec time;
		timespec_get(&time, TIME_UTC);
		if (((time.tv_sec - start_time.tv_sec) + (time.tv_nsec - start_time.tv_nsec) / 1e9) > DURATION) {
			break;
		}
	}
	int test_n = i;

	if (r == 0) {
		if (test_n != 2048) {
			printf("! ");
		}
#ifdef DETAILS
		printf("%s", CRYPTO_ALGNAME);
#else
		printf("%s & %d & %d & %d", CRYPTO_ALGNAME, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, CRYPTO_BYTES);
#endif
	} else {
		printf("\nFAIL!!\n\n%s & %d & %d & %d", CRYPTO_ALGNAME, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, CRYPTO_BYTES);
	}

#ifdef AVERAGE
	printf(" & avg.");
#endif

	printf(" & ");
	print_number(SUMMARY(t0, test_n));
	printf(" & ");
	print_number(SUMMARY(t1, test_n));
	printf(" & ");
	print_number(SUMMARY(t2, test_n));

#ifdef DETAILS
	printf(" & ");
	print_number(SUMMARY(t3, test_n));
	printf(" & ");
	print_number(SUMMARY(t4, test_n));
#endif

	if (cycles0) {
		printf(" &      ");
		print_number(cycles0 / test_n);
		printf("  &  ");
		print_number(cycles1 / test_n);
		printf("  &  ");
		print_number(cycles2 / test_n);
		printf("  &  ");
		print_number(cycles3 / test_n);
		printf("  &  ");
		print_number(cycles4 / test_n);
		printf("  &  ");
		print_number(cycles5 / test_n);
		printf("  &  ");
		print_number(cycles6 / test_n);
		printf("  &  ");
		print_number(cycles7 / test_n);
	}

	printf(" \\\\\n");

	return 0;
}
