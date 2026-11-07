#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "timing.h"
#include "utils.h"
#include "api.h"
#include "benchmark.h"

#define B_KEY_GENERATION 0
#define B_SIGN_ALGO 1
#define B_VERIFY_ALGO 2
#define NUMBER_OF_ALGO_BENCHES 3

int randombytes(unsigned char* x, unsigned long long xlen) {
	for (unsigned long long j = 0; j < xlen; j++) {
		x[j] = (uint8_t) rand();
	}
	return 0;
}

#ifdef BENCHMARK
btimer_t timers[NUMBER_OF_BENCHES];

#ifdef BENCHMARK_CYCLES
#define display_timer(label,num) printf("   - " label ": %f ms (%f cycles)\n", btimer_get(&timers[num]), btimer_get_cycles(&timers[num]))
#else
#define display_timer(label,num) printf("   - " label ": %f ms\n", btimer_get(&timers[num]))
#endif
#endif

int main(int argc, char *argv[]) {
	srand((unsigned int) time(NULL));

	int nb_tests = get_number_of_tests(argc, argv, 1);
	if (nb_tests < 0) {
		exit(EXIT_FAILURE);
	}

	print_configuration();
	printf("\n");

	btimer_t timers_algos[NUMBER_OF_ALGO_BENCHES];
	double std_timer[NUMBER_OF_ALGO_BENCHES];

	// Initialisation
	double timer_pow2[NUMBER_OF_ALGO_BENCHES];
	for (int j = 0; j < NUMBER_OF_ALGO_BENCHES; j++) {
		btimer_init(&timers_algos[j]);
		timer_pow2[j] = 0;
	}
#ifdef BENCHMARK
	for (int num = 0; num < NUMBER_OF_BENCHES; num++) {
		btimer_init(&timers[num]);
	}
#endif
	double mean_of_sig_size = 0;
	double sig_size_pow2 = 0;

	// Execution
	int score = 0;
	int ret;
	for (int i = 0; i < nb_tests; i++) {
#ifdef BENCHMARK
		for (int num = 0; num < NUMBER_OF_BENCHES; num++) {
			btimer_count(&timers[num]);
		}
#endif

		// Generate the keys
		uint8_t pk[CRYPTO_PUBLICKEYBYTES];
		uint8_t sk[CRYPTO_SECRETKEYBYTES];
		btimer_start(&timers_algos[B_KEY_GENERATION]);
		ret = crypto_sign_keypair(pk, sk);
		btimer_end(&timers_algos[B_KEY_GENERATION]);
		btimer_count(&timers_algos[B_KEY_GENERATION]);
		timer_pow2[B_KEY_GENERATION] += pow(btimer_diff(&timers_algos[B_KEY_GENERATION]), 2) / nb_tests;
		if (ret) {
			printf("Failure (num %d): crypto_sign_keypair\n", i);
			continue;
		}

		// Select the message
#define MLEN 32
		uint8_t m[MLEN] = {1, 2, 3, 4};
		uint8_t m2[MLEN] = {0};
		unsigned long long m2len = 0;

		// Sign the message
		uint8_t sm[MLEN + CRYPTO_BYTES];
		unsigned long long smlen = 0;
		btimer_start(&timers_algos[B_SIGN_ALGO]);
		ret = crypto_sign(sm, &smlen, m, MLEN, sk);
		btimer_end(&timers_algos[B_SIGN_ALGO]);
		btimer_count(&timers_algos[B_SIGN_ALGO]);
		timer_pow2[B_SIGN_ALGO] += pow(btimer_diff(&timers_algos[B_SIGN_ALGO]), 2) / nb_tests;
		// Update statistics
		size_t signature_len = smlen - MLEN;
		mean_of_sig_size += (double) signature_len / nb_tests;
		sig_size_pow2 += pow(signature_len, 2) / nb_tests;
		if (ret) {
			printf("Failure (num %d): crypto_sign\n", i);
			continue;
		}

		// Verify/Open the signature
		btimer_start(&timers_algos[B_VERIFY_ALGO]);
		ret = crypto_sign_open(m2, &m2len, sm, smlen, pk);
		btimer_end(&timers_algos[B_VERIFY_ALGO]);
		btimer_count(&timers_algos[B_VERIFY_ALGO]);
		timer_pow2[B_VERIFY_ALGO] += pow(btimer_diff(&timers_algos[B_VERIFY_ALGO]), 2) / nb_tests;
		if (ret) {
			printf("Failure (num %d): crypto_sign_open\n", i);
			continue;
		}

		// Test of correction of the primitives
		if (m2len != MLEN) {
			printf("Failure (num %d): message size does not match\n", i);
			continue;
		}
		for (int h = 0; h < MLEN; h++)
			if (m[h] != m2[h]) {
				printf("Failure (num %d): message does not match (char %d)\n", i, h);
				continue;
			}

		score++;
	}

	// Compute some statistics
	std_timer[B_KEY_GENERATION] = sqrt(timer_pow2[B_KEY_GENERATION] - pow(btimer_get(&timers_algos[B_KEY_GENERATION]), 2));
	std_timer[B_SIGN_ALGO] = sqrt(timer_pow2[B_SIGN_ALGO] - pow(btimer_get(&timers_algos[B_SIGN_ALGO]), 2));
	std_timer[B_VERIFY_ALGO] = sqrt(timer_pow2[B_VERIFY_ALGO] - pow(btimer_get(&timers_algos[B_VERIFY_ALGO]), 2));
	double std_sig_size = sqrt(sig_size_pow2 - pow(mean_of_sig_size, 2));

	// Display Infos
	printf("===== SUMMARY =====\n");
	printf("Correctness: %d/%d\n", score, nb_tests);
	printf("\n");

	print_alloc_usage("keygen+sign+verif");
	printf("\n");

	printf("Timing in ms:\n");
	printf(" - Key Gen: %.2f ms (std=%.2f)\n",
	       btimer_get(&timers_algos[B_KEY_GENERATION]),
	       std_timer[B_KEY_GENERATION]
	      );
	printf(" - Sign:    %.2f ms (std=%.2f)\n",
	       btimer_get(&timers_algos[B_SIGN_ALGO]),
	       std_timer[B_SIGN_ALGO]
	      );
	printf(" - Verify:  %.2f ms (std=%.2f)\n",
	       btimer_get(&timers_algos[B_VERIFY_ALGO]),
	       std_timer[B_VERIFY_ALGO]
	      );
	printf("\n");

#ifdef BENCHMARK_CYCLES
	printf("Timing in cycles:\n");
	printf(" - Key Gen: %.2f cycles\n", btimer_get_cycles(&timers_algos[B_KEY_GENERATION]));
	printf(" - Sign:    %.2f cycles\n", btimer_get_cycles(&timers_algos[B_SIGN_ALGO]));
	printf(" - Verify:  %.2f cycles\n", btimer_get_cycles(&timers_algos[B_VERIFY_ALGO]));
	printf("\n");
#endif

	printf("Communication cost:\n");
	printf(" - PK size: %ld B\n", CRYPTO_PUBLICKEYBYTES);
	printf(" - SK size: %ld B\n", CRYPTO_SECRETKEYBYTES);
	printf(" - Signature size (MAX): %ld B\n", CRYPTO_BYTES);
	printf(" - Signature size: %.0f B (std=%.0f)\n", mean_of_sig_size, std_sig_size);
	printf("\n");

#ifdef BENCHMARK
	printf("\n===== DETAILED BENCHMARK =====\n");
	printf(" - Signing\n");
	display_timer("BLC.Commit", BS_BLC_COMMIT);
	display_timer("[BLC.Commit] Expand Trees", BS_BLC_EXPAND_TREE);
	display_timer("[BLC.Commit] Seed Commit", BS_BLC_SEED_COMMIT);
	display_timer("[BLC.Commit] PRG", BS_BLC_PRG);
	display_timer("[BLC.Commit] XOF", BS_BLC_XOF);
	display_timer("[BLC.Commit] Arithm", BS_BLC_ARITH);
	display_timer("PIOP.Compute", BS_PIOP_COMPUTE);
	display_timer("[PIOP.Compute] ExpandMQ", BS_PIOP_EXPAND_MQ);
	display_timer("[PIOP.Compute] Expand Batching Mat", BS_PIOP_EXPAND_BATCHING_MAT);
	display_timer("[PIOP.Compute] Matrix Mul Ext", BS_PIOP_MAT_MUL_EXT);
	display_timer("[PIOP.Compute] Compute t1", BS_PIOP_COMPUTE_T1);
	display_timer("[PIOP.Compute] Compute P_zi", BS_PIOP_COMPUTE_PZI);
	display_timer("[PIOP.Compute] Batch and Mask", BS_PIOP_BATCH_AND_MASK);
	display_timer("Sample Challenge", BS_SAMPLE_CHALLENGE);
	display_timer("BLC.Open", BS_BLC_OPEN);
	printf(" - Others\n");
	display_timer("Pin A", B_PIN_A);
	display_timer("Pin B", B_PIN_B);
	display_timer("Pin C", B_PIN_C);
	display_timer("Pin D", B_PIN_D);
#endif

	return 0;
}
