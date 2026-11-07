#ifndef MQOM_BENCH_UTILS_H
#define MQOM_BENCH_UTILS_H

#include <stdio.h>
#include "api.h"
#include "rijndael/rijndael.h"
#include "fields.h"
#include "fields_bitsliced.h"
#include "blc.h"

#ifndef STR
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#endif

static inline int get_number_of_tests(int argc, char *argv[], int default_value) {
	int nb_tests = default_value;
	if (argc == 2) {
		if ( sscanf(argv[1], "%d", &nb_tests) != 1) {
			printf("Integer awaited.\n");
			return -1;
		} else if ( nb_tests <= 0 ) {
			printf("Need to positive integer.\n");
			return -1;
		}
	}
	return nb_tests;
}

/* Architecture detection */
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64) || defined(_M_AMD64)
#  define COMPILER_PLATFORM "x86_64"
#elif defined(__i386__) || defined(_M_IX86)
#  define COMPILER_PLATFORM "x86 (32-bit)"
#elif defined(__aarch64__)
#  define COMPILER_PLATFORM "ARM64"
#elif defined(__arm__)
#  define COMPILER_PLATFORM "ARM (32-bit)"
#elif defined(__powerpc64__) || defined(__ppc64__)
#  define COMPILER_PLATFORM "PowerPC 64"
#elif defined(__powerpc__) || defined(__ppc__)
#  define COMPILER_PLATFORM "PowerPC 32"
#elif defined(__riscv) && (__riscv_xlen == 64)
#  define COMPILER_PLATFORM "RISC-V 64"
#elif defined(__riscv) && (__riscv_xlen == 32)
#  define COMPILER_PLATFORM "RISC-V 32"
#else
#  define COMPILER_PLATFORM "unknown"
#endif

static inline void print_configuration(void) {
	printf("===== SCHEME CONFIG =====\r\n");
	printf("[API] Algo Name: " CRYPTO_ALGNAME "\r\n");
	printf("[API] Algo Version: " CRYPTO_VERSION "\r\n");
	printf("Instruction Sets:");
	printf(" [Platform = %s]", COMPILER_PLATFORM);
#ifdef __SSE__
	printf(" SSE");
#endif
#ifdef __AVX__
	printf(" AVX");
#endif
#ifdef __AVX2__
	printf(" AVX2");
#endif
#if defined(__AVX512BW__) && defined(__AVX512F__) && defined(__AVX512VL__) && defined(__AVX512VPOPCNTDQ__) && defined(__AVX512VBMI__)
	printf(" AVX512BW AVX512F AVX512VL AVX512VPOPCNTDQ AVX512VBMI");
#endif
#if defined(__GFNI__) && !defined(NO_GFNI)
	printf(" GFNI");
#endif
#ifdef __AES__
	printf(" AES-NI");
#endif
	printf("\r\n");

	printf("Configuration elements:\r\n");
#if defined(KECCAK_PLATFORM)
	printf("  Keccak implementation: %s\r\n", STR(KECCAK_PLATFORM));
#else
	printf("  Keccak implementation: unkown\r\n");
#endif

#ifdef MEMORY_EFFICIENT_KEYGEN
	printf("  Keygen: memopt\r\n");
#else
	printf("  Keygen: default\r\n");
#endif
#ifdef VERIFY_MEMOPT
	printf("  Verify: memopt\r\n");
#else
	printf("  Verify: default\r\n");
#endif
#if defined(PIOP_BITSLICE)
	printf("  PIOP: bitslice (see below for specific implementation type)\r\n");
#elif defined(MEMORY_EFFICIENT_PIOP)
	printf("  PIOP: memopt\r\n");
#ifdef PIOP_NB_PARALLEL_REPETITIONS_SIGN
	printf("    PIOP_NB_PARALLEL_REPETITIONS_SIGN %d\r\n", PIOP_NB_PARALLEL_REPETITIONS_SIGN);
#else
	printf("    PIOP_NB_PARALLEL_REPETITIONS_SIGN %d (default)\r\n", MQOM2_PARAM_TAU);
#endif
#ifdef PIOP_NB_PARALLEL_REPETITIONS_VERIFY
	printf("    PIOP_NB_PARALLEL_REPETITIONS_VERIFY %d\r\n", PIOP_NB_PARALLEL_REPETITIONS_VERIFY);
#else
	printf("    PIOP_NB_PARALLEL_REPETITIONS_VERIFY %d (default)\r\n", MQOM2_PARAM_TAU);
#endif

#else
	printf("  PIOP: default\r\n");
#endif

#if defined(USE_PIOP_CACHE) && !defined(MEMORY_EFFICIENT_PIOP)
	printf("    PIOP cache ON\r\n");
#else
	printf("    PIOP cache OFF\r\n");
#endif

	BLC_PrintConfig();

	printf("  Rijndael implementation: %s\r\n", rijndael_conf);
	printf("  Rijndael public implementation: %s\r\n", rijndael_conf_pub);
	printf("  Fields implementation: %s\r\n", fields_conf);
#if defined(PIOP_BITSLICE)
	printf("  Fields bitslice implementation: %s\r\n", fields_bitslice_conf);
#endif

	printf("  MISC options:\r\n");
#ifdef USE_ENC_X8
	printf("    - USE_ENC_X8: ON\r\n");
#else
	printf("    - USE_ENC_X8: OFF\r\n");
#endif
#ifdef USE_XOF_X4
	printf("    - USE_XOF_X4: ON\r\n");
#else
	printf("    - USE_XOF_X4: OFF\r\n");
#endif

#ifdef PRG_ONE_RIJNDAEL_CTX
	printf("    - PRG_ONE_RIJNDAEL_CTX: ON (forcing only one Rijndael ctx for PRG_x1)\r\n");
#endif

#ifdef USE_PRG_CACHE
#ifdef NO_EXPANDMQ_PRG_CACHE
	printf("    - NO_EXPANDMQ_PRG_CACHE: ON (i.e. NO cache used for expand MQ)\r\n");
#else
	printf("    - NO_EXPANDMQ_PRG_CACHE: OFF (i.e. cache is used for expand MQ)\r\n");
#endif
#endif

#ifdef USE_PRG_CACHE
#ifdef NO_BLC_PRG_CACHE
	printf("    - NO_BLC_PRG_CACHE: ON (i.e. NO cache used for BLC)\r\n");
#else
	printf("    - NO_BLC_PRG_CACHE: OFF (i.e. cache is used for BLC)\r\n");
#endif
#endif

#ifdef USE_ENC_CTX_CLEANSING
	printf("    - USE_ENC_CTX_CLEANSING: ON\r\n");
#else
	printf("    - USE_ENC_CTX_CLEANSING: OFF\r\n");
#endif

#ifdef USE_SIGNATURE_BUFFER_AS_TEMP
	printf("    - USE_SIGNATURE_BUFFER_AS_TEMP active (i.e. use signature buffer as temporary storage)\r\n");
#endif

#ifndef NDEBUG
	printf("Debug: On\r\n");
#else
	printf("Debug: Off\r\n");
#endif
}

#endif /* MQOM_BENCH_UTILS_H */
