// SPDX-License-Identifier: MIT

/**
 * Tool to create ABQ source files.
 *
 * Copyright (c) 2026 SNOVA TEAM
 */

#define SNOVA_OPT 0

#include "snova_opt_16.c"
#include "symmetric.c"
#include "aes.c"
#include "api.h"

int main(void) {
	// SPDX-License-Identifier: MIT

	printf("// SPDX-License-Identifier: MIT\n\n");
	printf("/**\n * ABQ data %s  (o=%d, l=%d, r=%d)\n", CRYPTO_ALGNAME, SNOVA_o, SNOVA_l, SNOVA_r);
	printf(" *\n * Copyright (c) 2026 SNOVA TEAM\n */\n\n");

	snova_init();

#if 1
	printf("static const uint16_t fixedAm[SNOVA_o * SNOVA_alpha * SNOVA_r2] = {\n");
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_r2; ++idx) {
		if (idx % 16 == 15) {
			printf(",\n\t");
		} else if (idx) {
			printf(", ");
		} else {
			printf("\t");
		}
		printf("0x%x", fixedAm[idx]);
	}
	printf("\n};\n\n");

	printf("static const uint16_t fixedBm[SNOVA_o * SNOVA_alpha * SNOVA_lr] = {\n");
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_lr; ++idx) {
		if (idx % 16 == 15) {
			printf(",\n\t");
		} else if (idx) {
			printf(", ");
		} else {
			printf("\t");
		}
		printf("0x%x", fixedBm[idx]);
	}
	printf("\n};\n\n");

	printf("static const uint16_t fixedq1[SNOVA_o * SNOVA_alpha * SNOVA_l] = {\n");
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_l; ++idx) {
		if (idx % 16 == 15) {
			printf(",\n\t");
		} else if (idx) {
			printf(", ");
		} else {
			printf("\t");
		}
		printf("0x%x", fixedq1[idx]);
	}
	printf("\n};\n\n");

	printf("static const uint16_t fixedq2[SNOVA_o * SNOVA_alpha * SNOVA_l] = {\n");
	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha * SNOVA_l; ++idx) {
		if (idx % 16 == 15) {
			printf(",\n\t");
		} else if (idx) {
			printf(", ");
		} else {
			printf("\t");
		}
		printf("0x%x", fixedq2[idx]);
	}
	printf("\n};\n");
#else
	printf("static uint16_t gf_S[SNOVA_l * SNOVA_l2] = {\n");
	for (int idx = 0; idx < SNOVA_l * SNOVA_l2; ++idx) {
		if (idx % 16 == 15) {
			printf(",\n\t");
		} else if (idx) {
			printf(", ");
		} else {
			printf("\t");
		}
		printf("0x%x", gf16_compress(gf_Sx[idx]));
	}
	printf("\n};\n\n");

	printf("static uint16_t gf_Sx[SNOVA_l * SNOVA_l2] = {\n");
	for (int idx = 0; idx < SNOVA_l * SNOVA_l2; ++idx) {
		if (idx % 16 == 15) {
			printf(",\n\t");
		} else if (idx) {
			printf(", ");
		} else {
			printf("\t");
		}
		printf("0x%x", gf_Sx[idx]);
	}
	printf("\n};\n\n");
#endif
}
