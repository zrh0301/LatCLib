/* SPDX-License-Identifier: MIT */
#ifndef SNOVA_R3_PROFILE_H
#define SNOVA_R3_PROFILE_H

#include <stdint.h>

#ifdef SNOVA_PROFILE_R3

enum snova_r3_profile_stage {
	SNOVA_R3_PROFILE_HASH,
	SNOVA_R3_PROFILE_VINEGAR,
	SNOVA_R3_PROFILE_WHIP_SIG,
	SNOVA_R3_PROFILE_SUM_T0_MEMSET,
	SNOVA_R3_PROFILE_SUM_T0_AXPY,
	SNOVA_R3_PROFILE_SUM_T0_COMPRESS,
	SNOVA_R3_PROFILE_SUM_T1,
	SNOVA_R3_PROFILE_CENTRAL_E,
	SNOVA_R3_PROFILE_F21,
	SNOVA_R3_PROFILE_F12,
	SNOVA_R3_PROFILE_F12_TRANSPOSE,
	SNOVA_R3_PROFILE_XO_F21_COEFF,
	SNOVA_R3_PROFILE_XO_F21_AXPY,
	SNOVA_R3_PROFILE_XO_F21_REORDER,
	SNOVA_R3_PROFILE_XO_F12_COEFF,
	SNOVA_R3_PROFILE_XO_F12_AXPY,
	SNOVA_R3_PROFILE_XO_F12_REORDER,
	SNOVA_R3_PROFILE_GAUSS_INPUT,
	SNOVA_R3_PROFILE_GAUSS_FORWARD,
	SNOVA_R3_PROFILE_BACKSUB,
	SNOVA_R3_PROFILE_T12,
	SNOVA_R3_PROFILE_PACK,
	SNOVA_R3_PROFILE_STAGE_COUNT
};

struct snova_r3_profile {
	uint64_t calls;
	uint64_t attempts;
	uint64_t cycles[SNOVA_R3_PROFILE_STAGE_COUNT];
};

extern struct snova_r3_profile snova_r3_profile_data;

static inline uint64_t snova_r3_profile_read_cycle(void) {
#if defined(__riscv)
	uint64_t value;
	__asm__ __volatile__("rdcycle %0" : "=r"(value));
	return value;
#else
	return 0;
#endif
}

void snova_r3_profile_reset(void);
const struct snova_r3_profile *snova_r3_profile_get(void);

#define SNOVA_R3_PROFILE_BEGIN(name) \
	uint64_t snova_r3_profile_start_##name = snova_r3_profile_read_cycle()
#define SNOVA_R3_PROFILE_END(name, stage) \
	snova_r3_profile_data.cycles[stage] += \
		snova_r3_profile_read_cycle() - snova_r3_profile_start_##name
#define SNOVA_R3_PROFILE_CALL() (++snova_r3_profile_data.calls)
#define SNOVA_R3_PROFILE_ATTEMPT() (++snova_r3_profile_data.attempts)

#else

#define SNOVA_R3_PROFILE_BEGIN(name) ((void)0)
#define SNOVA_R3_PROFILE_END(name, stage) ((void)0)
#define SNOVA_R3_PROFILE_CALL() ((void)0)
#define SNOVA_R3_PROFILE_ATTEMPT() ((void)0)

#endif

#endif
