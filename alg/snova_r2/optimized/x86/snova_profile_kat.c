#include "snova_profile.h"

/*
 * KAT links the same profiled implementation as the benchmark, but it does
 * not measure or print timings. Keep the profile symbols local to this KAT
 * target so benchmark timing remains owned by bench_official.c, including its
 * perf_event_open state.
 */
#if defined(SNOVA_DEEP_PROFILE)

uint64_t snova_profile_cycles_total[SNOVA_PROFILE_PHASE_COUNT];
uint64_t snova_profile_calls[SNOVA_PROFILE_PHASE_COUNT];

uint64_t snova_profile_now(void) {
	return 0;
}

void snova_profile_add(int phase, uint64_t cycles) {
	if (phase < 0 || phase >= SNOVA_PROFILE_PHASE_COUNT) {
		return;
	}
	snova_profile_cycles_total[phase] += cycles;
	snova_profile_calls[phase]++;
}

#endif
