#if defined(BENCHMARK_CYCLES) && defined(__linux__)
#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
static inline void set_cpu_affinity(int cpu) {
	cpu_set_t set;

	CPU_ZERO(&set);
	CPU_SET(cpu, &set);
	if (sched_setaffinity(getpid(), sizeof(set), &set) == -1) {
		fprintf(stderr, "Error: error when setting affinity to CPU 0 ...\n");
		exit(-1);
	}
}
#endif

#include "timing.h"

#ifdef BENCHMARK_CYCLES
/* ====================================================== */
/* Getting cycles primitives depending on the platform */
#if defined(__linux__) && (defined(__amd64__) || defined(__x86_64__))
/* On Linux specifically, we use the perf events */
/* NOTE: stolen and adapted from https://cpucycles.cr.yp.to/libcpucycles-20240318/cpucycles/amd64-pmc.c.html */
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>

struct perf_event_attr attr;
int fdperf = -1;
struct perf_event_mmap_page *buf = 0;

long long ticks(void) {
	long long result;
	unsigned int seq;
	long long index;
	long long offset;

	if (buf == 0) {
		return 0;
	}
	do {
		seq = buf->lock;
		asm volatile("" ::: "memory");
		index = buf->index;
		offset = buf->offset;
		asm volatile("rdpmc;shlq $32,%%rdx;orq %%rdx,%%rax"
		             : "=a"(result) : "c"(index-1) : "%rdx");
		asm volatile("" ::: "memory");
	} while (buf->lock != seq);

	result += offset;
	result &= 0xffffffffffff;
	return result;
}

/* NOTE: constructor attribute to be executed first */
__attribute__((constructor)) void ticks_setup(void) {
	/* First of all, set the CPU affinity to CPU 0 to have stable measurements
	 * and avoid issues with P and E Cores on recent Intel CPUs */
	set_cpu_affinity(0);

	if (fdperf == -1) {
		attr.type = PERF_TYPE_HARDWARE;
		attr.config = PERF_COUNT_HW_CPU_CYCLES;
		attr.exclude_kernel = 1;
		attr.exclude_hv = 1;
		fdperf = syscall(__NR_perf_event_open, &attr, 0, -1, -1, 0);
		if (fdperf == -1) {
			fprintf(stderr, "Error: performance counters configuration failed ...\n");
			fprintf(stderr, "  => Please configure RDPMC access with (as superuser) 'echo 2 > /proc/sys/kernel/perf_event_paranoid' (i.e. allow access from userland)\n");
			exit(-1);
		}
		buf = mmap(NULL, sysconf(_SC_PAGESIZE), PROT_READ, MAP_SHARED, fdperf, 0);
	}

	return;
}


long long platform_get_cycles(void) {
	return ticks();
}

#elif defined(__amd64__) || defined(__x86_64__) || defined(__i386__)
/* On other platforms with x86, use rdtscp intrinsics */
#include <x86intrin.h>
long long platform_get_cycles(void) {
	unsigned int garbage;
	return __rdtscp(&garbage);
}
#elif defined(__aarch64__)
long long platform_get_cycles(void) {
	uint64_t result;
	asm volatile ("isb \n mrs %0, CNTVCT_EL0" : "=r" (result));
	return result;
}
#else
/* For other unknown platforms, the platform_get_cycles is externally defined by the user */
extern long long platform_get_cycles(void);
#endif
#endif

void btimer_init(btimer_t* timer) {
	if (timer != NULL) {
		timer->counter = 0;
		timer->nb_milliseconds = 0.;
		timer->nb_cycles = 0;
		timer->start.tv_sec = timer->start.tv_usec = 0;
		timer->stop.tv_sec = timer->stop.tv_usec = 0;
	}
}
void btimer_count(btimer_t *timer) {
	if (timer != NULL) {
		timer->counter++;
	}
}

void btimer_start(btimer_t *timer) {
	if (timer != NULL) {
#ifdef BENCHMARK_TIME
#if defined(CLOCK_MONOTONIC_COARSE) && !defined(BENCHMARK_USE_GETTIMEOFDAY)
		/* NOTE: when available, we use CLOCK_MONOTONIC_COARSE
		 * as it does not require a costly system call */
		struct timespec t;
		clock_gettime(CLOCK_MONOTONIC_COARSE, &t);
		timer->start.tv_sec  = t.tv_sec;
		timer->start.tv_usec = (double)t.tv_nsec * 0.001;
#else
		/* NOTE: on POSIX like systems, this usually requires a syscall, so this can
		 * incur a perfomance hit and perturb the measurements */
		gettimeofday(&timer->start, NULL);
#endif
#else
		(void)timer;
#endif /* BENCHMARK_TIME */
#ifdef BENCHMARK_CYCLES
		timer->cstart = platform_get_cycles();
#endif
	}
}
double btimer_diff(btimer_t *timer) {
	return ( (timer->stop.tv_sec - timer->start.tv_sec) * 1000000 + (timer->stop.tv_usec - timer->start.tv_usec) ) / 1000.;
}
uint64_t btimer_diff_cycles(btimer_t *timer) {
	return (timer->cstop - timer->cstart);
}
void btimer_end(btimer_t *timer) {
	if (timer != NULL) {
#ifdef BENCHMARK_TIME
#if defined(CLOCK_MONOTONIC_COARSE) && !defined(BENCHMARK_USE_GETTIMEOFDAY)
		/* NOTE: when available, we use CLOCK_MONOTONIC_COARSE
		 * as it does not require a costly system call */
		struct timespec t;
		clock_gettime(CLOCK_MONOTONIC_COARSE, &t);
		timer->stop.tv_sec  = t.tv_sec;
		timer->stop.tv_usec = (double)t.tv_nsec * 0.001;
#else
		/* NOTE: on POSIX like systems, this usually requires a syscall, so this can
		 * incur a perfomance hit and perturb the measurements */
		gettimeofday(&timer->stop, NULL);
#endif
		timer->nb_milliseconds += btimer_diff(timer);
#else
		(void)timer;
#endif /* BENCHMARK_TIME */
#ifdef BENCHMARK_CYCLES
		timer->cstop = platform_get_cycles();
		timer->nb_cycles += btimer_diff_cycles(timer);
#endif
	}
}
double btimer_get(btimer_t *timer) {
	return timer->nb_milliseconds / timer->counter;
}
double btimer_get_cycles(btimer_t *timer) {
	return (double)timer->nb_cycles / timer->counter;
}
