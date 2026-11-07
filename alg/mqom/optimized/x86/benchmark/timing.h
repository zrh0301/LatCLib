#ifndef MQOM_TIMING_H
#define MQOM_TIMING_H

#include <stdint.h>
#include <stddef.h>
/* On POSIX platforms, use the timeval structure */
#if defined(__unix__) || (defined (__APPLE__) && defined (__MACH__)) || defined(_WIN32)
#include <time.h>
#include <sys/time.h>
#else
/* Other non-POSIX platforms */
/* Define our timeval structure */
struct timeval {
	long int tv_sec;
	long int tv_usec;
};
/* The gettimeofday API is external and should be provided by the user */
extern void gettimeofday(struct timeval*, void*);
#endif

/* Namespacing with the appropriate prefix */
#ifndef MQOM_NAMESPACE
#ifdef APPLY_NAMESPACE
#ifndef concat2
#define _concat2(a, b) a ## b
#define concat2(a, b) _concat2(a, b)
#endif
#define MQOM_NAMESPACE(s) concat2(APPLY_NAMESPACE, s)
#else
#define MQOM_NAMESPACE(s) s
#endif
#endif

/* Deal with namespacing */
#define btimer_init MQOM_NAMESPACE(btimer_init)
#define btimer_start MQOM_NAMESPACE(btimer_start)
#define btimer_count MQOM_NAMESPACE(btimer_count)
#define btimer_end MQOM_NAMESPACE(btimer_end)
#define btimer_diff MQOM_NAMESPACE(btimer_diff)
#define btimer_diff_cycles MQOM_NAMESPACE(btimer_diff_cycles)
#define btimer_get MQOM_NAMESPACE(btimer_get)
#define btimer_get_cycles MQOM_NAMESPACE(btimer_get_cycles)

typedef struct btimer_t {
	unsigned int counter;
	// gettimeofday
	double nb_milliseconds;
	struct timeval start, stop;
	// rdtscp or RDPMC
	uint64_t nb_cycles;
	unsigned int garbage;
	uint64_t cstart, cstop;
} btimer_t;

void btimer_init(btimer_t* timer);
void btimer_start(btimer_t *timer);
void btimer_count(btimer_t *timer);
void btimer_end(btimer_t *timer);
double btimer_diff(btimer_t *timer);
uint64_t btimer_diff_cycles(btimer_t *timer);
double btimer_get(btimer_t *timer);
double btimer_get_cycles(btimer_t *timer);

#endif /* MQOM_TIMING_H */
