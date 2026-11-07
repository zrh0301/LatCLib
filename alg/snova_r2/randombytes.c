#include "rand.h"
#include <stdlib.h>
#include <stdint.h>
int pqc_snova_r2_randombytes(unsigned char *out,unsigned long long n){
 if(n>SIZE_MAX) return -1;
 int rc=pqc_randombytes(out,(size_t)n,NULL,NULL);
 /* Some imported call sites ignore the status. Fail closed on RNG failure. */
 if(rc) abort();
 return rc;
}
