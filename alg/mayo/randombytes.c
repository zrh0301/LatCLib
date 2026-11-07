#include "rand.h"
#include <stdlib.h>
#include <stdint.h>
int pqc_mayo_randombytes(unsigned char *out,size_t n){
 if(n>SIZE_MAX) return -1;
 int rc=pqc_randombytes(out,(size_t)n,NULL,NULL);
 /* Some imported call sites ignore the status. Fail closed on RNG failure. */
 if(rc) abort();
 return rc;
}
