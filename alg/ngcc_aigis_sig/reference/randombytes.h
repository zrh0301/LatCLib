#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

/* Aigis-Sig+ draws its key-generation seed from the ICCS DRNG so that the
 * published KAT schedule is reproducible. The reference declared this as a
 * void function; it now reports DRNG failure (sign-01-10). */

#include <stddef.h>

int randombytes(unsigned char *x, size_t xlen);

#endif
