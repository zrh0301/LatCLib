#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include <memory.h>

#define WBITS 32
#define WBYTS (WBITS>>3)
#define LQBITS 32
#define LQBYTS (LQBITS>>3)

//  WBIT = 32
typedef uint32_t dig_t;
typedef uint64_t udi_t;

#define WMASK 0xfffffff

//LQBIT = 32
typedef int32_t lssi_t;//lattice signed single 
typedef uint32_t lusi_t;
typedef int64_t lsdi_t;
typedef uint64_t ludi_t;

#endif