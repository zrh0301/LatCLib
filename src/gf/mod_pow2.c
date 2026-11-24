#include "gf/mod_pow2.h"

uint16_t pqc_mod_pow2(uint32_t value, unsigned int logq)
{
    return (uint16_t)(value & ((UINT32_C(1) << logq) - 1U));
}

int32_t pqc_round_pow2(uint16_t value, unsigned int logq, unsigned int bits)
{
    const unsigned int shift = logq - bits;
    const uint32_t v = value & ((UINT32_C(1) << logq) - 1U);
    return (int32_t)((v + (UINT32_C(1) << (shift - 1U))) >> shift);
}
