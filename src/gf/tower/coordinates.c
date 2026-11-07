/* Coordinate tables from GF_cal/gf256_tower.c; provenance in sources.json. */
#include "gf256_tower.h"
#include "gf16.h"
static const uint8_t TOWER_PHI_LO[16] = {
    0x00, 0x01, 0x20, 0x21, 0x46, 0x47, 0x66, 0x67,
    0x4c, 0x4d, 0x6c, 0x6d, 0x0a, 0x0b, 0x2a, 0x2b
};
static const uint8_t TOWER_PHI_HI[16] = {
    0x00, 0x3c, 0xd5, 0xe9, 0x34, 0x08, 0xe1, 0xdd,
    0xe5, 0xd9, 0x30, 0x0c, 0xd1, 0xed, 0x04, 0x38
};
static const uint8_t TOWER_PSI_LO[16] = {
    0x00, 0x01, 0x5c, 0x5d, 0xe0, 0xe1, 0xbc, 0xbd,
    0x50, 0x51, 0x0c, 0x0d, 0xb0, 0xb1, 0xec, 0xed
};
static const uint8_t TOWER_PSI_HI[16] = {
    0x00, 0xa2, 0x02, 0xa0, 0xb8, 0x1a, 0xba, 0x18,
    0xdb, 0x79, 0xd9, 0x7b, 0x63, 0xc1, 0x61, 0xc3
};
/* Fixed-index scans avoid secret-indexed scalar memory lookups. */
static uint8_t select16(const uint8_t *t,uint8_t x){uint8_t v=0;for(unsigned i=0;i<16;i++)v^=t[i]&(uint8_t)(0u-(unsigned)(i==x));return v;}
uint8_t gf256_tower_encode(uint8_t x){return select16(TOWER_PHI_LO,x&15)^select16(TOWER_PHI_HI,x>>4);}
uint8_t gf256_tower_decode(uint8_t x){return select16(TOWER_PSI_LO,x&15)^select16(TOWER_PSI_HI,x>>4);}
uint8_t gf256_tower_product(uint8_t a,uint8_t b){uint8_t m1=gf16_mul(a&15,b&15),m2=gf16_mul(a>>4,b>>4),m3=gf16_mul((a&15)^(a>>4),(b&15)^(b>>4));return (m1^gf16_mul(8,m2))|((m1^m3)<<4);}
