#ifndef GF16_PLATFORM_H
#define GF16_PLATFORM_H

#include <arm_neon.h>
#include <stdint.h>

#define GF16_BACKEND_NAME "armv8-neon"
#define GF16_VEC_BYTES 16

typedef uint8x16_t pqc_kernel_gf16_vec;

static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_load(const void *p) { return vld1q_u8((const uint8_t *)p); }
static inline void pqc_kernel_gf16_vec_store(void *p, pqc_kernel_gf16_vec v) { vst1q_u8((uint8_t *)p, v); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_zero(void) { return vdupq_n_u8(0); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_set1(uint8_t x) { return vdupq_n_u8(x); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return veorq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_and(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return vandq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return vorrq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_eq(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return vceqq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_gt_u8_small(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return vcgtq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_shuffle(pqc_kernel_gf16_vec tab, pqc_kernel_gf16_vec idx) { return vqtbl1q_u8(tab, idx); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_srli4(pqc_kernel_gf16_vec a) { return vandq_u8(vshrq_n_u8(a, 4), vdupq_n_u8(0x0f)); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_slli4(pqc_kernel_gf16_vec a) { return vshlq_n_u8(a, 4); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_add_u8(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return vaddq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_sub_u8(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) { return vsubq_u8(a, b); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_mask_u8(pqc_kernel_gf16_vec v, uint8_t m) { return vandq_u8(v, vdupq_n_u8(m)); }
static inline pqc_kernel_gf16_vec pqc_kernel_gf16_vec_not(pqc_kernel_gf16_vec v) { return vmvnq_u8(v); }

#endif
