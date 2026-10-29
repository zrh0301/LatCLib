/* SPDX-License-Identifier: MIT */
#ifndef SNOVA_R3_K230_RVV_COMPAT_H
#define SNOVA_R3_K230_RVV_COMPAT_H

/* The K230 SDK GCC 12 toolchain uses the older unprefixed RVV spellings. */
#if defined(__riscv_vector) && defined(__GNUC__) && !defined(__clang__) && (__GNUC__ < 13)
#define __riscv_vsetvl_e16m2 vsetvl_e16m2
#define __riscv_vle16_v_u16m2 vle16_v_u16m2
#define __riscv_vse16_v_u16m2 vse16_v_u16m2
#define __riscv_vxor_vv_u16m2 vxor_vv_u16m2
#define __riscv_vmul_vx_u16m2 vmul_vx_u16m2
#define __riscv_vand_vx_u16m2 vand_vx_u16m2
#define __riscv_vsrl_vx_u16m2 vsrl_vx_u16m2
#endif

#endif
