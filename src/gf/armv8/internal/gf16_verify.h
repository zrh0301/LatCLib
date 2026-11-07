#ifndef GF16_VERIFY_H
#define GF16_VERIFY_H

#include "gf16.h"

#include <stddef.h>
#include <stdint.h>

typedef enum {
    GF16_VERIFY_TABLE = 0,
    GF16_VERIFY_OFFICIAL_EXPANDED = 1,
    GF16_VERIFY_SHUFFLE_BITPLANE = 2,
    GF16_VERIFY_PACKED_SPLIT_FULL = 3,
    GF16_VERIFY_PACKED_SPLIT_COMPUTE = 4
} pqc_kernel_gf16_verify_method;

typedef struct pqc_kernel_gf16_verify_case pqc_kernel_gf16_verify_case;

typedef struct {
    uint64_t ps_ns;
    uint64_t variable_ns;
    uint64_t secondary_ns;
    uint64_t pack_ns;
    uint64_t unpack_ns;
    uint64_t total_ns;
} pqc_kernel_gf16_verify_metrics;

/*
 * 创建一个 l=2 或 l=4 的完整 verify workload。
 * rows/cols/lanes 是抽象收缩尺寸；l 表示并行的独立 GF(16) 流数量。
 * 该 workload 不保存密钥或跨调用缓存，只保存本次测试的输入和工作区。
 */
pqc_kernel_gf16_verify_case *pqc_kernel_gf16_verify_create(size_t l, size_t rows, size_t cols,
                                     size_t lanes, uint32_t seed);
void pqc_kernel_gf16_verify_destroy(pqc_kernel_gf16_verify_case *c);

/* 执行：P*S -> S^T*T -> secondary，并生成最终输出。 */
int pqc_kernel_gf16_verify_run(pqc_kernel_gf16_verify_case *c, pqc_kernel_gf16_verify_method method,
                    pqc_kernel_gf16_verify_metrics *metrics);

/* 对所有方法的输出做同一参考结果对拍。 */
int pqc_kernel_gf16_verify_check(pqc_kernel_gf16_verify_case *c, pqc_kernel_gf16_verify_method method);

const char *pqc_kernel_gf16_verify_method_name(pqc_kernel_gf16_verify_method method);
const uint8_t *pqc_kernel_gf16_verify_output(const pqc_kernel_gf16_verify_case *c);
size_t pqc_kernel_gf16_verify_output_bytes(const pqc_kernel_gf16_verify_case *c);
size_t pqc_kernel_gf16_verify_l(const pqc_kernel_gf16_verify_case *c);

#endif
