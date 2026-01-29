#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// 引入头文件 (假设编译时指定了 -I../inc)
#include "../inc/types.h"
#include "../inc/poly.h"

// =========================================================
// 0. CPU Cycle 计数器 (x86/x64)
// =========================================================
static inline uint64_t cpucycles(void) {
    uint64_t result;
    __asm__ volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
    : "=a" (result) : : "%rdx"); 
    return result;
}

// =========================================================
// 1. 核心技巧：包含 poly_ntt.c 两次以生成两套函数
// =========================================================

/* 
 * 实例化 LAYER 8 (Full NTT)
 * 我们通过宏将函数名和静态变量 zetas 重命名，防止符号冲突
 */
#define LAYER 8
#define poly_ntt        poly_ntt_full
#define poly_intt       poly_intt_full
#define poly_pointwise  poly_pointwise_full
#define zetas           zetas_full  
// 直接包含源文件
#include "poly_ntt.c" 
// 清理宏
#undef LAYER
#undef poly_ntt
#undef poly_intt
#undef poly_pointwise
#undef zetas

/* 
 * 实例化 LAYER 7 (Partial NTT)
 */
#define LAYER 7
#define poly_ntt        poly_ntt_part
#define poly_intt       poly_intt_part
#define poly_pointwise  poly_pointwise_part
#define zetas           zetas_part
#include "poly_ntt.c"
#undef LAYER
#undef poly_ntt
#undef poly_intt
#undef poly_pointwise
#undef zetas


// =========================================================
// 2. 主测试程序
// =========================================================
int main() {
    srand((unsigned int)time(NULL));

    // 初始化 Dilithium 参数 (Context)
    lp_ctx ctx;
    ctx.n = 256;
    ctx.q = 8380417;
    ctx.qc = 58728449; // -q^{-1} mod 2^32

    // 分配内存
    lssi_t a[MAX_POLY_COEFS];
    lssi_t b[MAX_POLY_COEFS];
    lssi_t t_full[MAX_POLY_COEFS]; // Full 模式的工作区
    lssi_t t_part[MAX_POLY_COEFS]; // Partial 模式的工作区

    printf("==================================================\n");
    printf("     Dilithium NTT Validator (Full vs Partial)    \n");
    printf("==================================================\n");

    // -----------------------------------------------------
    // 测试 1: 正确性对比 (Cross Validation)
    // -----------------------------------------------------
    printf("\n[1] Correctness Check (Random Inputs)...\n");
    
    int rounds = 1000;
    int passed_count = 0;

    for (int i = 0; i < rounds; i++) {
        // 1. 使用 src/poly.c 中已有的函数生成数据
        poly_rand(a, ctx.n); 
        poly_rand(b, ctx.n);

        // 2. 路径 A: Full NTT (8层)
        // 计算: INTT( NTT(a) * NTT(b) )
        poly_copy(t_full, a, ctx.n); // 复制 a 到工作区
        lssi_t b_tmp[MAX_POLY_COEFS];
        poly_copy(b_tmp, b, ctx.n);  // 复制 b 到临时区

        poly_ntt_full(t_full, &ctx);
        poly_ntt_full(b_tmp, &ctx);
        poly_pointwise_full(t_full, t_full, b_tmp, &ctx); // In-place mul
        poly_intt_full(t_full, &ctx);

        // 3. 路径 B: Partial NTT (7层)
        poly_copy(t_part, a, ctx.n);
        poly_copy(b_tmp, b, ctx.n);

        poly_ntt_part(t_part, &ctx);
        poly_ntt_part(b_tmp, &ctx);
        poly_pointwise_part(t_part, t_part, b_tmp, &ctx); // In-place mul (Basecase)
        poly_intt_part(t_part, &ctx);

        // 4. 比较结果
        // 使用 src/poly.c 中已有的比较函数
        if (!poly_eq(t_full, t_part, ctx.n)) {
            printf("[FAIL] Mismatch at round %d!\n", i);
            // 打印前几个系数以便调试
            printf("Expected: %d %d %d ...\n", t_full[0], t_full[1], t_full[2]);
            printf("Got     : %d %d %d ...\n", t_part[0], t_part[1], t_part[2]);
            return 1;
        } else {
            passed_count++;
        }
    }
    printf("[PASS] %d rounds matched perfectly.\n", passed_count);


    // -----------------------------------------------------
    // 测试 2: 性能测试 (CPU Cycles)
    // -----------------------------------------------------
    printf("\n[2] Performance Benchmark (Average Cycles)...\n");
    
    int bench_loops = 50000;
    uint64_t start, end;
    uint64_t sum_full = 0, sum_part = 0;

    // 准备数据 (预先转到 NTT 域，模拟实际乘法场景)
    poly_rand(a, ctx.n);
    poly_rand(b, ctx.n);
    
    // 预热 Full
    poly_copy(t_full, a, ctx.n);
    poly_ntt_full(t_full, &ctx); 
    poly_copy(t_part, b, ctx.n);
    poly_ntt_full(t_part, &ctx);

    // --- Bench Full (NTT + Mul + INTT) ---
    // 为了测得准，我们包含 NTT, Mul, INTT 整个链路
    // 注意：每次循环我们不重置数据，只让数据不断变化，避免 memcpy 干扰测速
    for(int i=0; i<bench_loops; i++) {
        start = cpucycles();
        poly_ntt_full(t_full, &ctx);
        poly_pointwise_full(t_full, t_full, t_part, &ctx);
        poly_intt_full(t_full, &ctx);
        end = cpucycles();
        sum_full += (end - start);
    }

    // 预热 Part
    poly_copy(t_full, a, ctx.n);
    poly_ntt_part(t_full, &ctx); 
    poly_copy(t_part, b, ctx.n);
    poly_ntt_part(t_part, &ctx);

    // --- Bench Part (NTT + Mul + INTT) ---
    for(int i=0; i<bench_loops; i++) {
        start = cpucycles();
        poly_ntt_part(t_full, &ctx);
        poly_pointwise_part(t_full, t_full, t_part, &ctx);
        poly_intt_part(t_full, &ctx);
        end = cpucycles();
        sum_part += (end - start);
    }

    uint64_t avg_full = sum_full / bench_loops;
    uint64_t avg_part = sum_part / bench_loops;

    printf("Op: NTT + Mul + INTT\n");
    printf("Full Layer-8 : %lu cycles\n", avg_full);
    printf("Part Layer-7 : %lu cycles\n", avg_part);
    
    int64_t diff = (int64_t)avg_part - (int64_t)avg_full;
    printf("Diff (7L - 8L): %ld cycles (%.2f%%)\n", diff, 100.0 * diff / avg_full);

    return 0;
}