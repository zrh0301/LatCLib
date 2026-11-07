#include "gfext.h"

#include <string.h>

/*
 * GF(256) 矩阵级场景（正文 5.6.5 节）。
 *
 *   M1 常数系数收缩:  out[r][lane] = XOR_c coeff[r][c] * sig[c][lane]
 *   M2 变量矩阵收缩:  out[r][lane] = XOR_c left[(r*cols+c)][lane] * right[c][lane]
 *
 * 计时口径与 3.0 节一致：一次“逻辑 GF(256) 乘法”指一个 (r,c,lane) 三元组，
 * 归一化时除以 rows*cols*lanes；独立收缩测试用 rows=cols=lanes=64。
 *
 * 本文件只放通用乘法族的矩阵实现（标量、SWAR、SIMD shift/reduce、GFNI）；
 * 塔分解的矩阵内核在 pqc_kernel_gf256_tower.c，避免两种口径混在一个文件里。
 * GF256_SIMD_MUL 由平台头文件给出（x86 是 AVX2、ARMv8 是 NEON）。
 */

/* --------------------------------------------------------------- M1 常数收缩 */

void pqc_kernel_gf256_ps_scalar(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                     size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        for (size_t lane = 0; lane < lanes; ++lane) out[r * lanes + lane] = 0;
    }
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t k = coeff[r * cols + c];
            const uint8_t *s = sig + c * lanes;
            for (size_t lane = 0; lane < lanes; ++lane) {
                o[lane] ^= pqc_kernel_gf256_mul_ref(s[lane], k);
            }
        }
    }
}

void pqc_kernel_gf256_ps_swar(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                   size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        memset(o, 0, lanes);
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t k = coeff[r * cols + c];
            const uint32_t kv = (uint32_t)k * 0x01010101u;
            const uint8_t *s = sig + c * lanes;
            size_t lane = 0;
            for (; lane + 4u <= lanes; lane += 4u) {
                uint32_t sv, ov;
                memcpy(&sv, s + lane, sizeof(sv));
                memcpy(&ov, o + lane, sizeof(ov));
                ov ^= pqc_kernel_gf256_mul4_swar_word(sv, kv);
                memcpy(o + lane, &ov, sizeof(ov));
            }
            for (; lane < lanes; ++lane) o[lane] ^= pqc_kernel_gf256_mul_ref(s[lane], k);
        }
    }
}

void pqc_kernel_gf256_ps_simd(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                   size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES)
            pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_zero());
        for (; i < lanes; ++i) o[i] = 0;
        for (size_t c = 0; c < cols; ++c) {
            const pqc_kernel_gf16_vec kv = pqc_kernel_gf16_vec_set1(coeff[r * cols + c]);
            const uint8_t *s = sig + c * lanes;
            i = 0;
            for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES) {
                const pqc_kernel_gf16_vec p = GF256_SIMD_MUL(pqc_kernel_gf16_vec_load(s + i), kv);
                pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(o + i), p));
            }
            for (; i < lanes; ++i) o[i] ^= pqc_kernel_gf256_mul_ref(s[i], coeff[r * cols + c]);
        }
    }
}

#if defined(GF256_MATRIX_HAS_GFNI)
void pqc_kernel_gf256_ps_gfni(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                   size_t rows, size_t cols, size_t lanes) {
    if (!pqc_kernel_gf256_gfni_available()) return;
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES)
            pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_zero());
        for (; i < lanes; ++i) o[i] = 0;
        for (size_t c = 0; c < cols; ++c) {
            const pqc_kernel_gf16_vec kv = pqc_kernel_gf16_vec_set1(coeff[r * cols + c]);
            const uint8_t *s = sig + c * lanes;
            i = 0;
            for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES) {
                const pqc_kernel_gf16_vec p = pqc_kernel_gf256_vec_mul_gfni(pqc_kernel_gf16_vec_load(s + i), kv);
                pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(o + i), p));
            }
            for (; i < lanes; ++i) o[i] ^= pqc_kernel_gf256_mul_ref(s[i], coeff[r * cols + c]);
        }
    }
}
#endif /* GF256_MATRIX_HAS_GFNI */

/* --------------------------------------------------------------- M2 变量收缩 */

void pqc_kernel_gf256_sts_scalar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                      size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        for (size_t lane = 0; lane < lanes; ++lane) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf256_mul_ref(left[(r * cols + c) * lanes + lane],
                                     right[c * lanes + lane]);
            }
            out[r * lanes + lane] = acc;
        }
    }
}

void pqc_kernel_gf256_sts_swar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                    size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        memset(o, 0, lanes);
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *s = left + (r * cols + c) * lanes;
            const uint8_t *t = right + c * lanes;
            size_t lane = 0;
            for (; lane + 4u <= lanes; lane += 4u) {
                uint32_t sv, tv, ov;
                memcpy(&sv, s + lane, sizeof(sv));
                memcpy(&tv, t + lane, sizeof(tv));
                memcpy(&ov, o + lane, sizeof(ov));
                ov ^= pqc_kernel_gf256_mul4_swar_word(sv, tv);
                memcpy(o + lane, &ov, sizeof(ov));
            }
            for (; lane < lanes; ++lane) o[lane] ^= pqc_kernel_gf256_mul_ref(s[lane], t[lane]);
        }
    }
}

void pqc_kernel_gf256_sts_simd(uint8_t *out, const uint8_t *left, const uint8_t *right,
                    size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES)
            pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_zero());
        for (; i < lanes; ++i) o[i] = 0;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *s = left + (r * cols + c) * lanes;
            const uint8_t *t = right + c * lanes;
            i = 0;
            for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES) {
                const pqc_kernel_gf16_vec p = GF256_SIMD_MUL(pqc_kernel_gf16_vec_load(s + i),
                                                      pqc_kernel_gf16_vec_load(t + i));
                pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(o + i), p));
            }
            for (; i < lanes; ++i) o[i] ^= pqc_kernel_gf256_mul_ref(s[i], t[i]);
        }
    }
}

#if defined(GF256_MATRIX_HAS_GFNI)
void pqc_kernel_gf256_sts_gfni(uint8_t *out, const uint8_t *left, const uint8_t *right,
                    size_t rows, size_t cols, size_t lanes) {
    if (!pqc_kernel_gf256_gfni_available()) return;
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        size_t i = 0;
        for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES)
            pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_zero());
        for (; i < lanes; ++i) o[i] = 0;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *s = left + (r * cols + c) * lanes;
            const uint8_t *t = right + c * lanes;
            i = 0;
            for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES) {
                const pqc_kernel_gf16_vec p = pqc_kernel_gf256_vec_mul_gfni(pqc_kernel_gf16_vec_load(s + i),
                                                      pqc_kernel_gf16_vec_load(t + i));
                pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(o + i), p));
            }
            for (; i < lanes; ++i) o[i] ^= pqc_kernel_gf256_mul_ref(s[i], t[i]);
        }
    }
}
#endif /* GF256_MATRIX_HAS_GFNI */
