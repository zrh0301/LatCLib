#include "gf16.h"

/* SNOVA 风格的 expand。
 * 输入：a 是普通 byte GF16 元素。
 * 输出：一个 uint16_t 扩展态，把 GF16 非零信息摊到 bit-plane 上。
 * 场景：模拟verify / open 中“先展开再延迟压缩”的表示。
 */
uint16_t pqc_kernel_gf16_official_expand(uint8_t a) {
    uint16_t v = (uint16_t)((a & 0x0f) | ((a & 0x0f) << 3) | ((a & 0x0f) << 6) | ((a & 0x0f) << 9));
    return (uint16_t)(v & 0x1111u);
}

/* SNOVA 风格的 compress。
 * 输入：a 是扩展态累加器。
 * 输出：压回普通 byte GF16。
 * 场景：把延迟积累的结果还原成普通域元素。
 */
uint8_t pqc_kernel_gf16_official_compress(uint16_t a) {
    uint16_t v = (uint16_t)((a & 0x000f) ^ ((a & 0x00f0) >> 3) ^ ((a & 0x0f00) >> 6) ^ ((a & 0xf000) >> 9));
    return (uint8_t)((v ^ ((v & 0x00f0) >> 3) ^ (v >> 4)) & 0x000f);
}

/* 扩展态数组版本。
 * 输入：in 为普通 byte GF16 数组。
 * 输出：out 为逐元素扩展态。
 */
void pqc_kernel_gf16_official_expand_array(uint16_t *out, const uint8_t *in, size_t n) {
    size_t i = 0;
#if defined(__AVX2__)
    const __m256i nibble = _mm256_set1_epi16(0x0f);
    const __m256i mask = _mm256_set1_epi16(0x1111);
    for (; i + 8 <= n; i += 8) {
        __m128i r8 = _mm_loadl_epi64((const __m128i *)(in + i));
        __m256i x = _mm256_and_si256(_mm256_cvtepu8_epi16(r8), nibble);
        __m256i v = _mm256_or_si256(x, _mm256_slli_epi16(x, 3));
        v = _mm256_or_si256(v, _mm256_slli_epi16(x, 6));
        v = _mm256_or_si256(v, _mm256_slli_epi16(x, 9));
        _mm_storeu_si128((__m128i *)(out + i),
                         _mm256_castsi256_si128(_mm256_and_si256(v, mask)));
    }
#endif
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_official_expand(in[i]);
    }
}

/* 压缩态数组版本。
 * 输入：in 为扩展态数组。
 * 输出：out 为普通 byte GF16 数组。
 */
void pqc_kernel_gf16_official_compress_array(uint8_t *out, const uint16_t *in, size_t n) {
    size_t i = 0;
#if defined(__AVX2__)
    const __m256i m0 = _mm256_set1_epi16(0x000f);
    const __m256i m1 = _mm256_set1_epi16(0x00f0);
    const __m256i m2 = _mm256_set1_epi16(0x0f00);
    const __m256i m3 = _mm256_set1_epi16((short)0xf000);
    const __m256i m4 = _mm256_set1_epi16(0x00f0);
    for (; i + 8 <= n; i += 8) {
        __m256i a = _mm256_castsi128_si256(
            _mm_loadu_si128((const __m128i *)(in + i)));
        __m256i v = _mm256_xor_si256(_mm256_and_si256(a, m0),
                                     _mm256_srli_epi16(_mm256_and_si256(a, m1), 3));
        v = _mm256_xor_si256(v, _mm256_srli_epi16(_mm256_and_si256(a, m2), 6));
        v = _mm256_xor_si256(v, _mm256_srli_epi16(_mm256_and_si256(a, m3), 9));
        v = _mm256_xor_si256(v, _mm256_srli_epi16(_mm256_and_si256(v, m4), 3));
        v = _mm256_xor_si256(v, _mm256_srli_epi16(v, 4));
        v = _mm256_and_si256(v, m0);
        __m256i packed = _mm256_packus_epi16(v, v);
        _mm_storel_epi64((__m128i *)(out + i), _mm256_castsi256_si128(packed));
    }
#endif
    for (; i < n; ++i) {
        out[i] = pqc_kernel_gf16_official_compress(in[i]);
    }
}

/* raw × expanded 的累加。
 * 输入：raw 是普通 byte GF16；expanded 是展开后的 uint16_t 乘数；acc 是扩展态累加器。
 * 输出：acc[i] ^= raw[i] * expanded[i]。
 * 场景：verify 热路径里的核心乘加。
 */
void pqc_kernel_gf16_official_raw_expanded_muladd(uint16_t *acc, const uint8_t *raw, const uint16_t *expanded, size_t n) {
    size_t i = 0;
#if defined(__AVX2__)
    for (; i + 8 <= n; i += 8) {
        __m128i r8 = _mm_loadl_epi64((const __m128i *)(raw + i));
        __m256i r16 = _mm256_cvtepu8_epi16(r8);
        __m256i e = _mm256_castsi128_si256(
            _mm_loadu_si128((const __m128i *)(expanded + i)));
        __m256i old = _mm256_castsi128_si256(
            _mm_loadu_si128((const __m128i *)(acc + i)));
        __m256i prod = _mm256_mullo_epi16(r16, e);
        _mm_storeu_si128((__m128i *)(acc + i),
                         _mm256_castsi256_si128(_mm256_xor_si256(old, prod)));
    }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    for (; i + 8 <= n; i += 8) {
        uint8x8_t r8 = vld1_u8(raw + i);
        uint16x8_t r16 = vmovl_u8(r8);
        uint16x8_t e = vld1q_u16(expanded + i);
        uint16x8_t old = vld1q_u16(acc + i);
        uint16x8_t prod = vmulq_u16(r16, e);
        vst1q_u16(acc + i, veorq_u16(old, prod));
    }
#endif
    for (; i < n; ++i) {
        acc[i] ^= (uint16_t)((raw[i] & 0x0f) * expanded[i]);
    }
}

void pqc_kernel_gf16_official_const_expanded_muladd(uint16_t *acc, uint8_t k, const uint16_t *expanded, size_t n) {
    size_t i = 0;
#if defined(__AVX2__)
    __m256i kv = _mm256_set1_epi16((short)(k & 0x0f));
    for (; i + 16 <= n; i += 16) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(acc + i));
        __m256i e = _mm256_loadu_si256((const __m256i *)(expanded + i));
        _mm256_storeu_si256((__m256i *)(acc + i),
                            _mm256_xor_si256(a, _mm256_mullo_epi16(kv, e)));
    }
#endif
    for (; i < n; ++i) {
        acc[i] ^= (uint16_t)((k & 0x0f) * expanded[i]);
    }
}
