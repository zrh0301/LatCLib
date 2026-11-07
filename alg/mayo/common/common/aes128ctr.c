// SPDX-License-Identifier: Apache-2.0 and MIT and Public Domain

#ifdef ENABLE_AESNI

#include <mem.h>
#include <stdint.h>
#include <string.h>
#include <tmmintrin.h>
#include <wmmintrin.h>
#ifdef MAYO_VAES_CTR
#include <immintrin.h>
#endif

// Adapted from liboqs/src/common/aes which in turn takes it from:
// crypto_core/aes128ncrypt/dolbeau/aesenc-int
// (https://bench.cr.yp.to/supercop.html)
static inline void aes128ni_setkey_encrypt(const unsigned char *key,
        __m128i rkeys[11]) {
    __m128i key0 = _mm_loadu_si128((const __m128i *)(key + 0));
    __m128i temp0, temp1, temp4;
    int idx = 0;

    temp0 = key0;

#define BLOCK1(IMM)                                                            \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                               \
  rkeys[idx++] = temp0;                                                        \
  temp4 = _mm_slli_si128(temp0, 4);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp4 = _mm_slli_si128(temp0, 8);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp1 = _mm_shuffle_epi32(temp1, 0xff);                                      \
  temp0 = _mm_xor_si128(temp0, temp1)

    BLOCK1(0x01);
    BLOCK1(0x02);
    BLOCK1(0x04);
    BLOCK1(0x08);
    BLOCK1(0x10);
    BLOCK1(0x20);
    BLOCK1(0x40);
    BLOCK1(0x80);
    BLOCK1(0x1b);
    BLOCK1(0x36);
    rkeys[idx++] = temp0;
}

void oqs_aes128_load_schedule_ni(const uint8_t *key, void **_schedule) {
    *_schedule = malloc(11 * sizeof(__m128i));
    // assert(*_schedule != NULL);
    __m128i *schedule = (__m128i *)*_schedule;
    aes128ni_setkey_encrypt(key, schedule);
}

void oqs_aes128_free_schedule_ni(void *schedule) {
    if (schedule != NULL) {
        mayo_secure_free(schedule, 11 * sizeof(__m128i));
    }
}

// Single encryption
static inline void aes128ni_encrypt(const __m128i rkeys[11], __m128i nv,
                                    unsigned char *out) {
    __m128i temp = _mm_xor_si128(nv, rkeys[0]);
    temp = _mm_aesenc_si128(temp, rkeys[1]);
    temp = _mm_aesenc_si128(temp, rkeys[2]);
    temp = _mm_aesenc_si128(temp, rkeys[3]);
    temp = _mm_aesenc_si128(temp, rkeys[4]);
    temp = _mm_aesenc_si128(temp, rkeys[5]);
    temp = _mm_aesenc_si128(temp, rkeys[6]);
    temp = _mm_aesenc_si128(temp, rkeys[7]);
    temp = _mm_aesenc_si128(temp, rkeys[8]);
    temp = _mm_aesenc_si128(temp, rkeys[9]);
    temp = _mm_aesenclast_si128(temp, rkeys[10]);
    _mm_storeu_si128((__m128i *)(out), temp);
}

// 4x interleaved encryption
static inline void aes128ni_encrypt_x4(const __m128i rkeys[11], __m128i n0,
                                       __m128i n1, __m128i n2, __m128i n3,
                                       unsigned char *out) {
    __m128i temp0 = _mm_xor_si128(n0, rkeys[0]);
    __m128i temp1 = _mm_xor_si128(n1, rkeys[0]);
    __m128i temp2 = _mm_xor_si128(n2, rkeys[0]);
    __m128i temp3 = _mm_xor_si128(n3, rkeys[0]);

#define AESNENCX4(IDX)                                                         \
  temp0 = _mm_aesenc_si128(temp0, rkeys[IDX]);                                 \
  temp1 = _mm_aesenc_si128(temp1, rkeys[IDX]);                                 \
  temp2 = _mm_aesenc_si128(temp2, rkeys[IDX]);                                 \
  temp3 = _mm_aesenc_si128(temp3, rkeys[IDX])

    AESNENCX4(1);
    AESNENCX4(2);
    AESNENCX4(3);
    AESNENCX4(4);
    AESNENCX4(5);
    AESNENCX4(6);
    AESNENCX4(7);
    AESNENCX4(8);
    AESNENCX4(9);

    temp0 = _mm_aesenclast_si128(temp0, rkeys[10]);
    temp1 = _mm_aesenclast_si128(temp1, rkeys[10]);
    temp2 = _mm_aesenclast_si128(temp2, rkeys[10]);
    temp3 = _mm_aesenclast_si128(temp3, rkeys[10]);

    _mm_storeu_si128((__m128i *)(out + 0), temp0);
    _mm_storeu_si128((__m128i *)(out + 16), temp1);
    _mm_storeu_si128((__m128i *)(out + 32), temp2);
    _mm_storeu_si128((__m128i *)(out + 48), temp3);
}

#ifdef MAYO_VAES_CTR
static inline void aes128ni_encrypt_x8_vaes(
    const __m256i rkeys[11], __m256i temp01, __m256i temp23,
    __m256i temp45, __m256i temp67, unsigned char *out) {
    temp01 ^= rkeys[0];
    temp23 ^= rkeys[0];
    temp45 ^= rkeys[0];
    temp67 ^= rkeys[0];

#define VAESENCX8(IDX)                                                        \
    temp01 = _mm256_aesenc_epi128(temp01, rkeys[IDX]);                        \
    temp23 = _mm256_aesenc_epi128(temp23, rkeys[IDX]);                        \
    temp45 = _mm256_aesenc_epi128(temp45, rkeys[IDX]);                        \
    temp67 = _mm256_aesenc_epi128(temp67, rkeys[IDX])

    VAESENCX8(1);
    VAESENCX8(2);
    VAESENCX8(3);
    VAESENCX8(4);
    VAESENCX8(5);
    VAESENCX8(6);
    VAESENCX8(7);
    VAESENCX8(8);
    VAESENCX8(9);

    temp01 = _mm256_aesenclast_epi128(temp01, rkeys[10]);
    temp23 = _mm256_aesenclast_epi128(temp23, rkeys[10]);
    temp45 = _mm256_aesenclast_epi128(temp45, rkeys[10]);
    temp67 = _mm256_aesenclast_epi128(temp67, rkeys[10]);
    _mm256_storeu_si256((__m256i *)(out + 0), temp01);
    _mm256_storeu_si256((__m256i *)(out + 32), temp23);
    _mm256_storeu_si256((__m256i *)(out + 64), temp45);
    _mm256_storeu_si256((__m256i *)(out + 96), temp67);
#undef VAESENCX8
}
#endif

// Not for general use: IV = 0, nonce = 0
static void oqs_aes128_ctr_enc_sch_ni(const void *schedule, uint8_t *out,
                                      size_t out_len) {
    __m128i mask =
        _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
    __m128i block = _mm_set_epi64x(0, 0);
    // block = _mm_xor_si128(block, block); // set to zero

#ifdef MAYO_VAES_CTR
    const __m128i *rkeys = (const __m128i *)schedule;
    __m256i vrkeys[11];
    for (size_t i = 0; i < 11; ++i) {
        vrkeys[i] = _mm256_broadcastsi128_si256(rkeys[i]);
    }

    const __m256i mask256 = _mm256_broadcastsi128_si256(mask);
    const __m256i counter_step = _mm256_set_epi64x(8, 0, 8, 0);
    __m256i counter01 = _mm256_set_epi64x(1, 0, 0, 0);
    __m256i counter23 = _mm256_set_epi64x(3, 0, 2, 0);
    __m256i counter45 = _mm256_set_epi64x(5, 0, 4, 0);
    __m256i counter67 = _mm256_set_epi64x(7, 0, 6, 0);
    while (out_len >= 128) {
        aes128ni_encrypt_x8_vaes(
            vrkeys, _mm256_shuffle_epi8(counter01, mask256),
            _mm256_shuffle_epi8(counter23, mask256),
            _mm256_shuffle_epi8(counter45, mask256),
            _mm256_shuffle_epi8(counter67, mask256), out);
        counter01 += counter_step;
        counter23 += counter_step;
        counter45 += counter_step;
        counter67 += counter_step;
        out += 128;
        out_len -= 128;
    }
    block = _mm256_castsi256_si128(_mm256_shuffle_epi8(counter01, mask256));
#endif
    while (out_len >= 64) {
        __m128i nv0 = block;
        __m128i nv1 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                          mask);
        __m128i nv2 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(2, 0)),
                          mask);
        __m128i nv3 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(3, 0)),
                          mask);
        aes128ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(4, 0)),
                    mask);
        out += 64;
        out_len -= 64;
    }
    while (out_len >= 16) {
        aes128ni_encrypt(schedule, block, out);
        out += 16;
        out_len -= 16;
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                    mask);
    }
    if (out_len > 0) {
        uint8_t tmp[16];
        aes128ni_encrypt(schedule, block, tmp);
        memcpy(out, tmp, out_len);
    }
}

int AES_128_CTR_NI(unsigned char *output, size_t outputByteLen,
                   const unsigned char *input, size_t inputByteLen) {
#ifdef MAYO_VAES_CTR
    __m128i schedule[11];
    (void)inputByteLen;
    aes128ni_setkey_encrypt(input, schedule);
    oqs_aes128_ctr_enc_sch_ni(schedule, output, outputByteLen);
    mayo_secure_clear(schedule, sizeof(schedule));
#else
    void *schedule = NULL;
    oqs_aes128_load_schedule_ni(input, &schedule);
    oqs_aes128_ctr_enc_sch_ni(schedule, output, outputByteLen);
    oqs_aes128_free_schedule_ni(schedule);
#endif
    return (int)outputByteLen;
}

// 4-Round AES...

// From crypto_core/aes128ncrypt/dolbeau/aesenc-int
static inline void aes128r4ni_setkey_encrypt(const unsigned char *key,
        __m128i rkeys[5]) {
    __m128i key0 = _mm_loadu_si128((const __m128i *)(key + 0));
    __m128i temp0, temp1, temp4;
    int idx = 0;

    temp0 = key0;

    /* blockshift-based block by Cedric Bourrasset */
#define BLOCK1(IMM)                                                            \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                               \
  rkeys[idx++] = temp0;                                                        \
  temp4 = _mm_slli_si128(temp0, 4);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp4 = _mm_slli_si128(temp0, 8);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp1 = _mm_shuffle_epi32(temp1, 0xff);                                      \
  temp0 = _mm_xor_si128(temp0, temp1)

    BLOCK1(0x01);
    BLOCK1(0x02);
    BLOCK1(0x04);
    BLOCK1(0x08);
    rkeys[idx++] = temp0;
}

void oqs_aes128r4_load_schedule_ni(const uint8_t *key, void **_schedule) {
    *_schedule = malloc(5 * sizeof(__m128i));
    // assert(*_schedule != NULL);
    __m128i *schedule = (__m128i *)*_schedule;
    aes128r4ni_setkey_encrypt(key, schedule);
}

void oqs_aes128r4_free_schedule_ni(void *schedule) {
    if (schedule != NULL) {
        mayo_secure_free(schedule, 5 * sizeof(__m128i));
    }
}

// Single encryption
static inline void aes128r4ni_encrypt(const __m128i rkeys[5], __m128i nv,
                                      unsigned char *out) {
    __m128i temp = _mm_xor_si128(nv, rkeys[0]);
    temp = _mm_aesenc_si128(temp, rkeys[1]);
    temp = _mm_aesenc_si128(temp, rkeys[2]);
    temp = _mm_aesenc_si128(temp, rkeys[3]);
    temp = _mm_aesenclast_si128(temp, rkeys[4]);
    _mm_storeu_si128((__m128i *)(out), temp);
}

// 4x interleaved encryption
static inline void aes128r4ni_encrypt_x4(const __m128i rkeys[5], __m128i n0,
        __m128i n1, __m128i n2, __m128i n3,
        unsigned char *out) {
    __m128i temp0 = _mm_xor_si128(n0, rkeys[0]);
    __m128i temp1 = _mm_xor_si128(n1, rkeys[0]);
    __m128i temp2 = _mm_xor_si128(n2, rkeys[0]);
    __m128i temp3 = _mm_xor_si128(n3, rkeys[0]);

#define AESNENCX4(IDX)                                                         \
  temp0 = _mm_aesenc_si128(temp0, rkeys[IDX]);                                 \
  temp1 = _mm_aesenc_si128(temp1, rkeys[IDX]);                                 \
  temp2 = _mm_aesenc_si128(temp2, rkeys[IDX]);                                 \
  temp3 = _mm_aesenc_si128(temp3, rkeys[IDX])

    AESNENCX4(1);
    AESNENCX4(2);
    AESNENCX4(3);

    temp0 = _mm_aesenclast_si128(temp0, rkeys[4]);
    temp1 = _mm_aesenclast_si128(temp1, rkeys[4]);
    temp2 = _mm_aesenclast_si128(temp2, rkeys[4]);
    temp3 = _mm_aesenclast_si128(temp3, rkeys[4]);

    _mm_storeu_si128((__m128i *)(out + 0), temp0);
    _mm_storeu_si128((__m128i *)(out + 16), temp1);
    _mm_storeu_si128((__m128i *)(out + 32), temp2);
    _mm_storeu_si128((__m128i *)(out + 48), temp3);
}

// Not for general use: IV = 0, nonce = 0
static void oqs_aes128r4_ctr_enc_sch_ni(const void *schedule, uint8_t *out,
                                        size_t out_len) {
    __m128i mask =
        _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
    __m128i block = _mm_set_epi64x(0, 0);

    while (out_len >= 64) {
        __m128i nv0 = block;
        __m128i nv1 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                          mask);
        __m128i nv2 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(2, 0)),
                          mask);
        __m128i nv3 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(3, 0)),
                          mask);
        aes128r4ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(4, 0)),
                    mask);
        out += 64;
        out_len -= 64;
    }
    while (out_len >= 16) {
        aes128r4ni_encrypt(schedule, block, out);
        out += 16;
        out_len -= 16;
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                    mask);
    }
    if (out_len > 0) {
        uint8_t tmp[16];
        aes128r4ni_encrypt(schedule, block, tmp);
        memcpy(out, tmp, out_len);
    }
}

int AES_128_CTR_4R_NI(unsigned char *output, size_t outputByteLen,
                      const unsigned char *input, size_t inputByteLen) {
    void *schedule = NULL;
    oqs_aes128r4_load_schedule_ni(input, &schedule);
    oqs_aes128r4_ctr_enc_sch_ni(schedule, output, outputByteLen);
    oqs_aes128r4_free_schedule_ni(schedule);
    return (int)outputByteLen;
}
#endif

