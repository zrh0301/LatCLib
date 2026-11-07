#include "gfext.h"

#include <string.h>

/*
 * GF(256) 的 GF(16) 塔分解层（正文 4.5 节）。
 *
 *   GF(256) = GF(16)[y] / (y^2 + y + c),  c = 0x08
 *   塔坐标字节 = (a1 << 4) | a0，低 nibble 是 1 的系数、高 nibble 是 y 的系数
 *   phi: Rijndael 0x11B 字节 -> 塔坐标;  psi: 塔坐标 -> 0x11B
 *
 * 接口分两层，用来把“乘法成本”和“表示转换成本”分开计时：
 *   *_array        输入/输出保持塔坐标，内核里没有 phi/psi；
 *   *_11b_array    输入/输出都是 0x11B 字节，转换在调用边界发生。
 *
 * 乘法恒等式（Karatsuba，c 在 high 中消失）：
 *   low  = m1 ^ (c * m2),   high = m1 ^ m3
 *   m1 = a0*b0, m2 = a1*b1, m3 = (a0^a1)*(b0^b1)
 * 等价地 high = m2 ^ (a0*b1) ^ (a1*b0)。SIMD 走四个叉积这一支，因为 packed
 * 内核一次 pass 就能同时给出“直配”和“交叉”两路结果。
 *
 * 静态数据 = 4 张 16 项坐标表 + 1 张 16 项逆元表（80 byte）；其余向量常量
 * （含 16 张常数乘表）在首次调用时惰性建立，因此本层不依赖 pqc_kernel_gf16_tables_init()。
 */

#define TOWER_C 0x08u

/* 0x11B -> 塔坐标，按 nibble 查表后异或（4.5.1 节，已穷举校验）。 */
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

/* 16 项逆元表：约定 INV[0] = 0，零元素自动落到零，无需额外非零掩码。 */
static const uint8_t TOWER_INV16[16] = {
    0x00, 0x01, 0x09, 0x0e, 0x0d, 0x0b, 0x07, 0x06,
    0x0f, 0x02, 0x0c, 0x05, 0x0a, 0x04, 0x03, 0x08
};

#if defined(_MSC_VER)
#define TOWER_ALIGN(N) __declspec(align(N))
#else
#define TOWER_ALIGN(N) __attribute__((aligned(N)))
#endif

static pqc_kernel_gf16_vec vec_phi_lo;
static pqc_kernel_gf16_vec vec_phi_hi;
static pqc_kernel_gf16_vec vec_psi_lo;
static pqc_kernel_gf16_vec vec_psi_hi;
static pqc_kernel_gf16_vec vec_cmul_c;
static pqc_kernel_gf16_vec vec_inv16;
static pqc_kernel_gf16_vec vec_sqr16;
static pqc_kernel_gf16_vec vec_mk[16]; /* vec_mk[k][x] = x * k */
static int tower_ready = 0;

/* 本层自带的 GF(16) 乘法，避免依赖 pqc_kernel_gf16_tables_init()。 */
static uint8_t pqc_kernel_gf16_mul_local(uint8_t a, uint8_t b) {
    uint8_t result = 0;
    a &= 0x0fu;
    b &= 0x0fu;
    for (unsigned i = 0; i < 4; ++i) {
        result ^= (uint8_t)(a & (uint8_t)-(int)(b & 1u));
        const uint8_t carry = (uint8_t)(a & 0x08u);
        a = (uint8_t)(a << 1);
        if (carry) a ^= 0x13u;
        b >>= 1;
    }
    return (uint8_t)(result & 0x0fu);
}

static uint8_t pqc_kernel_gf16_cmul_c_local(uint8_t x) { return pqc_kernel_gf16_mul_local(x, TOWER_C); }

static void build_vec16(pqc_kernel_gf16_vec *dst, const uint8_t tab16[16]) {
    TOWER_ALIGN(32) uint8_t tmp[GF16_VEC_BYTES];
    for (size_t i = 0; i < GF16_VEC_BYTES; ++i) tmp[i] = tab16[i & 0x0fu];
    *dst = pqc_kernel_gf16_vec_load(tmp);
}

void pqc_kernel_gf256_tower_init(void) {
    if (tower_ready) return;
    /* 塔层复用 pqc_kernel_gf16_bitplane 的 packed 内核，而该内核依赖 shuffle 家族
     * 的全局向量常量（pqc_kernel_gf16_mul2/4/8_vec 等），所以这里必须先初始化它们。 */
    pqc_kernel_gf16_tables_init();
    uint8_t tab[16];
    build_vec16(&vec_phi_lo, TOWER_PHI_LO);
    build_vec16(&vec_phi_hi, TOWER_PHI_HI);
    build_vec16(&vec_psi_lo, TOWER_PSI_LO);
    build_vec16(&vec_psi_hi, TOWER_PSI_HI);
    build_vec16(&vec_inv16, TOWER_INV16);
    for (unsigned x = 0; x < 16; ++x) tab[x] = pqc_kernel_gf16_cmul_c_local((uint8_t)x);
    build_vec16(&vec_cmul_c, tab);
    for (unsigned x = 0; x < 16; ++x) tab[x] = pqc_kernel_gf16_mul_local((uint8_t)x, (uint8_t)x);
    build_vec16(&vec_sqr16, tab);
    for (unsigned k = 0; k < 16; ++k) {
        for (unsigned x = 0; x < 16; ++x) tab[x] = pqc_kernel_gf16_mul_local((uint8_t)x, (uint8_t)k);
        build_vec16(&vec_mk[k], tab);
    }
    tower_ready = 1;
}

/* ---------------------------------------------------------------- 标量路径 */

uint8_t pqc_kernel_gf256_tower_phi(uint8_t a) {
    return (uint8_t)(TOWER_PHI_LO[a & 0x0fu] ^ TOWER_PHI_HI[a >> 4]);
}

uint8_t pqc_kernel_gf256_tower_psi(uint8_t a) {
    return (uint8_t)(TOWER_PSI_LO[a & 0x0fu] ^ TOWER_PSI_HI[a >> 4]);
}

uint8_t pqc_kernel_gf256_tower_mul(uint8_t a, uint8_t b) {
    const uint8_t a0 = (uint8_t)(a & 0x0fu), a1 = (uint8_t)(a >> 4);
    const uint8_t b0 = (uint8_t)(b & 0x0fu), b1 = (uint8_t)(b >> 4);
    const uint8_t m1 = pqc_kernel_gf16_mul_local(a0, b0);
    const uint8_t m2 = pqc_kernel_gf16_mul_local(a1, b1);
    const uint8_t m3 = pqc_kernel_gf16_mul_local((uint8_t)(a0 ^ a1), (uint8_t)(b0 ^ b1));
    const uint8_t low = (uint8_t)(m1 ^ pqc_kernel_gf16_cmul_c_local(m2));
    const uint8_t high = (uint8_t)(m1 ^ m3);
    return (uint8_t)(low | (uint8_t)(high << 4));
}

uint8_t pqc_kernel_gf256_tower_mul_11b(uint8_t a, uint8_t b) {
    return pqc_kernel_gf256_tower_psi(pqc_kernel_gf256_tower_mul(pqc_kernel_gf256_tower_phi(a), pqc_kernel_gf256_tower_phi(b)));
}

/* 平方：(a0 + a1 y)^2 = (a0^2 + c a1^2) + a1^2 y，只要 2 次 GF(16) 平方。 */
uint8_t pqc_kernel_gf256_tower_square(uint8_t a) {
    const uint8_t a0 = (uint8_t)(a & 0x0fu), a1 = (uint8_t)(a >> 4);
    const uint8_t m1 = pqc_kernel_gf16_mul_local(a0, a0);
    const uint8_t m2 = pqc_kernel_gf16_mul_local(a1, a1);
    const uint8_t low = (uint8_t)(m1 ^ pqc_kernel_gf16_cmul_c_local(m2));
    return (uint8_t)(low | (uint8_t)(m2 << 4));
}

uint8_t pqc_kernel_gf256_tower_square_11b(uint8_t a) {
    return pqc_kernel_gf256_tower_psi(pqc_kernel_gf256_tower_square(pqc_kernel_gf256_tower_phi(a)));
}

/* Frobenius：sigma(a) = a^16 = (a0 ^ a1) + a1 y，一次 nibble 异或。 */
uint8_t pqc_kernel_gf256_tower_frobenius(uint8_t a) {
    const uint8_t a0 = (uint8_t)(a & 0x0fu), a1 = (uint8_t)(a >> 4);
    return (uint8_t)((a0 ^ a1) | (uint8_t)(a1 << 4));
}

uint8_t pqc_kernel_gf256_tower_frobenius_11b(uint8_t a) {
    return pqc_kernel_gf256_tower_psi(pqc_kernel_gf256_tower_frobenius(pqc_kernel_gf256_tower_phi(a)));
}

/* 求逆：a^-1 = sigma(a) * N(a)^-1，N(a) = a0^2 + a0 a1 + c a1^2 落在 GF(16) 中。
 * 因此不需要 256 项表，只要 3.9 节的 16 项逆元表。 */
uint8_t pqc_kernel_gf256_tower_inv(uint8_t a) {
    const uint8_t a0 = (uint8_t)(a & 0x0fu), a1 = (uint8_t)(a >> 4);
    const uint8_t s0 = (uint8_t)(a0 ^ a1);
    const uint8_t s1 = a1;
    const uint8_t m1 = pqc_kernel_gf16_mul_local(a0, s0);
    const uint8_t m2 = pqc_kernel_gf16_mul_local(s1, s1);
    const uint8_t norm = (uint8_t)(m1 ^ pqc_kernel_gf16_cmul_c_local(m2));
    const uint8_t ninv = TOWER_INV16[norm];
    const uint8_t low = pqc_kernel_gf16_mul_local(s0, ninv);
    const uint8_t high = pqc_kernel_gf16_mul_local(s1, ninv);
    return (uint8_t)(low | (uint8_t)(high << 4));
}

uint8_t pqc_kernel_gf256_tower_inv_11b(uint8_t a) {
    return pqc_kernel_gf256_tower_psi(pqc_kernel_gf256_tower_inv(pqc_kernel_gf256_tower_phi(a)));
}

uint8_t pqc_kernel_gf256_tower_pow(uint8_t a, uint32_t e) {
    uint8_t result = 0x01u; /* 塔坐标下的乘法单位元仍是 0x01 */
    uint8_t base = a;
    while (e) {
        if (e & 1u) result = pqc_kernel_gf256_tower_mul(result, base);
        e >>= 1;
        if (e) base = pqc_kernel_gf256_tower_square(base);
    }
    return result;
}

uint8_t pqc_kernel_gf256_tower_pow_11b(uint8_t a, uint32_t e) {
    return pqc_kernel_gf256_tower_psi(pqc_kernel_gf256_tower_pow(pqc_kernel_gf256_tower_phi(a), e));
}

/* 常数乘（塔坐标）：k = k0 + k1 y 时
 *   low  = k0*a0 ^ (c k1)*a1
 *   high = k1*a0 ^ (k0^k1)*a1
 * 四路都是 GF(16) 常数乘，因此可以全部走 shuffle 表。 */
static uint8_t tower_const_mul_byte(uint8_t a, uint8_t k0, uint8_t k1, uint8_t ck1) {
    const uint8_t a0 = (uint8_t)(a & 0x0fu), a1 = (uint8_t)(a >> 4);
    const uint8_t low = (uint8_t)(pqc_kernel_gf16_mul_local(k0, a0) ^ pqc_kernel_gf16_mul_local(ck1, a1));
    const uint8_t high = (uint8_t)(pqc_kernel_gf16_mul_local(k1, a0)
                                   ^ pqc_kernel_gf16_mul_local((uint8_t)(k0 ^ k1), a1));
    return (uint8_t)(low | (uint8_t)(high << 4));
}

/* ---------------------------------------------------------- SIMD 基本内核 */

static inline pqc_kernel_gf16_vec tower_phi_block(pqc_kernel_gf16_vec a) {
    const pqc_kernel_gf16_vec lo = pqc_kernel_gf16_vec_shuffle(vec_phi_lo, pqc_kernel_gf16_vec_mask_u8(a, 0x0f));
    const pqc_kernel_gf16_vec hi = pqc_kernel_gf16_vec_shuffle(vec_phi_hi, pqc_kernel_gf16_vec_srli4(a));
    return pqc_kernel_gf16_vec_xor(lo, hi);
}

static inline pqc_kernel_gf16_vec tower_psi_block(pqc_kernel_gf16_vec a) {
    const pqc_kernel_gf16_vec lo = pqc_kernel_gf16_vec_shuffle(vec_psi_lo, pqc_kernel_gf16_vec_mask_u8(a, 0x0f));
    const pqc_kernel_gf16_vec hi = pqc_kernel_gf16_vec_shuffle(vec_psi_hi, pqc_kernel_gf16_vec_srli4(a));
    return pqc_kernel_gf16_vec_xor(lo, hi);
}

/* 塔坐标乘法：1 次 packed 双路乘（直配 + 交叉）＋ 1 次 CMUL_C 查表。
 * P = [m2 | m1]，X = [a1*b0 | a0*b1]。 */
static inline pqc_kernel_gf16_vec tower_mul_block(pqc_kernel_gf16_vec A, pqc_kernel_gf16_vec B) {
    pqc_kernel_gf16_vec P, X;
    pqc_kernel_gf16_vec_mul_packed_ab_dual_packed(A, B, &P, &X);
    const pqc_kernel_gf16_vec m1 = pqc_kernel_gf16_vec_mask_u8(P, 0x0f);
    const pqc_kernel_gf16_vec m2 = pqc_kernel_gf16_vec_srli4(P);
    const pqc_kernel_gf16_vec low = pqc_kernel_gf16_vec_xor(m1, pqc_kernel_gf16_vec_shuffle(vec_cmul_c, m2));
    const pqc_kernel_gf16_vec high = pqc_kernel_gf16_vec_xor(m2, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_mask_u8(X, 0x0f),
                                                        pqc_kernel_gf16_vec_srli4(X)));
    return pqc_kernel_gf16_vec_or(low, pqc_kernel_gf16_vec_slli4(high));
}

/* 乘加：low / high 两个坐标分别累加，省掉每一轮的重新打包。 */
static inline void tower_mulacc_block(pqc_kernel_gf16_vec A, pqc_kernel_gf16_vec B,
                                      pqc_kernel_gf16_vec *lo_acc, pqc_kernel_gf16_vec *hi_acc) {
    pqc_kernel_gf16_vec P, X;
    pqc_kernel_gf16_vec_mul_packed_ab_dual_packed(A, B, &P, &X);
    const pqc_kernel_gf16_vec m1 = pqc_kernel_gf16_vec_mask_u8(P, 0x0f);
    const pqc_kernel_gf16_vec m2 = pqc_kernel_gf16_vec_srli4(P);
    const pqc_kernel_gf16_vec low = pqc_kernel_gf16_vec_xor(m1, pqc_kernel_gf16_vec_shuffle(vec_cmul_c, m2));
    const pqc_kernel_gf16_vec high = pqc_kernel_gf16_vec_xor(m2, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_mask_u8(X, 0x0f),
                                                        pqc_kernel_gf16_vec_srli4(X)));
    *lo_acc = pqc_kernel_gf16_vec_xor(*lo_acc, low);
    *hi_acc = pqc_kernel_gf16_vec_xor(*hi_acc, high);
}

static inline pqc_kernel_gf16_vec tower_square_block(pqc_kernel_gf16_vec A) {
    const pqc_kernel_gf16_vec m1 = pqc_kernel_gf16_vec_shuffle(vec_sqr16, pqc_kernel_gf16_vec_mask_u8(A, 0x0f));
    const pqc_kernel_gf16_vec m2 = pqc_kernel_gf16_vec_shuffle(vec_sqr16, pqc_kernel_gf16_vec_srli4(A));
    const pqc_kernel_gf16_vec low = pqc_kernel_gf16_vec_xor(m1, pqc_kernel_gf16_vec_shuffle(vec_cmul_c, m2));
    return pqc_kernel_gf16_vec_or(low, pqc_kernel_gf16_vec_slli4(m2));
}

static inline pqc_kernel_gf16_vec tower_frobenius_block(pqc_kernel_gf16_vec A) {
    const pqc_kernel_gf16_vec a0 = pqc_kernel_gf16_vec_mask_u8(A, 0x0f);
    const pqc_kernel_gf16_vec a1 = pqc_kernel_gf16_vec_srli4(A);
    return pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_xor(a0, a1), pqc_kernel_gf16_vec_slli4(a1));
}

static inline pqc_kernel_gf16_vec tower_inv_block(pqc_kernel_gf16_vec A) {
    const pqc_kernel_gf16_vec a0 = pqc_kernel_gf16_vec_mask_u8(A, 0x0f);
    const pqc_kernel_gf16_vec a1 = pqc_kernel_gf16_vec_srli4(A);
    const pqc_kernel_gf16_vec sigma = pqc_kernel_gf16_vec_or(pqc_kernel_gf16_vec_xor(a0, a1), pqc_kernel_gf16_vec_slli4(a1));
    /* P = [a1*a1 | a0*(a0^a1)]，低 nibble 正好是范数 N(a) 的一项。 */
    const pqc_kernel_gf16_vec P = pqc_kernel_gf16_vec_mul_packed_ab_straight(A, sigma);
    const pqc_kernel_gf16_vec m1 = pqc_kernel_gf16_vec_mask_u8(P, 0x0f);
    const pqc_kernel_gf16_vec m2 = pqc_kernel_gf16_vec_srli4(P);
    const pqc_kernel_gf16_vec norm = pqc_kernel_gf16_vec_xor(m1, pqc_kernel_gf16_vec_shuffle(vec_cmul_c, m2));
    const pqc_kernel_gf16_vec ninv = pqc_kernel_gf16_vec_shuffle(vec_inv16, norm);
    const pqc_kernel_gf16_vec ndup = pqc_kernel_gf16_vec_or(ninv, pqc_kernel_gf16_vec_slli4(ninv));
    return pqc_kernel_gf16_vec_mul_packed_ab_straight(sigma, ndup);
}

/* 常数乘：4 次 16 项 shuffle + 3 次异或，表在循环外取好。 */
static inline pqc_kernel_gf16_vec tower_const_mul_block(pqc_kernel_gf16_vec A, uint8_t k0, uint8_t k1,
                                             uint8_t ck1) {
    const pqc_kernel_gf16_vec a0 = pqc_kernel_gf16_vec_mask_u8(A, 0x0f);
    const pqc_kernel_gf16_vec a1 = pqc_kernel_gf16_vec_srli4(A);
    const pqc_kernel_gf16_vec low = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_shuffle(vec_mk[k0 & 0x0f], a0),
                                      pqc_kernel_gf16_vec_shuffle(vec_mk[ck1 & 0x0f], a1));
    const pqc_kernel_gf16_vec high = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_shuffle(vec_mk[k1 & 0x0f], a0),
                                       pqc_kernel_gf16_vec_shuffle(vec_mk[(k0 ^ k1) & 0x0f], a1));
    return pqc_kernel_gf16_vec_or(low, pqc_kernel_gf16_vec_slli4(high));
}

static uint8_t tower_reduce_xor(pqc_kernel_gf16_vec v) {
    TOWER_ALIGN(32) uint8_t tmp[GF16_VEC_BYTES];
    pqc_kernel_gf16_vec_store(tmp, v);
    uint8_t r = 0;
    for (size_t i = 0; i < GF16_VEC_BYTES; ++i) r ^= tmp[i];
    return r;
}

/* --------------------------------------------------------- 坐标转换（SIMD） */

void pqc_kernel_gf256_tower_phi_array(uint8_t *dst, const uint8_t *src, size_t n) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(dst + i, tower_phi_block(pqc_kernel_gf16_vec_load(src + i)));
    for (; i < n; ++i) dst[i] = pqc_kernel_gf256_tower_phi(src[i]);
}

void pqc_kernel_gf256_tower_psi_array(uint8_t *dst, const uint8_t *src, size_t n) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(dst + i, tower_psi_block(pqc_kernel_gf16_vec_load(src + i)));
    for (; i < n; ++i) dst[i] = pqc_kernel_gf256_tower_psi(src[i]);
}

/* ------------------------------------------------------- 标量参考数组接口 */

#define TOWER_SCALAR_ARRAY(name, op)                                          \
    void name(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {   \
        for (size_t i = 0; i < n; ++i) out[i] = (op)(a[i], b[i]);             \
    }

TOWER_SCALAR_ARRAY(pqc_kernel_gf256_tower_mul_ref_array, pqc_kernel_gf256_tower_mul)
TOWER_SCALAR_ARRAY(pqc_kernel_gf256_tower_mul_11b_ref_array, pqc_kernel_gf256_tower_mul_11b)

void pqc_kernel_gf256_tower_mulacc_ref_array(uint8_t *acc, const uint8_t *a, const uint8_t *b,
                                  size_t n) {
    for (size_t i = 0; i < n; ++i) acc[i] ^= pqc_kernel_gf256_tower_mul(a[i], b[i]);
}

void pqc_kernel_gf256_tower_mulacc_11b_ref_array(uint8_t *acc, const uint8_t *a, const uint8_t *b,
                                      size_t n) {
    for (size_t i = 0; i < n; ++i) acc[i] ^= pqc_kernel_gf256_tower_mul_11b(a[i], b[i]);
}

uint8_t pqc_kernel_gf256_tower_dot_ref_array(const uint8_t *a, const uint8_t *b, size_t n) {
    uint8_t r = 0;
    for (size_t i = 0; i < n; ++i) r ^= pqc_kernel_gf256_tower_mul(a[i], b[i]);
    return r;
}

uint8_t pqc_kernel_gf256_tower_dot_11b_ref_array(const uint8_t *a, const uint8_t *b, size_t n) {
    uint8_t r = 0;
    for (size_t i = 0; i < n; ++i) r ^= pqc_kernel_gf256_tower_mul_11b(a[i], b[i]);
    return r;
}

void pqc_kernel_gf256_tower_const_mul_ref_array(uint8_t *out, const uint8_t *a, uint8_t k,
                                     size_t n) {
    const uint8_t k0 = (uint8_t)(k & 0x0fu), k1 = (uint8_t)(k >> 4);
    const uint8_t ck1 = pqc_kernel_gf16_cmul_c_local(k1);
    for (size_t i = 0; i < n; ++i) out[i] = tower_const_mul_byte(a[i], k0, k1, ck1);
}

void pqc_kernel_gf256_tower_const_mul_11b_ref_array(uint8_t *out, const uint8_t *a, uint8_t k,
                                         size_t n) {
    const uint8_t kt = pqc_kernel_gf256_tower_phi(k);
    const uint8_t k0 = (uint8_t)(kt & 0x0fu), k1 = (uint8_t)(kt >> 4);
    const uint8_t ck1 = pqc_kernel_gf16_cmul_c_local(k1);
    for (size_t i = 0; i < n; ++i) {
        out[i] = pqc_kernel_gf256_tower_psi(tower_const_mul_byte(pqc_kernel_gf256_tower_phi(a[i]), k0, k1, ck1));
    }
}

#define TOWER_SCALAR_UNARY_ARRAY(name, scalar_op)                             \
    void name(uint8_t *out, const uint8_t *a, size_t n) {                      \
        for (size_t i = 0; i < n; ++i) out[i] = (scalar_op)(a[i]);             \
    }

TOWER_SCALAR_UNARY_ARRAY(pqc_kernel_gf256_tower_square_ref_array, pqc_kernel_gf256_tower_square)
TOWER_SCALAR_UNARY_ARRAY(pqc_kernel_gf256_tower_square_11b_ref_array, pqc_kernel_gf256_tower_square_11b)
TOWER_SCALAR_UNARY_ARRAY(pqc_kernel_gf256_tower_frobenius_ref_array, pqc_kernel_gf256_tower_frobenius)
TOWER_SCALAR_UNARY_ARRAY(pqc_kernel_gf256_tower_frobenius_11b_ref_array, pqc_kernel_gf256_tower_frobenius_11b)
TOWER_SCALAR_UNARY_ARRAY(pqc_kernel_gf256_tower_inv_ref_array, pqc_kernel_gf256_tower_inv)
TOWER_SCALAR_UNARY_ARRAY(pqc_kernel_gf256_tower_inv_11b_ref_array, pqc_kernel_gf256_tower_inv_11b)

/* ---------------------------------------------------------- 尾部向量化处理
 * 长度不是向量宽度整数倍时（例如 MQOM 的 n=116），尾部若走标量塔内核，
 * 单字节成本比向量 lane 高一个数量级，会吃掉整条链的收益。
 * 这里把尾部补零成一个完整向量再过同一内核：塔内核把 0 映射到 0，
 * 多出来的 lane 结果直接丢弃，只写回 r 个字节。
 */
static pqc_kernel_gf16_vec tower_tail_load(const uint8_t *p, size_t r) {
    uint8_t tmp[GF16_VEC_BYTES];
    memset(tmp, 0, sizeof(tmp));
    memcpy(tmp, p, r);
    return pqc_kernel_gf16_vec_load(tmp);
}

static void tower_tail_store(uint8_t *out, pqc_kernel_gf16_vec v, size_t r) {
    uint8_t tmp[GF16_VEC_BYTES];
    pqc_kernel_gf16_vec_store(tmp, v);
    memcpy(out, tmp, r);
}

/* ------------------------------------------------------------ SIMD 数组接口 */

void pqc_kernel_gf256_tower_mul_array(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(out + i, tower_mul_block(pqc_kernel_gf16_vec_load(a + i),
                                                pqc_kernel_gf16_vec_load(b + i)));
    if (i < n)
        tower_tail_store(out + i,
                         tower_mul_block(tower_tail_load(a + i, n - i),
                                         tower_tail_load(b + i, n - i)),
                         n - i);
}

void pqc_kernel_gf256_tower_mulacc_array(uint8_t *acc, const uint8_t *a, const uint8_t *b,
                              size_t n) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec lo_acc = pqc_kernel_gf16_vec_zero(), hi_acc = pqc_kernel_gf16_vec_zero();
        tower_mulacc_block(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i), &lo_acc, &hi_acc);
        const pqc_kernel_gf16_vec t = pqc_kernel_gf16_vec_or(lo_acc, pqc_kernel_gf16_vec_slli4(hi_acc));
        pqc_kernel_gf16_vec_store(acc + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(acc + i), t));
    }
    if (i < n) {
        pqc_kernel_gf16_vec lo_t = pqc_kernel_gf16_vec_zero(), hi_t = pqc_kernel_gf16_vec_zero();
        tower_mulacc_block(tower_tail_load(a + i, n - i),
                           tower_tail_load(b + i, n - i), &lo_t, &hi_t);
        const pqc_kernel_gf16_vec t = pqc_kernel_gf16_vec_or(lo_t, pqc_kernel_gf16_vec_slli4(hi_t));
        tower_tail_store(acc + i,
                         pqc_kernel_gf16_vec_xor(tower_tail_load(acc + i, n - i), t), n - i);
    }
}

uint8_t pqc_kernel_gf256_tower_dot_array(const uint8_t *a, const uint8_t *b, size_t n) {
    pqc_kernel_gf256_tower_init();
    pqc_kernel_gf16_vec lo_acc = pqc_kernel_gf16_vec_zero(), hi_acc = pqc_kernel_gf16_vec_zero();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        tower_mulacc_block(pqc_kernel_gf16_vec_load(a + i), pqc_kernel_gf16_vec_load(b + i), &lo_acc, &hi_acc);
    if (i < n)
        tower_mulacc_block(tower_tail_load(a + i, n - i),
                           tower_tail_load(b + i, n - i), &lo_acc, &hi_acc);
    return tower_reduce_xor(pqc_kernel_gf16_vec_or(lo_acc, pqc_kernel_gf16_vec_slli4(hi_acc)));
}

void pqc_kernel_gf256_tower_const_mul_array(uint8_t *out, const uint8_t *a, uint8_t k, size_t n) {
    pqc_kernel_gf256_tower_init();
    const uint8_t k0 = (uint8_t)(k & 0x0fu), k1 = (uint8_t)(k >> 4);
    const uint8_t ck1 = pqc_kernel_gf16_cmul_c_local(k1);
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        pqc_kernel_gf16_vec_store(out + i, tower_const_mul_block(pqc_kernel_gf16_vec_load(a + i), k0, k1, ck1));
    if (i < n)
        tower_tail_store(out + i,
                         tower_const_mul_block(tower_tail_load(a + i, n - i), k0, k1, ck1),
                         n - i);
}

#define TOWER_SIMD_UNARY_ARRAY(name, block_op, scalar_op, pre, post)          \
    void name(uint8_t *out, const uint8_t *a, size_t n) {                      \
        pqc_kernel_gf256_tower_init();                                                    \
        size_t i = 0;                                                          \
        for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)                   \
            pqc_kernel_gf16_vec_store(out + i, (post)((block_op)((pre)(                   \
                                        pqc_kernel_gf16_vec_load(a + i)))));              \
        if (i < n)                                                             \
            tower_tail_store(out + i, (post)((block_op)((pre)(                 \
                tower_tail_load(a + i, n - i)))), n - i);                      \
    }

static inline pqc_kernel_gf16_vec tower_ident_block(pqc_kernel_gf16_vec v) { return v; }

TOWER_SIMD_UNARY_ARRAY(pqc_kernel_gf256_tower_square_array, tower_square_block,
                       pqc_kernel_gf256_tower_square, tower_ident_block, tower_ident_block)
TOWER_SIMD_UNARY_ARRAY(pqc_kernel_gf256_tower_frobenius_array, tower_frobenius_block,
                       pqc_kernel_gf256_tower_frobenius, tower_ident_block, tower_ident_block)
TOWER_SIMD_UNARY_ARRAY(pqc_kernel_gf256_tower_inv_array, tower_inv_block,
                       pqc_kernel_gf256_tower_inv, tower_ident_block, tower_ident_block)
TOWER_SIMD_UNARY_ARRAY(pqc_kernel_gf256_tower_square_11b_array, tower_square_block,
                       pqc_kernel_gf256_tower_square_11b, tower_phi_block, tower_psi_block)
TOWER_SIMD_UNARY_ARRAY(pqc_kernel_gf256_tower_frobenius_11b_array, tower_frobenius_block,
                       pqc_kernel_gf256_tower_frobenius_11b, tower_phi_block, tower_psi_block)
TOWER_SIMD_UNARY_ARRAY(pqc_kernel_gf256_tower_inv_11b_array, tower_inv_block,
                       pqc_kernel_gf256_tower_inv_11b, tower_phi_block, tower_psi_block)

/* 平方-乘链的向量版：每个 block 内独立完成整条幂链。 */
static void tower_pow_array_impl(uint8_t *out, const uint8_t *a, uint32_t e, size_t n,
                                 int convert) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec base = pqc_kernel_gf16_vec_load(a + i);
        pqc_kernel_gf16_vec result = pqc_kernel_gf16_vec_set1(0x01);
        if (convert) base = tower_phi_block(base);
        uint32_t bits = e;
        while (bits) {
            if (bits & 1u) result = tower_mul_block(result, base);
            bits >>= 1;
            if (bits) base = tower_square_block(base);
        }
        pqc_kernel_gf16_vec_store(out + i, convert ? tower_psi_block(result) : result);
    }
    if (i < n) {
        pqc_kernel_gf16_vec base = tower_tail_load(a + i, n - i);
        pqc_kernel_gf16_vec result = pqc_kernel_gf16_vec_set1(0x01);
        uint32_t bits = e;
        if (convert) base = tower_phi_block(base);
        while (bits) {
            if (bits & 1u) result = tower_mul_block(result, base);
            bits >>= 1;
            if (bits) base = tower_square_block(base);
        }
        tower_tail_store(out + i, convert ? tower_psi_block(result) : result, n - i);
    }
}

void pqc_kernel_gf256_tower_pow_array(uint8_t *out, const uint8_t *a, uint32_t e, size_t n) {
    tower_pow_array_impl(out, a, e, n, 0);
}

void pqc_kernel_gf256_tower_pow_11b_array(uint8_t *out, const uint8_t *a, uint32_t e, size_t n) {
    tower_pow_array_impl(out, a, e, n, 1);
}

/* ---------------------------------------------------- 0x11B 表示数组接口 */

void pqc_kernel_gf256_tower_mul_11b_array(uint8_t *out, const uint8_t *a, const uint8_t *b,
                               size_t n) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf16_vec t = tower_mul_block(tower_phi_block(pqc_kernel_gf16_vec_load(a + i)),
                                           tower_phi_block(pqc_kernel_gf16_vec_load(b + i)));
        pqc_kernel_gf16_vec_store(out + i, tower_psi_block(t));
    }
    if (i < n)
        tower_tail_store(out + i,
                         tower_psi_block(tower_mul_block(
                             tower_phi_block(tower_tail_load(a + i, n - i)),
                             tower_phi_block(tower_tail_load(b + i, n - i)))),
                         n - i);
}

void pqc_kernel_gf256_tower_mulacc_11b_array(uint8_t *acc, const uint8_t *a, const uint8_t *b,
                                  size_t n) {
    pqc_kernel_gf256_tower_init();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        pqc_kernel_gf16_vec lo_acc = pqc_kernel_gf16_vec_zero(), hi_acc = pqc_kernel_gf16_vec_zero();
        tower_mulacc_block(tower_phi_block(pqc_kernel_gf16_vec_load(a + i)),
                           tower_phi_block(pqc_kernel_gf16_vec_load(b + i)), &lo_acc, &hi_acc);
        const pqc_kernel_gf16_vec t = pqc_kernel_gf16_vec_or(lo_acc, pqc_kernel_gf16_vec_slli4(hi_acc));
        pqc_kernel_gf16_vec_store(acc + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(acc + i),
                                             tower_psi_block(t)));
    }
    if (i < n) {
        pqc_kernel_gf16_vec lo_t = pqc_kernel_gf16_vec_zero(), hi_t = pqc_kernel_gf16_vec_zero();
        tower_mulacc_block(tower_phi_block(tower_tail_load(a + i, n - i)),
                           tower_phi_block(tower_tail_load(b + i, n - i)),
                           &lo_t, &hi_t);
        const pqc_kernel_gf16_vec t = tower_psi_block(pqc_kernel_gf16_vec_or(lo_t, pqc_kernel_gf16_vec_slli4(hi_t)));
        tower_tail_store(acc + i,
                         pqc_kernel_gf16_vec_xor(tower_tail_load(acc + i, n - i), t), n - i);
    }
}

uint8_t pqc_kernel_gf256_tower_dot_11b_array(const uint8_t *a, const uint8_t *b, size_t n) {
    pqc_kernel_gf256_tower_init();
    pqc_kernel_gf16_vec lo_acc = pqc_kernel_gf16_vec_zero(), hi_acc = pqc_kernel_gf16_vec_zero();
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES)
        tower_mulacc_block(tower_phi_block(pqc_kernel_gf16_vec_load(a + i)),
                           tower_phi_block(pqc_kernel_gf16_vec_load(b + i)), &lo_acc, &hi_acc);
    if (i < n)
        tower_mulacc_block(tower_phi_block(tower_tail_load(a + i, n - i)),
                           tower_phi_block(tower_tail_load(b + i, n - i)),
                           &lo_acc, &hi_acc);
    /* psi 是 GF(2)-线性的，因此“先合并再归约”与“先归约再 psi”等价。 */
    return tower_reduce_xor(tower_psi_block(
        pqc_kernel_gf16_vec_or(lo_acc, pqc_kernel_gf16_vec_slli4(hi_acc))));
}

void pqc_kernel_gf256_tower_const_mul_11b_array(uint8_t *out, const uint8_t *a, uint8_t k,
                                     size_t n) {
    pqc_kernel_gf256_tower_init();
    const uint8_t kt = pqc_kernel_gf256_tower_phi(k);
    const uint8_t k0 = (uint8_t)(kt & 0x0fu), k1 = (uint8_t)(kt >> 4);
    const uint8_t ck1 = pqc_kernel_gf16_cmul_c_local(k1);
    size_t i = 0;
    for (; i + GF16_VEC_BYTES <= n; i += GF16_VEC_BYTES) {
        const pqc_kernel_gf16_vec t = tower_const_mul_block(tower_phi_block(pqc_kernel_gf16_vec_load(a + i)),
                                                 k0, k1, ck1);
        pqc_kernel_gf16_vec_store(out + i, tower_psi_block(t));
    }
    if (i < n)
        tower_tail_store(out + i,
                         tower_psi_block(tower_const_mul_block(
                             tower_phi_block(tower_tail_load(a + i, n - i)),
                             k0, k1, ck1)),
                         n - i);
}

/* ------------------------------------------------------ 矩阵级（塔坐标内核） */

/* M1 常数系数收缩：Y[r][lane] = XOR_c P[r][c] * S[c][lane]（塔坐标）。
 * 系数是常数：每个 (r,c) 的 4 张 shuffle 表在 lane 循环外取一次。 */
void pqc_kernel_gf256_ps_tower(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                    size_t rows, size_t cols, size_t lanes) {
    pqc_kernel_gf256_tower_init();
    for (size_t r = 0; r < rows; ++r) {
        for (size_t lane = 0; lane < lanes; ++lane) out[r * lanes + lane] = 0;
    }
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        for (size_t c = 0; c < cols; ++c) {
            const uint8_t *s = sig + c * lanes;
            const uint8_t k = coeff[r * cols + c];
            const uint8_t k0 = (uint8_t)(k & 0x0fu), k1 = (uint8_t)(k >> 4);
            const uint8_t ck1 = pqc_kernel_gf16_cmul_c_local(k1);
            const pqc_kernel_gf16_vec t0 = vec_mk[k0 & 0x0f];
            const pqc_kernel_gf16_vec t1 = vec_mk[ck1 & 0x0f];
            const pqc_kernel_gf16_vec t2 = vec_mk[k1 & 0x0f];
            const pqc_kernel_gf16_vec t3 = vec_mk[(k0 ^ k1) & 0x0f];
            size_t i = 0;
            for (; i + GF16_VEC_BYTES <= lanes; i += GF16_VEC_BYTES) {
                const pqc_kernel_gf16_vec A = pqc_kernel_gf16_vec_load(s + i);
                const pqc_kernel_gf16_vec a0 = pqc_kernel_gf16_vec_mask_u8(A, 0x0f);
                const pqc_kernel_gf16_vec a1 = pqc_kernel_gf16_vec_srli4(A);
                const pqc_kernel_gf16_vec low = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_shuffle(t0, a0),
                                                  pqc_kernel_gf16_vec_shuffle(t1, a1));
                const pqc_kernel_gf16_vec high = pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_shuffle(t2, a0),
                                                   pqc_kernel_gf16_vec_shuffle(t3, a1));
                const pqc_kernel_gf16_vec t = pqc_kernel_gf16_vec_or(low, pqc_kernel_gf16_vec_slli4(high));
                pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_xor(pqc_kernel_gf16_vec_load(o + i), t));
            }
            for (; i < lanes; ++i) o[i] ^= tower_const_mul_byte(s[i], k0, k1, ck1);
        }
    }
}

/* M1（0x11B 入/出）：系数与信号各转换一次，输出统一转回，转换被所有行列复用。
 * scratch 需要 max(rows*cols, cols*lanes, rows*lanes) 字节。 */
void pqc_kernel_gf256_ps_tower_11b(uint8_t *out, const uint8_t *coeff, const uint8_t *sig,
                        size_t rows, size_t cols, size_t lanes, uint8_t *scratch) {
    pqc_kernel_gf256_tower_init();
    uint8_t *coeff_t = scratch;
    uint8_t *sig_t = scratch + rows * cols;
    uint8_t *out_t = scratch + rows * cols + cols * lanes;
    for (size_t i = 0; i < rows * cols; ++i) coeff_t[i] = pqc_kernel_gf256_tower_phi(coeff[i]);
    pqc_kernel_gf256_tower_phi_array(sig_t, sig, cols * lanes);
    pqc_kernel_gf256_ps_tower(out_t, coeff_t, sig_t, rows, cols, lanes);
    pqc_kernel_gf256_tower_psi_array(out, out_t, rows * lanes);
}

static void tower_sts_scalar(uint8_t *out, const uint8_t *left, const uint8_t *right,
                             size_t rows, size_t cols, size_t lanes) {
    for (size_t r = 0; r < rows; ++r) {
        for (size_t lane = 0; lane < lanes; ++lane) {
            uint8_t acc = 0;
            for (size_t c = 0; c < cols; ++c) {
                acc ^= pqc_kernel_gf256_tower_mul(left[(r * cols + c) * lanes + lane],
                                       right[c * lanes + lane]);
            }
            out[r * lanes + lane] = acc;
        }
    }
}

/* M2 变量矩阵收缩：Y[r][lane] = XOR_c S[(r*cols+c)][lane] * T[c][lane]（塔坐标）。
 * low / high 两个坐标分别累加，循环里不做重新打包。 */
void pqc_kernel_gf256_sts_tower(uint8_t *out, const uint8_t *left, const uint8_t *right,
                     size_t rows, size_t cols, size_t lanes) {
    pqc_kernel_gf256_tower_init();
    if (lanes % GF16_VEC_BYTES) {
        tower_sts_scalar(out, left, right, rows, cols, lanes);
        return;
    }
    for (size_t r = 0; r < rows; ++r) {
        uint8_t *o = out + r * lanes;
        /* lane 之间互相独立，因此每个 lane 向量各有一对累加器。 */
        for (size_t i = 0; i < lanes; i += GF16_VEC_BYTES) {
            pqc_kernel_gf16_vec lo_acc = pqc_kernel_gf16_vec_zero(), hi_acc = pqc_kernel_gf16_vec_zero();
            for (size_t c = 0; c < cols; ++c) {
                const uint8_t *s = left + (r * cols + c) * lanes + i;
                const uint8_t *t = right + c * lanes + i;
                tower_mulacc_block(pqc_kernel_gf16_vec_load(s), pqc_kernel_gf16_vec_load(t), &lo_acc, &hi_acc);
            }
            pqc_kernel_gf16_vec_store(o + i, pqc_kernel_gf16_vec_or(lo_acc, pqc_kernel_gf16_vec_slli4(hi_acc)));
        }
    }
}

void pqc_kernel_gf256_sts_tower_11b(uint8_t *out, const uint8_t *left, const uint8_t *right,
                         size_t rows, size_t cols, size_t lanes, uint8_t *scratch) {
    pqc_kernel_gf256_tower_init();
    const size_t left_n = rows * cols * lanes;
    uint8_t *left_t = scratch;
    uint8_t *right_t = scratch + left_n;
    uint8_t *out_t = scratch + left_n + cols * lanes;
    pqc_kernel_gf256_tower_phi_array(left_t, left, left_n);
    pqc_kernel_gf256_tower_phi_array(right_t, right, cols * lanes);
    pqc_kernel_gf256_sts_tower(out_t, left_t, right_t, rows, cols, lanes);
    pqc_kernel_gf256_tower_psi_array(out, out_t, rows * lanes);
}

/* --------------------------------------------------------------- 自检 */

/* 穷举校验：phi/psi 互逆、塔乘法等价于 0x11B 乘法、平方/求逆/Frobenius 自洽。
 * 返回 1 表示全部通过。 */
int pqc_kernel_gf256_tower_selftest(void) {
    pqc_kernel_gf256_tower_init();
    if (sizeof(TOWER_INV16) != 16) return 0;
    for (unsigned a = 0; a < 256; ++a) {
        const uint8_t pa = pqc_kernel_gf256_tower_phi((uint8_t)a);
        if (pqc_kernel_gf256_tower_psi(pa) != (uint8_t)a) return 0;
        if (pqc_kernel_gf256_tower_square(pa) != pqc_kernel_gf256_tower_mul(pa, pa)) return 0;
        if (pqc_kernel_gf256_tower_frobenius(pa) != pqc_kernel_gf256_tower_pow(pa, 16)) return 0;
        if (a == 0) {
            if (pqc_kernel_gf256_tower_inv(pa) != 0) return 0;
        } else if (pqc_kernel_gf256_tower_mul(pqc_kernel_gf256_tower_inv(pa), pa) != 0x01) {
            return 0;
        }
        for (unsigned b = 0; b < 256; ++b) {
            const uint8_t pb = pqc_kernel_gf256_tower_phi((uint8_t)b);
            if (pqc_kernel_gf256_tower_psi(pqc_kernel_gf256_tower_mul(pa, pb)) !=
                pqc_kernel_gf256_mul_ref((uint8_t)a, (uint8_t)b)) {
                return 0;
            }
        }
    }
    return 1;
}

/* ---------------------------------------------------- 向量级公开内核（供嵌入）
 * 输入输出都是塔坐标，边界转换由调用方负责。第三方代码（例如 MQOM 的字段层）
 * 可以把这几条内核直接编进内循环，而不必走逐数组的 0x11B 边界接口。
 */
pqc_kernel_gf16_vec pqc_kernel_gf256_tower_phi_block(pqc_kernel_gf16_vec v) {
    pqc_kernel_gf256_tower_init();
    return tower_phi_block(v);
}

pqc_kernel_gf16_vec pqc_kernel_gf256_tower_psi_block(pqc_kernel_gf16_vec v) {
    pqc_kernel_gf256_tower_init();
    return tower_psi_block(v);
}

pqc_kernel_gf16_vec pqc_kernel_gf256_tower_mul_block(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b) {
    pqc_kernel_gf256_tower_init();
    return tower_mul_block(a, b);
}

void pqc_kernel_gf256_tower_mulacc_block(pqc_kernel_gf16_vec a, pqc_kernel_gf16_vec b, pqc_kernel_gf16_vec *lo, pqc_kernel_gf16_vec *hi) {
    pqc_kernel_gf256_tower_init();
    tower_mulacc_block(a, b, lo, hi);
}

pqc_kernel_gf16_vec pqc_kernel_gf256_tower_const_mul_block(pqc_kernel_gf16_vec a, uint8_t k) {
    pqc_kernel_gf256_tower_init();
    const uint8_t kt = pqc_kernel_gf256_tower_phi(k);
    const uint8_t k0 = (uint8_t)(kt & 0x0fu);
    const uint8_t k1 = (uint8_t)(kt >> 4);
    return tower_const_mul_block(a, k0, k1, pqc_kernel_gf16_cmul_c_local(k1));
}
