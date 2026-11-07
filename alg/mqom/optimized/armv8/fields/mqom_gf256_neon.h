/*
 * mqom_gf256_neon.h —— MOQM ARMv8 常数路径 NEON 内核
 *
 * 目标：替换 MOQM 现有 fields_neon.h 里"用通用 8 轮位串行去乘一个常数"的地方。
 *
 * 关键恒等式（GF(256) = Rijndael 0x11B，普通多项式基）：
 *   乘法对 GF(2) 是线性的，且字节 a = (a & 0x0f) XOR ((a & 0xf0))
 *   把 a 拆成低 nibble 与高 nibble 两部分，则
 *       k*a = k*a_lo  XOR  k*(a_hi * x^4)
 *   其中 a_lo 只取 16 种值，a_hi 也只取 16 种值，于是
 *       k*a = T0[a_lo] XOR T1[a_hi]
 *       T0[v] = k*v            （v 是低 nibble 对应的域元素）
 *       T1[v] = (k*x^4)*v      （v 是高 nibble，乘完整体左移 4 位）
 *
 *   ⇒ 每个 16 字节向量只要 2 次 vqtbl1q_u8 + 1 次异或。
 *   vqtbl1q_u8 是"寄存器内查表"，不碰内存，索引不参与任何分支/地址计算，
 *   因此仍然满足常数时间要求（这也是 GFNI 之所以能用在密码实现里的同一道理）。
 *
 * 对比（GF_cal 在 ARMv8 上用 perf_event CPU cycles 实测，单位 cycles/逻辑乘）：
 *   - 通用变量乘（MQOM 现有 mqom_neon_mul，8 轮位串行）      2.478
 *   - 塔分解常数乘（8 次 shuffle，含 phi/psi）                 0.815
 *   - 本文件的 nibble 拆分常数乘（2 次 shuffle）              ≈ 0.19  ← 推导值
 *   实测 vqtbl1q_u8 吞吐 ≈ 1.5 周期/条（GF_cal 两条独立路径互证）
 *
 * 本文件不依赖 MOQM 任何头文件，只依赖 <arm_neon.h>，便于单独交叉编译验证。
 */

#ifndef MQOM_GF256_NEON_H
#define MQOM_GF256_NEON_H

/* vqtbl1q_u8 是 AArch64 专属 intrinsic（ARMv7 的 vtbl 只支持 8 字节表）。
 * 显式挡在这里，免得跑到 ARM 机器上才发现 "implicit declaration of vqtbl1q_u8"。 */
#if !defined(__aarch64__)
#error "mqom_gf256_neon.h requires AArch64: vqtbl1q_u8 is an AArch64-only intrinsic."
#endif

#include <arm_neon.h>
#include <stdint.h>
#include <string.h>

/* ------------------------------------------------------------------ 标量参考
 * 只用于"建表"和自检对拍，不在热路径上。GF(2^8)/0x11B。
 */
static inline uint8_t mqom_gf256_mul_scalar(uint8_t a, uint8_t b)
{
	uint8_t r = 0;
	for (int i = 0; i < 8; ++i) {
		r = (uint8_t)(r ^ (uint8_t)(-(int)(b & 1u) & (int)a));
		const uint8_t hi = (uint8_t)(a & 0x80u);
		a = (uint8_t)(a << 1);
		if (hi != 0u) {
			a = (uint8_t)(a ^ 0x1bu);
		}
		b = (uint8_t)(b >> 1);
	}
	return r;
}

/* GF(16) -> GF(256) 嵌入（与 MOQM fields_neon.h 的 gf16_gf256_mult_neon(x,1) 等价）。 */
static inline uint8_t mqom_gf16_embed_gf256(uint8_t x)
{
	return (uint8_t)(((x + 0u) & 0x01u)
		^ (uint8_t)(-((int)((x >> 1) & 1u)) & 0xe0)
		^ (uint8_t)(-((int)((x >> 2) & 1u)) & 0x5d)
		^ (uint8_t)(-((int)((x >> 3) & 1u)) & 0xb0));
}

/* GF(16) -> GF(256) 嵌入表（固定常量）。
 * 由 embed(x) = (x&1)·0x01 ^ ((x>>1)&1)·0xe0 ^ ((x>>2)&1)·0x5d ^ ((x>>3)&1)·0xb0 求出，
 * selftest 用例 [1] 逐项对拍 MOQM 的原始表达式。 */
static inline const uint8_t *mqom_gf16_embed_table(void)
{
	static const uint8_t EMB16[16] = {
		0x00, 0x01, 0xe0, 0xe1, 0x5d, 0x5c, 0xbd, 0xbc,
		0xb0, 0xb1, 0x50, 0x51, 0xed, 0xec, 0x0d, 0x0c
	};
	return EMB16;
}

/* -------------------------------------------------------- 常数乘：nibble 拆分
 *
 * 建表成本是关键：MQOM 的实际数组长度只有 32 / 58 / 116 字节（ETA / MQ_M/MU / MQ_N），
 * 而 API 是“一次调用一棵表、零复用”。用 16 次标量乘建表（每次 8 轮、含数据相关
 * 分支，约 320 周期）比直接用 NEON 变量乘还慢，所以必须把建表压到几十条指令。
 */

/* t[v] = k*v （v = 0..15），幂基 + XOR 树。
 *   k·v = ⊕_{i: v_i=1} (k·x^i),
 *   p_{i+1} = (p_i << 1) ^ (-(p_i>>7) & 0x1b)
 * 3 次移位归约 + 16 次异或，全 branchless（不用 ?:，避免密钥相关分支）。
 * v 是公开循环变量，编译器会完全展开并常量折叠。 */
static inline void mqom_gf256_mul_table16(uint8_t k, uint8_t t[16])
{
	const uint8_t p0 = k;
	const uint8_t p1 = (uint8_t)((uint8_t)(k << 1) ^ (uint8_t)(-(int)((k >> 7) & 1u) & 0x1bu));
	const uint8_t p2 = (uint8_t)((uint8_t)(p1 << 1) ^ (uint8_t)(-(int)((p1 >> 7) & 1u) & 0x1bu));
	const uint8_t p3 = (uint8_t)((uint8_t)(p2 << 1) ^ (uint8_t)(-(int)((p2 >> 7) & 1u) & 0x1bu));
	for (unsigned v = 0; v < 16u; ++v) {
		t[v] = (uint8_t)(((v & 8u) ? p3 : 0u) ^ ((v & 4u) ? p2 : 0u)
		               ^ ((v & 2u) ? p1 : 0u) ^ ((v & 1u) ? p0 : 0u));
	}
}

typedef struct {
	uint8x16_t t0;   /* T0[v] = k*v          */
	uint8x16_t t1;   /* T1[v] = (k*x^4)*v    */
} mqom_gf256_cmul_tab;

static inline mqom_gf256_cmul_tab mqom_gf256_cmul_tab_build(uint8_t k)
{
	mqom_gf256_cmul_tab tab;
	uint8_t t0[16];
	uint8_t t1[16];
	/* k*x^4 由同一条移位链再走一步得到 */
	const uint8_t p1 = (uint8_t)((uint8_t)(k << 1) ^ (uint8_t)(-(int)((k >> 7) & 1u) & 0x1bu));
	const uint8_t p2 = (uint8_t)((uint8_t)(p1 << 1) ^ (uint8_t)(-(int)((p1 >> 7) & 1u) & 0x1bu));
	const uint8_t p3 = (uint8_t)((uint8_t)(p2 << 1) ^ (uint8_t)(-(int)((p2 >> 7) & 1u) & 0x1bu));
	const uint8_t k4 = (uint8_t)((uint8_t)(p3 << 1) ^ (uint8_t)(-(int)((p3 >> 7) & 1u) & 0x1bu));
	mqom_gf256_mul_table16(k, t0);
	mqom_gf256_mul_table16(k4, t1);
	tab.t0 = vld1q_u8(t0);
	tab.t1 = vld1q_u8(t1);
	return tab;
}

/* 一个 16 字节向量乘以常数 k：2 次 shuffle + 1 次异或。 */
static inline uint8x16_t mqom_gf256_cmul_vec(uint8x16_t v, const mqom_gf256_cmul_tab *tab)
{
	const uint8x16_t idx_lo = vandq_u8(v, vdupq_n_u8(0x0fu));
	const uint8x16_t idx_hi = vshrq_n_u8(v, 4);
	return veorq_u8(vqtbl1q_u8(tab->t0, idx_lo), vqtbl1q_u8(tab->t1, idx_hi));
}

/* c[i] = k * a[i]。替换 fields_neon.h 的 gf256_constant_vect_mult_neon。 */
static inline void mqom_gf256_const_vect_mul_neon(uint8_t k, const uint8_t *a,
		uint8_t *c, uint32_t len)
{
	const mqom_gf256_cmul_tab tab = mqom_gf256_cmul_tab_build(k);
	uint32_t i = 0;
	for (; i + 16u <= len; i += 16u) {
		vst1q_u8(c + i, mqom_gf256_cmul_vec(vld1q_u8(a + i), &tab));
	}
	if (i < len) {
		uint8_t tmp[16];
		memset(tmp, 0, sizeof(tmp));
		memcpy(tmp, a + i, (size_t)(len - i));
		vst1q_u8(tmp, mqom_gf256_cmul_vec(vld1q_u8(tmp), &tab));
		memcpy(c + i, tmp, (size_t)(len - i));
	}
}

/* ------------------------------------------------- GF(16) packed -> GF(256) 提升 */

/* 8 个 packed 字节（低 nibble 在前）→ 16 个 GF(16) 元素，每元素一字节。 */
static inline uint8x16_t mqom_gf16_unpack16_low_first(const uint8_t *p)
{
	const uint8x8_t v  = vld1_u8(p);
	const uint8x8_t lo = vand_u8(v, vdup_n_u8(0x0fu));
	const uint8x8_t hi = vshr_n_u8(v, 4);
	return vcombine_u8(vzip1_u8(lo, hi), vzip2_u8(lo, hi));
}

/* 同前，但只保留前 count 个 lane（count<=16），其余 lane 归零。
 * 奇数 count 时最后一个字节的高 nibble 属于填充，必须屏蔽掉。 */
static inline uint8x16_t mqom_gf16_unpack_n_low_first(const uint8_t *p, uint32_t count)
{
	static const uint8_t LANE_IDX[16] = {
		0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u
	};
	const uint8x16_t nib = mqom_gf16_unpack16_low_first(p);
	/* 被屏蔽的 lane 索引变 0，而 embed(0) == 0，所以查表后自然是 0。 */
	return vandq_u8(nib, vcltq_u8(vld1q_u8(LANE_IDX), vdupq_n_u8((uint8_t)count)));
}

/* 16 项嵌入表向量（固定常量，不再每次调用重算）。 */
static inline uint8x16_t mqom_gf16_embed_table_vec(void)
{
	return vld1q_u8(mqom_gf16_embed_table());
}

/* out[i] = embed(nibble_i(packed))，n 是 GF(16) 元素个数（不是字节数）。
 * 替换 fields_neon.h 的 gf256_vect_lift_from_gf16_neon。 */
static inline void mqom_gf16_lift_vect_neon(const uint8_t *packed, uint8_t *out, uint32_t n)
{
	const uint8x16_t tab = mqom_gf16_embed_table_vec();
	uint32_t i = 0;
	for (; i + 16u <= n; i += 16u) {
		vst1q_u8(out + i, vqtbl1q_u8(tab, mqom_gf16_unpack16_low_first(packed + (i >> 1))));
	}
	if (i < n) {
		const uint32_t rem = n - i;
		uint8_t src[8];
		uint8_t tmp[16];
		memset(src, 0, sizeof(src));
		memcpy(src, packed + (i >> 1), (size_t)((rem + 1u) / 2u));
		vst1q_u8(tmp, vqtbl1q_u8(tab, mqom_gf16_unpack16_low_first(src)));
		memcpy(out + i, tmp, (size_t)rem);
	}
}

/* 与 MOQM mqom_gf16_lift_block 同签名：返回一个向量，前 count 个 lane 是提升结果，
 * 其余 lane 为 0（原实现用 `uint8_t lifted[16] = {0}` 保证这一点）。
 * 原实现每次调用都要跑 16 次标量乘重建表 + memcpy + 逐字节查表，这里是 1 次 shuffle。 */
static inline uint8x16_t mqom_gf16_lift_block_neon(const uint8_t *packed, uint32_t count)
{
	const uint8x16_t tab = mqom_gf16_embed_table_vec();
	uint8_t src[8];
	uint8_t tmp[16];
	if (count > 16u) {
		count = 16u;
	}
	memset(src, 0, sizeof(src));
	memcpy(src, packed, (size_t)((count + 1u) / 2u));
	vst1q_u8(tmp, vqtbl1q_u8(tab, mqom_gf16_unpack_n_low_first(src, count)));
	return vld1q_u8(tmp);
}

/* -------------------------------------- 常数乘 + 提升融合（M1 形状的主力内核）
 *
 * c[i] = k * embed(nibble_i(b))，n 是 GF(16) 元素个数。
 * 常数的乘法直接烘进表里 ⇒ 仍然只要 1 次 shuffle，没有任何逐元素乘法。
 * 替换 fields_neon.h 的 gf256_gf16_constant_vect_mult_neon。
 */
static inline void mqom_gf256_gf16_const_vect_mul_neon(uint8_t k, const uint8_t *b,
		uint8_t *c, uint32_t n)
{
	/* tab[x] = k * embed(x)。embed 是固定表，而 k*· 对 GF(2) 线性，
	 * 所以整张表可以直接从 k 的 nibble 表查出来 —— 不需要 16 次标量乘。 */
	uint8_t tab[16];
	uint8_t t0a[16];
	uint8_t t1a[16];
	uint32_t i = 0;
	const uint8_t *emb = mqom_gf16_embed_table();
	const mqom_gf256_cmul_tab ct = mqom_gf256_cmul_tab_build(k);
	vst1q_u8(t0a, ct.t0);
	vst1q_u8(t1a, ct.t1);
	for (unsigned x = 0; x < 16u; ++x) {
		const uint8_t e = emb[x];
		tab[x] = (uint8_t)(t0a[e & 0x0fu] ^ t1a[e >> 4]);
	}
	{
		const uint8x16_t tv = vld1q_u8(tab);
		for (; i + 16u <= n; i += 16u) {
			vst1q_u8(c + i, vqtbl1q_u8(tv, mqom_gf16_unpack16_low_first(b + (i >> 1))));
		}
		if (i < n) {
			const uint32_t rem = n - i;
			uint8_t src[8];
			uint8_t tmp[16];
			memset(src, 0, sizeof(src));
			memcpy(src, b + (i >> 1), (size_t)((rem + 1u) / 2u));
			vst1q_u8(tmp, vqtbl1q_u8(tv, mqom_gf16_unpack16_low_first(src)));
			memcpy(c + i, tmp, (size_t)rem);
		}
	}
}

/* ------------------------------------- GF(16) 常数 × GF(256) 向量（c[i]=embed(a)*b[i]）
 * 这就是一个普通的 GF(256) 常数乘，直接走 nibble 拆分。
 * 替换 fields_neon.h 的 gf16_gf256_constant_vect_mult_neon。
 */
static inline void mqom_gf16_gf256_const_vect_mul_neon(uint8_t a_gf16, const uint8_t *b,
		uint8_t *c, uint32_t n)
{
	mqom_gf256_const_vect_mul_neon(mqom_gf16_embed_gf256(a_gf16), b, c, n);
}

/* ------------------------------------------------- 打包 GF(16) 输入的点积（提升+乘+归约）
 * 供 gf16_gf256_vect_mult_neon 使用；这里只提供"提升"部分，乘加复用 MOQM 现成的
 * mqom_neon_mul，避免重复实现变量乘内核。
 *
 * 说明：调用方按 16 元素一块调用 mqom_gf16_lift_block_neon 即可，本函数不重复给出。
 */

#endif /* MQOM_GF256_NEON_H */
