/**
 * @file hash_aes_shake.c
 * @brief SHAKE 家族的增量适配与 Scloud+ 的 F/G/H/K 原语。
 *
 * Keccak-f[1600] 置换、海绵吸收/挤压和一次性 SHAKE 都由库内 src/sym/keccak.c
 * 提供，参考包自带的 TweetFips202 副本已删除。本文件只保留提交代码实际使用
 * 的块式接口（sample.c 依赖 squeezeblocks）以及按标签域分隔的 F/G/H/K。
 */

#include "hash.h"
#include "scloudplus_param_common.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

void shake128_absorb(keccak_state *state, const uint8_t *in, size_t inlen)
{
	(void)sha3_update(state, in, inlen);
}

void shake128_finalize(keccak_state *state)
{
	shake_xof_final(state);
}

void shake128_squeeze(uint8_t *out, size_t outlen, keccak_state *state)
{
	shake_squeeze(state, out, outlen);
}

void shake128_absorb_once(keccak_state *state, const uint8_t *in, size_t inlen)
{
	shake128_init(state);
	(void)sha3_update(state, in, inlen);
	shake_xof_final(state);
}

void shake128_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state)
{
	shake_squeeze(state, out, nblocks * SHAKE128_RATE);
}

void shake256_absorb(keccak_state *state, const uint8_t *in, size_t inlen)
{
	(void)sha3_update(state, in, inlen);
}

void shake256_finalize(keccak_state *state)
{
	shake_xof_final(state);
}

void shake256_squeeze(uint8_t *out, size_t outlen, keccak_state *state)
{
	shake_squeeze(state, out, outlen);
}

void shake256_absorb_once(keccak_state *state, const uint8_t *in, size_t inlen)
{
	shake256_init(state);
	(void)sha3_update(state, in, inlen);
	shake_xof_final(state);
}

void shake256_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state)
{
	shake_squeeze(state, out, nblocks * SHAKE256_RATE);
}

enum
{
	SCLOUDPLUS_HASH_LABEL_F = 0x46U,
	SCLOUDPLUS_HASH_LABEL_G = 0x47U,
	SCLOUDPLUS_HASH_LABEL_H = 0x48U,
	SCLOUDPLUS_HASH_LABEL_K = 0x4BU
};

static void scloudplus_xof_by_level(uint8_t *output, size_t outlen,
									uint8_t label, const uint8_t *input,
									size_t inlen)
{
	keccak_state state;

	shake256_init(&state);
	(void)sha3_update(&state, &label, 1U);
	(void)sha3_update(&state, input, inlen);
	shake_xof_final(&state);
	shake_squeeze(&state, output, outlen);
	memset(&state, 0, sizeof(state));
}

#if scloudplus_hash_bytes == 64
static void sha3_512_labeled(uint8_t h[64], uint8_t label,
							 const uint8_t *in, size_t inlen)
{
	KECCAK1600_CTX ctx;

	sha3_init(&ctx, 0x06, 512);
	(void)sha3_update(&ctx, &label, 1U);
	(void)sha3_update(&ctx, in, inlen);
	sha3_final(h, &ctx);
	memset(&ctx, 0, sizeof(ctx));
}
#endif

void scloudplus_F(unsigned char *output, unsigned long long outlen,
				  const unsigned char *input, unsigned long long inlen)
{
	scloudplus_xof_by_level(output, (size_t)outlen, SCLOUDPLUS_HASH_LABEL_F,
							input, (size_t)inlen);
}
void scloudplus_K(unsigned char *output, unsigned long long outlen,
				  const unsigned char *input, unsigned long long inlen)
{
	scloudplus_xof_by_level(output, (size_t)outlen, SCLOUDPLUS_HASH_LABEL_K,
							input, (size_t)inlen);
}
void scloudplus_H(unsigned char *output, const unsigned char *input,
				  unsigned long long inlen)
{
#if scloudplus_hash_bytes == 64
	sha3_512_labeled(output, SCLOUDPLUS_HASH_LABEL_H, input, (size_t)inlen);
#else
	/*
	 * SHA3 has fixed 224/256/384/512-bit variants.  For the 384- and
	 * 512-bit Scloud+ parameter sets the agreed instantiation requires
	 * 768- or 1024-bit H output, so AES/SHAKE builds use the same labeled
	 * SHAKE-256 expansion interface as F/G/K for those longer hash strings.
	 */
	scloudplus_xof_by_level(output, scloudplus_hash_bytes,
							SCLOUDPLUS_HASH_LABEL_H, input, (size_t)inlen);
#endif
}
void scloudplus_G(unsigned char *output, const unsigned char *input,
				  unsigned long long inlen)
{
	scloudplus_xof_by_level(output, (size_t)scloudplus_G_bytes,
							SCLOUDPLUS_HASH_LABEL_G, input, (size_t)inlen);
}
