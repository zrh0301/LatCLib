/**
 * @file util.c
 * @brief Scloud+ 的常量时间比较/选择与内存擦除。
 *
 * 这些算子与库内 util/ct.h、util/secure_clear.h 完全同义，这里只保留参考实现
 * 的函数名与返回值约定（比较失败返回 0xFF），实现直接复用共享算子。
 */

#include "util.h"

#include "util/ct.h"
#include "util/secure_clear.h"

#include <string.h>

/**
 * @brief Compare two byte strings without data-dependent early exit.
 *
 * Returns 0xFF when the strings differ (the "failure" mask used by the FO
 * transform) and 0x00 when they match.
 */
uint8_t scloudplus_verify(const uint8_t *a, const uint8_t *b, size_t len)
{
	return (uint8_t)(0U - (uint32_t)(pqc_ct_equal(a, b, len) ^ 1));
}

/**
 * @brief Select between byte strings with a byte mask and no secret branch.
 *
 * ``mask`` is 0x00 (keep ``a``) or any nonzero byte such as 0xFF (take ``b``).
 */
void scloudplus_cmov(uint8_t *r, const uint8_t *a, const uint8_t *b, size_t len,
					 uint8_t mask)
{
	/* Map any nonzero mask byte to 1 without a data-dependent branch. */
	const uint8_t select = (uint8_t)(0U - (uint32_t)((mask | (uint8_t)(0U - mask)) >> 7));

	if (r == b)
	{
		return;
	}
	memmove(r, a, len);
	pqc_ct_select(r, b, len, select);
}

/**
 * @brief Wipe memory without relying on libc-specific explicit_bzero APIs.
 */
void scloudplus_secure_zeroize(void *ptr, size_t len)
{
	pqc_secure_clear(ptr, len);
}
