#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include <stddef.h>
#include <stdint.h>
#include "keccak.h"
#include "aes.h"
#include "secure_clear.h"

typedef KECCAK1600_CTX shake_t;

/* 哈希与 AES 均直接使用 pqc 公共层 (inc/keccak.h, inc/aes.h), 不再单独包装。 */

void snova_pk_expand(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

// Incremental API
typedef struct {
	uint64_t state[56];
} snova_pk_expander_t;

void snova_pk_expander_init(snova_pk_expander_t *instance, const uint8_t *seed, size_t input_bytes);
void snova_pk_expander_squeeze(uint16_t *data, size_t num_bytes, snova_pk_expander_t *instance);
void snova_pk_expander_goto(snova_pk_expander_t *instance, size_t index);
void snova_pk_expander_free(snova_pk_expander_t *instance);

/* 敏感数据清零统一使用 pqc_secure_clear (src/secure_clear.c)。
 * 保留别名与外部覆盖点, 兼容原有调用方式和插桩测试。 */
#ifndef snova_secure_clear
#define snova_secure_clear(ptr, len) pqc_secure_clear(ptr, len)
#endif

#define SNOVA_SECURE_CLEAR(x) snova_secure_clear(x, sizeof(x));
#define SNOVA_SECURE_CLEARI(x) snova_secure_clear(&x, sizeof(x));

#endif
