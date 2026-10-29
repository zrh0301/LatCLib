#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include <stddef.h>
#include <stdint.h>
#include "keccak.h"

typedef KECCAK1600_CTX shake_t;

void snova_shake128_init(shake_t *instance);
void snova_shake256_init(shake_t *instance);
void snova_shake_absorb(shake_t *instance, const uint8_t *in, size_t inlen);
void snova_shake_finalize(shake_t *instance);
void snova_shake_squeeze(shake_t *instance, uint8_t *out, size_t outlen);

void snova_pk_expand(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

// Incremental API
typedef struct {
	uint64_t state[56];
} snova_pk_expander_t;

void snova_pk_expander_init(snova_pk_expander_t *instance, const uint8_t *seed, size_t input_bytes);
void snova_pk_expander_squeeze(uint8_t *data, size_t num_bytes, snova_pk_expander_t *instance);
void snova_pk_expander_goto(snova_pk_expander_t *instance, size_t index);
void snova_pk_expander_free(snova_pk_expander_t *instance);

void snova_secure_clear(void *ptr, size_t len);
// #define snova_secure_clear(ptr, len) do { } while (0)

#define SNOVA_SECURE_CLEAR(x) snova_secure_clear(x, sizeof(x));
#define SNOVA_SECURE_CLEARI(x) snova_secure_clear(&x, sizeof(x));

#endif
