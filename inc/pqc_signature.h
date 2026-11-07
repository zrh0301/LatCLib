#ifndef PQC_SIGNATURE_H
#define PQC_SIGNATURE_H
#include <stddef.h>
#include <stdint.h>
/* Signature algorithms share this signed-message interface. Each format stays
 * algorithm-specific. Buffers must be disjoint and sized using these lengths.
 * sign/open return 0 on success and nonzero on failure. open's output capacity
 * must be at least signed_length. Output is not authenticated on failure.
 * These descriptors use the system RNG. */
typedef struct pqc_signature {
 const char *name;
 size_t public_key_bytes, secret_key_bytes, signature_bytes;
 int (*keypair)(uint8_t *,uint8_t *);
 int (*sign)(uint8_t *,size_t *,const uint8_t *,size_t,const uint8_t *);
 int (*open)(uint8_t *,size_t *,const uint8_t *,size_t,const uint8_t *);
} pqc_signature;
const pqc_signature *pqc_snova_r2_algorithm(void);
const pqc_signature *pqc_snova_r3_algorithm(void);
const pqc_signature *pqc_mayo_algorithm(void);
const pqc_signature *pqc_mqom_algorithm(void);
#endif
