#ifndef PQC_KEM_H
#define PQC_KEM_H
#include <stddef.h>
#include <stdint.h>
/* KEM algorithms share this key-encapsulation interface. The ciphertext and
 * shared secret are separate buffers, unlike the signed-message API. All
 * buffers must be disjoint and sized using these lengths. Every function
 * returns 0 on success and nonzero on failure; on failure the outputs are not
 * meaningful. These descriptors use the system RNG. */
typedef struct pqc_kem {
 const char *name;
 size_t public_key_bytes, secret_key_bytes, ciphertext_bytes, shared_secret_bytes;
 int (*keypair)(uint8_t *,uint8_t *);
 int (*encaps)(uint8_t *,uint8_t *,const uint8_t *);
 int (*decaps)(uint8_t *,const uint8_t *,const uint8_t *);
} pqc_kem;
const pqc_kem *pqc_ngcc_aigis_algorithm(void);
const pqc_kem *pqc_ngcc_scloud_algorithm(void);
#endif
