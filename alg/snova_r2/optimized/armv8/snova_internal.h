#ifndef SNOVA_INTERNAL_H
#define SNOVA_INTERNAL_H

#include "snova_profile.h"
#include "snova_public_layout.h"
#include "symmetric.h"

#define Keccak_HashInstance shake_t
#define Keccak_HashInitialize_SHAKE256 shake256_init
#define Keccak_HashUpdate(a, b, c) shake_absorb(a, b, (c) / 8)
#define Keccak_HashFinal(a, b) shake_finalize(a)
#define Keccak_HashSqueeze(a, b, c) shake_squeeze(b, (c) / 8, a)

static inline void snova_pk_expand_seed(const uint8_t* seed, uint8_t* pk_bytes) {
	snova_pk_expand(pk_bytes, NUM_GEN_PUB_BYTES, seed, SEED_LENGTH_PUBLIC);
}

static inline void snova_signed_hash(const uint8_t* digest, uint64_t bytes_digest, const uint8_t* pt_public_key_seed,
                                      const uint8_t* array_salt, uint8_t* signed_hash_out) {
	Keccak_HashInstance hashInstance;

	Keccak_HashInitialize_SHAKE256(&hashInstance);
	Keccak_HashUpdate(&hashInstance, pt_public_key_seed, 8 * seed_length_public);
	Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
	Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
	Keccak_HashFinal(&hashInstance, NULL);
	Keccak_HashSqueeze(&hashInstance, signed_hash_out, 8 * bytes_hash);
}

void snova_gen_seeds_and_T12(T12_t T12, const uint8_t* seed);

void snova_gen_ABQ(map_group1* map, const uint8_t* pt_public_key_seed);

void snova_gen_F_vtl(map_group2* map2, map_group1* map1, T12_t T12);

void snova_gen_P22_vtl(P22_byte_t outP22, T12_t T12, P21_t P21, F12_t F12);

#if FIXED_ABQ
void snova_gen_fixed_ABQ(const char* abq_seed);
#endif

int snova_sign_core(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, const uint8_t* array_salt,
                     const sk_gf16* sk);

int snova_verify_core(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature,
                       const public_key_expand* pkx);

#endif
