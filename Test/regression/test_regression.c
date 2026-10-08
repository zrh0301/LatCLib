#include "dili.h"
#include "kyber.h"
#include "keccak.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Behavior regression values captured before the first-round refactor.
 * These are NOT standards known-answer tests or proof of conformance.
 * SHA3-256 covers serialized bytes and explicit little-endian coefficients;
 * unused slots, host endianness and C struct padding are excluded.
 */
static const char *const kyber_expected[3][3] = {
    {"ae6c97390ec2cfb8cf67262045325e6d7f7b5b3ef2ba76ad72074b3d5588a0df",
     "eb38702e8584738e22bc607ede417d9f2f284e3d7bc1835ff234ac1dee2e0956",
     "60af065e5c76dcabab943bec8bf64b942cd1d98071dcdfd635be38c45cd5eb63"},
    {"f47d9f024b754ac9ce0ecf54e1b3ca204d34e5ab8cb3ca9c7dfe8b8e55df9011",
     "8aee4b36be7655ed3a458e54c431b5ee524fcc1ce8be7bb694510c4b8a13c45f",
     "6524a3c7994cc1a6787e1a2b1a92109a19d94e4585cc4e196a8bb3f04f29670c"},
    {"65f05d3cd4752a6ca354a826d78be07f350afcc0ca0ac86292fa354b1ff34303",
     "0d2332e96d1dbf2c319f9328fe4722dde3f70a3b3a3d4f4603e16bff67803433",
     "6a535f9731c37b4aa480332c0b8d554622c236b88d047469758310b76a788401"},
};

static const char *const dilithium_expected[3][3] = {
    {"2f04e40b322b9b9e4d817e150418ba11327ec1856216dcd81588da0342c764d4",
     "fbd19faeaea8706d4b906041a56e9b0f6f6a025e327186b4d2851199beb390ae",
     "5fc8291d1570d265c91d1c28f36f5ed9d2cb32d14bf429c830123e744ba8a902"},
    {"176528ec5915bcb345142e8d7e51e2e51e3b8c1f27c4c2131bf0a80338479bf8",
     "b4346f3f42a11e4a7e6bb450a0f0b74c0ddae01a10353d2d36676a61808e0318",
     "437f1af29d7d750339b8e0eadd4c906d58871662da1a9c96f43eb42356848b3e"},
    {"c1421ed47c6ec4916eb15f361e42d6fdd936e17aeafe1c093355fed835851369",
     "b166cb5f71939b450d143eed76cef12ed1d5549490ebab0316092d29b2d328ab",
     "666b3d3d23aa41657aa0775bcf2062cb32b80f21de1d0652ac9db1d672228b0a"},
};

static const dili_param sig_params[] = {
    {0, 4, 4, 2, 39, 78, 80, 1 << 17, (DILI_Q - 1) / 88, 192, 32},
    {1, 6, 5, 4, 49, 196, 55, 1 << 19, (DILI_Q - 1) / 32, 128, 48},
    {2, 8, 7, 2, 60, 120, 75, 1 << 19, (DILI_Q - 1) / 32, 128, 64}
};
static const size_t sig_sizes[] = {2420, 3309, 4627};

static void fill(uint8_t *out, size_t len, unsigned int tag) {
    for (size_t i = 0; i < len; ++i)
        out[i] = (uint8_t)(tag + 29u * i + (i >> 2));
}

static void hash_start(KECCAK1600_CTX *hash) {
    sha3_init(hash, 0x06, 256);
}

/* Explicit little-endian coefficient encoding, never struct padding. */
static void hash_polys(KECCAK1600_CTX *hash, const poly *v, int count) {
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < 256; ++j) {
            uint32_t x = (uint32_t)v[i][j];
            uint8_t encoded[4] = {(uint8_t)x, (uint8_t)(x >> 8),
                                  (uint8_t)(x >> 16), (uint8_t)(x >> 24)};
            sha3_update(hash, encoded, sizeof(encoded));
        }
    }
}

static int check_hash(const char *name, int mode, const char *part,
                      KECCAK1600_CTX *hash, const char *expected) {
    uint8_t digest[32];
    char hex[65];
    sha3_final(digest, hash);
    for (size_t i = 0; i < sizeof(digest); ++i)
        snprintf(hex + 2 * i, 3, "%02x", digest[i]);
    if (strcmp(hex, expected) != 0) {
        fprintf(stderr, "%s mode %d %s regression mismatch\nexpected: %s\nactual:   %s\n",
                name, mode, part, expected, hex);
        return 0;
    }
    printf("[PASS] %s mode %d %s regression\n", name, mode, part);
    return 1;
}

static int check_bytes(const char *name, int mode, const char *part,
                       const uint8_t *bytes, size_t len, const char *expected) {
    KECCAK1600_CTX hash;
    hash_start(&hash);
    sha3_update(&hash, bytes, len);
    return check_hash(name, mode, part, &hash, expected);
}

static int test_kyber(void) {
    const int modes[] = {KYBER_MODE_2, KYBER_MODE_3, KYBER_MODE_5};
    uint8_t seed[32], coins[32], message[32], recovered[32];
    uint8_t pk[KYBER_MAX_EK_BYTES], sk[KYBER_MAX_DK_PKE_BYTES];
    uint8_t ct[KYBER_MAX_CT_BYTES];
    int ok = 1;
    for (size_t i = 0; i < 3; ++i) {
        kyber_param param;
        kyber_init_param(&param, modes[i]);
        fill(seed, sizeof(seed), 0x10u + (unsigned int)i);
        fill(coins, sizeof(coins), 0x80u + (unsigned int)i);
        fill(message, sizeof(message), 0x40u + (unsigned int)i);
        if (kyber_pke_keygen(pk, sk, seed, param) != 0 ||
            kyber_pke_encrypt(ct, pk, message, coins, param) != 0 ||
            kyber_pke_decrypt(recovered, sk, ct, param) != 0 ||
            memcmp(message, recovered, sizeof(message)) != 0) {
            fprintf(stderr, "Kyber mode %d key/encrypt/decrypt regression failed\n", modes[i]);
            return 0;
        }
        ok &= check_bytes("Kyber", modes[i], "pk", pk, param.ek_bytes, kyber_expected[i][0]);
        ok &= check_bytes("Kyber", modes[i], "sk", sk, param.dk_pke_bytes, kyber_expected[i][1]);
        ok &= check_bytes("Kyber", modes[i], "ct", ct, param.ct_bytes, kyber_expected[i][2]);
    }
    return ok;
}

static int test_dilithium(void) {
    static const uint8_t message[] = "first-round PQC refactor regression";
    uint8_t seed[32], rnd[32], signature[4627];
    dili_keypair_t keypair;
    int ok = 1;
    for (size_t i = 0; i < 3; ++i) {
        const dili_param param = sig_params[i];
        KECCAK1600_CTX hash;
        fill(seed, sizeof(seed), 0x20u + (unsigned int)i);
        fill(rnd, sizeof(rnd), 0x90u + (unsigned int)i);
        memset(&keypair, 0, sizeof(keypair));
        if (dili_keygen_with_seed(&keypair, param, seed) != 0 ||
            dili_sign_with_rnd(signature, message, sizeof(message) - 1,
                               &keypair.sk, param, rnd) != 0 ||
            dili_verify(signature, message, sizeof(message) - 1,
                        &keypair.pk, param) != 0) {
            fprintf(stderr, "Dilithium mode %d key/sign/verify regression failed\n", param.mode);
            return 0;
        }
        hash_start(&hash);
        sha3_update(&hash, keypair.pk.rho, sizeof(keypair.pk.rho));
        hash_polys(&hash, keypair.pk.t1, param.k);
        ok &= check_hash("Dilithium", param.mode, "pk", &hash, dilithium_expected[i][0]);
        hash_start(&hash);
        sha3_update(&hash, keypair.sk.rho, sizeof(keypair.sk.rho));
        sha3_update(&hash, keypair.sk.K, sizeof(keypair.sk.K));
        sha3_update(&hash, keypair.sk.tr, sizeof(keypair.sk.tr));
        hash_polys(&hash, keypair.sk.s1, param.l);
        hash_polys(&hash, keypair.sk.s2, param.k);
        hash_polys(&hash, keypair.sk.t0, param.k);
        ok &= check_hash("Dilithium", param.mode, "sk", &hash, dilithium_expected[i][1]);
        ok &= check_bytes("Dilithium", param.mode, "sig", signature, sig_sizes[i], dilithium_expected[i][2]);
    }
    return ok;
}

int main(void) {
    int ok = test_kyber();
    ok &= test_dilithium();
    return ok ? 0 : 1;
}
