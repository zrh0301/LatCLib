#include "test.h"

const test_enc_alg test_enc_kyber2_pke_alg = {
    .name = "Kyber2 K-PKE",
    .security_level = 2,
    .param = {KYBER_MODE_2, 2, 3, 2, 10, 4, 800, 768, 768},
    .pk_bytes = 800,
    .sk_bytes = 768,
    .ct_bytes = 768,
    .msg_bytes = KYBER_SYMBYTES,
    .seed_bytes = KYBER_SYMBYTES,
    .coins_bytes = KYBER_SYMBYTES,
    .check_pk = kyber_check_ek_modulus,
    .keygen_with_seed = kyber_pke_keygen,
    .encrypt_with_coins = kyber_pke_encrypt,
    .decrypt = kyber_pke_decrypt,
};

const test_enc_alg test_enc_kyber3_pke_alg = {
    .name = "Kyber3 K-PKE",
    .security_level = 3,
    .param = {KYBER_MODE_3, 3, 2, 2, 10, 4, 1184, 1152, 1088},
    .pk_bytes = 1184,
    .sk_bytes = 1152,
    .ct_bytes = 1088,
    .msg_bytes = KYBER_SYMBYTES,
    .seed_bytes = KYBER_SYMBYTES,
    .coins_bytes = KYBER_SYMBYTES,
    .check_pk = kyber_check_ek_modulus,
    .keygen_with_seed = kyber_pke_keygen,
    .encrypt_with_coins = kyber_pke_encrypt,
    .decrypt = kyber_pke_decrypt,
};

const test_enc_alg test_enc_kyber5_pke_alg = {
    .name = "Kyber5 K-PKE",
    .security_level = 5,
    .param = {KYBER_MODE_5, 4, 2, 2, 11, 5, 1568, 1536, 1568},
    .pk_bytes = 1568,
    .sk_bytes = 1536,
    .ct_bytes = 1568,
    .msg_bytes = KYBER_SYMBYTES,
    .seed_bytes = KYBER_SYMBYTES,
    .coins_bytes = KYBER_SYMBYTES,
    .check_pk = kyber_check_ek_modulus,
    .keygen_with_seed = kyber_pke_keygen,
    .encrypt_with_coins = kyber_pke_encrypt,
    .decrypt = kyber_pke_decrypt,
};

#if PQC_KYBER_MODE == 2
const test_enc_alg *test_enc_algs[] = {&test_enc_kyber2_pke_alg};
#elif PQC_KYBER_MODE == 3
const test_enc_alg *test_enc_algs[] = {&test_enc_kyber3_pke_alg};
#elif PQC_KYBER_MODE == 5
const test_enc_alg *test_enc_algs[] = {&test_enc_kyber5_pke_alg};
#else
const test_enc_alg *test_enc_algs[] = {&test_enc_kyber2_pke_alg, &test_enc_kyber3_pke_alg, &test_enc_kyber5_pke_alg};
#endif
const size_t test_enc_alg_count = sizeof(test_enc_algs) / sizeof(test_enc_algs[0]);
