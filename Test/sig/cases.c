#include "test.h"

const test_sig_alg test_sig_dilithium2_alg = {
    .name = "Dilithium2",
    .security_level = 2,
    .param = {0, 4, 4, 2, 39, 78, 80, (1 << 17), (DILI_Q - 1) / 88, 192, 32},
    .sig_bytes = 2420,
    .rnd_bytes = 32,
    .keygen = dili_keygen,
    .sign = dili_sign,
    .sign_with_rnd = dili_sign_with_rnd,
    .sign_with_rng = dili_sign_with_rng,
    .verify = dili_verify,
};

const test_sig_alg test_sig_dilithium3_alg = {
    .name = "Dilithium3",
    .security_level = 3,
    .param = {1, 6, 5, 4, 49, 196, 55, (1 << 19), (DILI_Q - 1) / 32, 128, 48},
    .sig_bytes = 3309,
    .rnd_bytes = 32,
    .keygen = dili_keygen,
    .sign = dili_sign,
    .sign_with_rnd = dili_sign_with_rnd,
    .sign_with_rng = dili_sign_with_rng,
    .verify = dili_verify,
};

const test_sig_alg test_sig_dilithium5_alg = {
    .name = "Dilithium5",
    .security_level = 5,
    .param = {2, 8, 7, 2, 60, 120, 75, (1 << 19), (DILI_Q - 1) / 32, 128, 64},
    .sig_bytes = 4627,
    .rnd_bytes = 32,
    .keygen = dili_keygen,
    .sign = dili_sign,
    .sign_with_rnd = dili_sign_with_rnd,
    .sign_with_rng = dili_sign_with_rng,
    .verify = dili_verify,
};

#if PQC_DILITHIUM_MODE == 2
const test_sig_alg *test_sig_algs[] = {&test_sig_dilithium2_alg};
#elif PQC_DILITHIUM_MODE == 3
const test_sig_alg *test_sig_algs[] = {&test_sig_dilithium3_alg};
#elif PQC_DILITHIUM_MODE == 5
const test_sig_alg *test_sig_algs[] = {&test_sig_dilithium5_alg};
#else
const test_sig_alg *test_sig_algs[] = {&test_sig_dilithium2_alg, &test_sig_dilithium3_alg, &test_sig_dilithium5_alg};
#endif
const size_t test_sig_alg_count = sizeof(test_sig_algs) / sizeof(test_sig_algs[0]);
