#include "test.h"

extern const test_pack_suite test_pack_dilithium_suite;
extern const test_pack_suite test_pack_kyber_suite;

const test_pack_suite *test_pack_suites[] = {&test_pack_dilithium_suite, &test_pack_kyber_suite};
const size_t test_pack_suite_count = sizeof(test_pack_suites) / sizeof(test_pack_suites[0]);
