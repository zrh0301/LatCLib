cmake_minimum_required(VERSION 3.10)
project(MyPQC C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
option(PQC_NATIVE_OPTIMIZATION "Optimize for the build host CPU" OFF)
option(PQC_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)
set(PQC_DILITHIUM_MODE 0 CACHE STRING "Dilithium parameter set: 0=all, 2, 3, or 5")
set(PQC_KYBER_MODE 0 CACHE STRING "Kyber parameter set: 0=all, 2, 3, or 5")
set_property(CACHE PQC_DILITHIUM_MODE PROPERTY STRINGS 0 2 3 5)
set_property(CACHE PQC_KYBER_MODE PROPERTY STRINGS 0 2 3 5)
if(NOT PQC_DILITHIUM_MODE MATCHES "^(0|2|3|5)$")
    message(FATAL_ERROR "PQC_DILITHIUM_MODE must be 0, 2, 3, or 5")
endif()
if(NOT PQC_KYBER_MODE MATCHES "^(0|2|3|5)$")
    message(FATAL_ERROR "PQC_KYBER_MODE must be 0, 2, 3, or 5")
endif()

add_library(pqc STATIC
    src/fp32.c src/ring.c src/poly.c src/poly_ntt.c src/polyvec.c
    src/ntt_table.c src/ntt_compat.c src/poly_sample.c src/rand.c src/keccak.c
    alg/dilithium/dili.c alg/dilithium/dili_sample.c alg/dilithium/pack.c
    alg/kyber/kyber.c alg/kyber/kyber_sample.c alg/kyber/kyber_pack.c)
target_include_directories(pqc PUBLIC inc alg/dilithium alg/kyber)
if(NOT MSVC)
    target_link_libraries(pqc PUBLIC m)
endif()
if(WIN32)
    target_link_libraries(pqc PUBLIC bcrypt)
endif()
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pqc PRIVATE -Wall -Wextra)
    if(PQC_NATIVE_OPTIMIZATION)
        target_compile_options(pqc PRIVATE -march=native -mtune=native)
    endif()
    if(PQC_ENABLE_SANITIZERS)
        target_compile_options(pqc PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(pqc PUBLIC -fsanitize=address,undefined)
    endif()
elseif(MSVC)
    target_compile_options(pqc PRIVATE /W4)
endif()

include(CTest)
if(BUILD_TESTING)
    function(add_pqc_test name source)
        add_executable(${name} ${source})
        target_include_directories(${name} PRIVATE Test)
        target_link_libraries(${name} PRIVATE pqc)
        add_test(NAME ${name} COMMAND ${name})
        set_tests_properties(${name} PROPERTIES TIMEOUT 60)
    endfunction()

    # Operator correctness.
    add_pqc_test(operator_fp32_correctness Test/field/test_fp32.c)
    add_pqc_test(operator_context_correctness Test/algebra/test_contexts.c)
    add_pqc_test(operator_poly_correctness Test/poly/test_poly.c)
    add_pqc_test(operator_ntt_correctness Test/ntt/test_ntt.c)
    add_pqc_test(operator_sampling_correctness Test/sampling/test_sampling.c)
    add_pqc_test(operator_pack_correctness Test/pack/test_pack.c)
    add_pqc_test(operator_hash_correctness Test/hash/test_hash.c)
    target_sources(operator_hash_correctness PRIVATE Test/hash/cases.c)
    add_pqc_test(operator_rng_correctness Test/rng/test_rng.c)

    # Dedicated operator performance executable.
    add_pqc_test(operator_performance Test/perf/test_operator_perf.c)

    target_sources(operator_pack_correctness PRIVATE
        Test/pack/cases.c Test/pack/dilithium.c Test/pack/kyber.c)
    target_compile_definitions(operator_rng_correctness PRIVATE PQC_TEST_WRAP_GETRANDOM)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_link_options(operator_rng_correctness PRIVATE -Wl,--wrap=getrandom)
    endif()

    # Algorithm correctness and dedicated performance executables.
    add_pqc_test(algorithm_dilithium_correctness Test/sig/test_sig.c)
    target_sources(algorithm_dilithium_correctness PRIVATE Test/sig/cases.c)
    target_compile_definitions(algorithm_dilithium_correctness PRIVATE PQC_DILITHIUM_MODE=${PQC_DILITHIUM_MODE})
    add_pqc_test(algorithm_kyber_correctness Test/enc/test_enc.c)
    target_sources(algorithm_kyber_correctness PRIVATE Test/enc/cases.c)
    target_compile_definitions(algorithm_kyber_correctness PRIVATE PQC_KYBER_MODE=${PQC_KYBER_MODE})
    add_pqc_test(algorithm_dilithium_performance Test/sig/test_sig.c)
    target_sources(algorithm_dilithium_performance PRIVATE Test/sig/cases.c)
    target_compile_definitions(algorithm_dilithium_performance PRIVATE PQC_DILITHIUM_MODE=${PQC_DILITHIUM_MODE} PQC_PERF_ONLY)
    add_pqc_test(algorithm_kyber_performance Test/enc/test_enc.c)
    target_sources(algorithm_kyber_performance PRIVATE Test/enc/cases.c)
    target_compile_definitions(algorithm_kyber_performance PRIVATE PQC_KYBER_MODE=${PQC_KYBER_MODE} PQC_PERF_ONLY)
    if(PQC_DILITHIUM_MODE EQUAL 0 AND PQC_KYBER_MODE EQUAL 0)
        add_pqc_test(algorithm_regression Test/regression/test_regression.c)
    endif()

    add_custom_target(run_test
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS operator_fp32_correctness operator_context_correctness
                operator_poly_correctness operator_ntt_correctness
                operator_sampling_correctness operator_pack_correctness
                operator_hash_correctness operator_rng_correctness operator_performance
                algorithm_dilithium_correctness algorithm_kyber_correctness
                algorithm_dilithium_performance algorithm_kyber_performance
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endif()
