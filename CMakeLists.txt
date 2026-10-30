cmake_minimum_required(VERSION 3.10)
project(MyPQC C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
option(PQC_NATIVE_OPTIMIZATION "Optimize for the build host CPU" OFF)
option(PQC_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)
option(PQC_ENABLE_PERFORMANCE_TESTS "Enable performance tests in CTest" OFF)
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
    src/fp32.c src/gf16.c src/gf16_matrix.c src/ring.c src/poly.c src/poly_ntt.c src/polyvec.c
    src/ntt_table.c src/ntt_compat.c src/poly_sample.c src/rand.c src/secure_clear.c src/keccak.c
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

add_library(pqc_snova STATIC
    alg/snova/snova_opt_16.c
    alg/snova/symmetric.c
    alg/snova/aes.c
    alg/snova/sign.c)
target_include_directories(pqc_snova PUBLIC inc alg/snova)
target_compile_definitions(pqc_snova PUBLIC SNOVA_OPT=10)
target_link_libraries(pqc_snova PUBLIC pqc)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pqc_snova PRIVATE -Wall -Wextra)
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
    add_pqc_test(operator_gf16_correctness Test/field/test_gf16.c)
    add_pqc_test(operator_context_correctness Test/algebra/test_contexts.c)
    add_pqc_test(operator_poly_correctness Test/poly/test_poly.c)
    add_pqc_test(operator_ntt_correctness Test/ntt/test_ntt.c)
    add_pqc_test(operator_sampling_correctness Test/sampling/test_sampling.c)
    add_pqc_test(operator_pack_correctness Test/pack/test_pack.c)
    add_pqc_test(operator_hash_correctness Test/hash/test_hash.c)
    target_sources(operator_hash_correctness PRIVATE Test/hash/cases.c)
    add_pqc_test(operator_rng_correctness Test/rng/test_rng.c)

    # Each operator runs in its own CTest entry and reports its own timing.
    add_executable(operator_performance Test/perf/test_operator_perf.c)
    target_include_directories(operator_performance PRIVATE Test)
    target_link_libraries(operator_performance PRIVATE pqc)
    set(PQC_PERFORMANCE_TESTS)
    set(PQC_OPERATOR_PERF_CASES
        gf16_add gf16_sub gf16_mul gf16_inv
        gf16_matrix_mul_4 gf16_matrix_mul_16 gf16_matrix_mul_add_4
        gf16_matrix_det_4 gf16_matrix_det_5
        fp32_add fp32_sub fp32_mul fp32_mont_reduce fp32_inv
        poly_add_256 poly_sub_256 poly_reduce_256 poly_ntt_256
        poly_invntt_256 poly_mul_ntt_256
        polyvec_add_4 polyvec_ntt_4 polyvec_inner_4 polyvec_matrix_4
        sample_cbd_256 randombytes_168 sha3_256_168)
    set(PQC_PERFORMANCE_COMMANDS)
    foreach(operator_case IN LISTS PQC_OPERATOR_PERF_CASES)
        set(test_name "operator_${operator_case}_performance")
        add_test(NAME ${test_name} COMMAND operator_performance ${operator_case})
        list(APPEND PQC_PERFORMANCE_TESTS ${test_name})
        list(APPEND PQC_PERFORMANCE_COMMANDS
             COMMAND $<TARGET_FILE:operator_performance> ${operator_case})
    endforeach()

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
    add_pqc_test(algorithm_snova_correctness Test/snova/test_snova.c)
    target_link_libraries(algorithm_snova_correctness PRIVATE pqc_snova)
    add_pqc_test(algorithm_dilithium_performance Test/sig/test_sig.c)
    target_sources(algorithm_dilithium_performance PRIVATE Test/sig/cases.c)
    target_compile_definitions(algorithm_dilithium_performance PRIVATE PQC_DILITHIUM_MODE=${PQC_DILITHIUM_MODE} PQC_PERF_ONLY)
    add_pqc_test(algorithm_kyber_performance Test/enc/test_enc.c)
    target_sources(algorithm_kyber_performance PRIVATE Test/enc/cases.c)
    target_compile_definitions(algorithm_kyber_performance PRIVATE PQC_KYBER_MODE=${PQC_KYBER_MODE} PQC_PERF_ONLY)
    add_executable(algorithm_snova_performance Test/perf/test_snova_perf.c)
    target_include_directories(algorithm_snova_performance PRIVATE Test)
    target_link_libraries(algorithm_snova_performance PRIVATE pqc_snova)
    foreach(snova_case keygen sign verify shake genkeys sk_expand pk_expand sign_core verify_core)
        set(test_name "algorithm_snova_${snova_case}_performance")
        add_test(NAME ${test_name} COMMAND algorithm_snova_performance ${snova_case})
        list(APPEND PQC_PERFORMANCE_TESTS ${test_name})
        list(APPEND PQC_PERFORMANCE_COMMANDS
             COMMAND $<TARGET_FILE:algorithm_snova_performance> ${snova_case})
    endforeach()
    list(APPEND PQC_PERFORMANCE_TESTS
         algorithm_dilithium_performance algorithm_kyber_performance)
    list(APPEND PQC_PERFORMANCE_COMMANDS
         COMMAND $<TARGET_FILE:algorithm_dilithium_performance>
         COMMAND $<TARGET_FILE:algorithm_kyber_performance>)
    if(PQC_ENABLE_PERFORMANCE_TESTS)
        set_tests_properties(${PQC_PERFORMANCE_TESTS}
                             PROPERTIES TIMEOUT 600 LABELS performance)
    else()
        set_tests_properties(${PQC_PERFORMANCE_TESTS}
                             PROPERTIES DISABLED TRUE LABELS performance)
    endif()
    if(PQC_DILITHIUM_MODE EQUAL 0 AND PQC_KYBER_MODE EQUAL 0)
        add_pqc_test(algorithm_regression Test/regression/test_regression.c)
    endif()

    add_custom_target(run_test
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -LE performance
        DEPENDS operator_fp32_correctness operator_context_correctness
                operator_gf16_correctness
                operator_poly_correctness operator_ntt_correctness
                operator_sampling_correctness operator_pack_correctness
                operator_hash_correctness operator_rng_correctness
                algorithm_dilithium_correctness algorithm_kyber_correctness
                algorithm_snova_correctness
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

    add_custom_target(run_performance
        ${PQC_PERFORMANCE_COMMANDS}
        DEPENDS operator_performance algorithm_dilithium_performance
                algorithm_kyber_performance algorithm_snova_performance
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

    if(PQC_ENABLE_PERFORMANCE_TESTS AND UNIX)
        configure_file(Test/perf/bench.sh.in "${CMAKE_BINARY_DIR}/bench.sh" @ONLY)
        configure_file(Test/perf/run_performance.cmake
                       "${CMAKE_BINARY_DIR}/run_performance.cmake" COPYONLY)
        execute_process(COMMAND chmod +x "${CMAKE_BINARY_DIR}/bench.sh"
                        RESULT_VARIABLE bench_chmod_result)
        if(NOT bench_chmod_result EQUAL 0)
            message(FATAL_ERROR "Could not make ${CMAKE_BINARY_DIR}/bench.sh executable")
        endif()
    endif()
endif()
