# Test and benchmark registration.

# ---- Operator correctness ----

pqc_add_test(operator_fp32_correctness
    SOURCES field/test_fp32.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_gf16_correctness
    SOURCES field/test_gf16.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_context_correctness
    SOURCES algebra/test_contexts.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_poly_correctness
    SOURCES poly/test_poly.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_ntt_correctness
    SOURCES ntt/test_ntt.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_hash_correctness
    SOURCES hash/test_hash.c hash/cases.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_sm3_correctness
    SOURCES hash/test_sm3.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_gf_vectors_correctness
    SOURCES field/test_gf_vectors.c LIBS ${PQC_LINK_pqc_gf})

# ---- Operator performance ----

add_executable(operator_performance perf/test_operator_perf.c)
target_include_directories(operator_performance PRIVATE "${PROJECT_SOURCE_DIR}/Test")
target_link_libraries(operator_performance PRIVATE ${PQC_LINK_pqc})
pqc_apply_compile_flags(operator_performance)
pqc_apply_link_flags(operator_performance)

set(PQC_OPERATOR_PERF_CASES
    gf16_add gf16_sub gf16_mul gf16_inv
    gf16_matrix_mul_4 gf16_matrix_mul_16 gf16_matrix_mul_add_4
    gf16_matrix_det_4 gf16_matrix_det_5
    fp32_add fp32_sub fp32_mul fp32_mont_reduce fp32_inv
    poly_add_256 poly_sub_256 poly_reduce_256 poly_ntt_256
    poly_invntt_256 poly_mul_ntt_256
    polyvec_add_4 polyvec_ntt_4 polyvec_inner_4 polyvec_matrix_4
    sample_cbd_256 randombytes_168 sha3_256_168)
foreach(operator_case IN LISTS PQC_OPERATOR_PERF_CASES)
    pqc_add_perf_case(operator_${operator_case}_performance
        operator_performance ${operator_case})
endforeach()

# ---- Dilithium / Kyber ----

if(PQC_WITH_DILITHIUM)
    pqc_add_test(algorithm_dilithium_correctness
        SOURCES sig/test_sig.c sig/cases.c
        LIBS ${PQC_LINK_pqc}
        DEFINITIONS PQC_DILITHIUM_MODE=${PQC_DILITHIUM_MODE})
    pqc_add_test(algorithm_dilithium_performance
        SOURCES sig/test_sig.c sig/cases.c
        LIBS ${PQC_LINK_pqc}
        DEFINITIONS PQC_DILITHIUM_MODE=${PQC_DILITHIUM_MODE} PQC_PERF_ONLY)
    set_property(GLOBAL APPEND PROPERTY PQC_PERF_TEST_NAMES
        algorithm_dilithium_performance)
endif()

if(PQC_WITH_KYBER)
    pqc_add_test(algorithm_kyber_correctness
        SOURCES enc/test_enc.c enc/cases.c
        LIBS ${PQC_LINK_pqc}
        DEFINITIONS PQC_KYBER_MODE=${PQC_KYBER_MODE})
    pqc_add_test(algorithm_kyber_performance
        SOURCES enc/test_enc.c enc/cases.c
        LIBS ${PQC_LINK_pqc}
        DEFINITIONS PQC_KYBER_MODE=${PQC_KYBER_MODE} PQC_PERF_ONLY)
    set_property(GLOBAL APPEND PROPERTY PQC_PERF_TEST_NAMES
        algorithm_kyber_performance)
endif()

# These tests exercise both lattice algorithms, so they need both modules.
if(PQC_WITH_DILITHIUM AND PQC_WITH_KYBER)
    pqc_add_test(operator_sampling_correctness
        SOURCES sampling/test_sampling.c LIBS ${PQC_LINK_pqc})

    pqc_add_test(operator_pack_correctness
        SOURCES pack/test_pack.c pack/cases.c
                pack/dilithium.c pack/kyber.c
        LIBS ${PQC_LINK_pqc})

    if(PQC_STLIB)
        pqc_add_test(operator_rng_correctness
            SOURCES rng/test_rng.c LIBS ${PQC_LINK_pqc})
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_definitions(operator_rng_correctness PRIVATE
                PQC_TEST_WRAP_GETRANDOM)
            target_link_options(operator_rng_correctness PRIVATE
                -Wl,--wrap=getrandom)
        endif()
    else()
        message(WARNING
            "operator_rng_correctness needs the static library for --wrap=getrandom; "
            "it is skipped because PQC_STLIB=OFF")
    endif()

    if(PQC_DILITHIUM_MODE EQUAL 0 AND PQC_KYBER_MODE EQUAL 0)
        pqc_add_test(algorithm_regression
            SOURCES regression/test_regression.c LIBS ${PQC_LINK_pqc})
    endif()
endif()

# ---- SNOVA R3 internal API test ----

if(PQC_WITH_SNOVA_R3)
    pqc_add_test(algorithm_snova_correctness
        SOURCES sig/internal/test_snova_r3.c
        LIBS ${PQC_LINK_pqc_snova_r3})

    add_executable(algorithm_snova_performance perf/test_snova_perf.c)
    target_include_directories(algorithm_snova_performance PRIVATE
        "${PROJECT_SOURCE_DIR}/Test")
    target_link_libraries(algorithm_snova_performance PRIVATE
        ${PQC_LINK_pqc_snova_r3})
    pqc_apply_compile_flags(algorithm_snova_performance)
    pqc_apply_link_flags(algorithm_snova_performance)
    foreach(snova_case keygen sign verify shake genkeys sk_expand
                       pk_expand sign_core verify_core)
        pqc_add_perf_case(algorithm_snova_${snova_case}_performance
            algorithm_snova_performance ${snova_case})
    endforeach()
endif()

# ---- NGCC round-1 candidates ----

if(PQC_WITH_NGCC_DOVE)
    pqc_add_test(ngcc_dove_kat
        SOURCES ngcc/test_dove_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_dove}
        DEFINITIONS
            PQC_DOVE_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_SIG_dove_pkc_skc_128.txt"
        TIMEOUT 900 LABELS ngcc)
    if(NOT PQC_NGCC_DOVE_PARAM STREQUAL "128" OR PQC_NGCC_DOVE_BIND_SALT)
        set_tests_properties(ngcc_dove_kat PROPERTIES DISABLED TRUE)
    endif()

    # The submitted construction appends a salt that its verifier never reads
    # (sign-09-1), so that trailer is excluded from the tampering checks.
    # PQC_NGCC_DOVE_BIND_SALT binds it into the target and covers it.
    if(PQC_NGCC_DOVE_BIND_SALT)
        set(pqc_dove_trailer 0)
    elseif(PQC_NGCC_DOVE_PARAM STREQUAL "128")
        set(pqc_dove_trailer 24)
    elseif(PQC_NGCC_DOVE_PARAM STREQUAL "256")
        set(pqc_dove_trailer 40)
    else()
        set(pqc_dove_trailer 72)
    endif()
    pqc_add_test(signature_ngcc_dove_correctness
        SOURCES sig/test_signature_api.c
        LIBS ${PQC_LINK_pqc_ngcc_dove}
        DEFINITIONS
            PQC_SIG_GETTER=pqc_ngcc_dove_algorithm
            PQC_SIG_TAMPER_TRAILER=${pqc_dove_trailer}
        TIMEOUT 600 LABELS signatures)
endif()

# ---- Unified signature API ----

foreach(algorithm snova_r2 snova_r3 mayo mqom)
    string(TOUPPER "${algorithm}" algorithm_upper)
    if(PQC_WITH_${algorithm_upper})
        pqc_add_test(signature_${algorithm}_correctness
            SOURCES sig/test_signature_api.c
            LIBS ${PQC_LINK_pqc_${algorithm}}
            DEFINITIONS PQC_SIG_GETTER=pqc_${algorithm}_algorithm
            TIMEOUT 600 LABELS signatures)
    endif()
endforeach()

if(PQC_WITH_SNOVA_R2 AND PQC_WITH_SNOVA_R3 AND
   PQC_WITH_MAYO AND PQC_WITH_MQOM)
    pqc_add_test(signature_coexistence
        SOURCES sig/test_signature_link.c
        LIBS ${PQC_LINK_pqc_snova_r2} ${PQC_LINK_pqc_snova_r3}
             ${PQC_LINK_pqc_mayo} ${PQC_LINK_pqc_mqom}
        TIMEOUT 600 LABELS signatures)
endif()

pqc_finalize_performance_tests()

# ---- Aggregate targets ----

get_property(pqc_correctness_targets GLOBAL PROPERTY PQC_CORRECTNESS_TARGETS)
add_custom_target(run_test
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -LE performance
    DEPENDS ${pqc_correctness_targets}
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}")

if(PQC_PERF_TESTS)
    add_custom_target(run_performance
        COMMAND ${CMAKE_COMMAND} -DBUILD_DIR="${PQC_RUNTIME_OUTPUT_DIRECTORY}"
                -P "${PROJECT_BINARY_DIR}/run_performance.cmake"
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}")

    if(UNIX)
        configure_file(${PROJECT_SOURCE_DIR}/Test/perf/bench.sh.in
            "${PROJECT_BINARY_DIR}/bench.sh" @ONLY)
        configure_file(${PROJECT_SOURCE_DIR}/Test/perf/run_performance.cmake
            "${PROJECT_BINARY_DIR}/run_performance.cmake" @ONLY)
        execute_process(COMMAND chmod +x "${PROJECT_BINARY_DIR}/bench.sh"
            RESULT_VARIABLE bench_chmod_result)
        if(NOT bench_chmod_result EQUAL 0)
            message(FATAL_ERROR
                "Could not make ${PROJECT_BINARY_DIR}/bench.sh executable")
        endif()
    endif()
endif()
