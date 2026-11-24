# Test and benchmark registration.

# ---- Operator correctness ----

pqc_add_test(operator_fp32_correctness
    SOURCES field/test_fp32.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_fp16_correctness
    SOURCES field/test_fp16.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_gf16_correctness
    SOURCES field/test_gf16.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_context_correctness
    SOURCES algebra/test_contexts.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_poly_correctness
    SOURCES poly/test_poly.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_poly16_correctness
    SOURCES poly/test_poly16.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_ntt_correctness
    SOURCES ntt/test_ntt.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_ntt16_correctness
    SOURCES ntt/test_ntt16.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_ring4171_correctness
    SOURCES poly/test_ring4171.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_hash_correctness
    SOURCES hash/test_hash.c hash/cases.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_sm3_correctness
    SOURCES hash/test_sm3.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_gf_vectors_correctness
    SOURCES field/test_gf_vectors.c LIBS ${PQC_LINK_pqc_gf})
pqc_add_test(operator_gf31_correctness
    SOURCES field/test_gf31.c LIBS ${PQC_LINK_pqc_gf})
pqc_add_test(operator_mod_pow2_correctness
    SOURCES field/test_mod_pow2.c LIBS ${PQC_LINK_pqc_gf})
pqc_add_test(operator_ngcc_gf_correctness
    SOURCES field/test_ngcc_gf.c LIBS ${PQC_LINK_pqc_gf})
pqc_add_test(operator_u256_correctness
    SOURCES util/test_u256.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_poly16_pack_correctness
    SOURCES poly/test_poly16_pack.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_merkle_correctness
    SOURCES tree/test_merkle.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_combinatorics_correctness
    SOURCES util/test_combinatorics.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_ct_correctness
    SOURCES util/test_ct.c LIBS ${PQC_LINK_pqc})
pqc_add_test(operator_bits_correctness
    SOURCES util/test_bits.c LIBS ${PQC_LINK_pqc})

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
    gf256_dot_256 gf256_solve_16 gf256_inv_16 gf256_matvec_16 gf256_quadratic_16
    gf_byte_upper_16 gf_byte_batch_upper_16
    gf16_composite_mul gf256_composite_mul gf16_composite_inv gf256_composite_inv
    gf31_add gf31_mul gf31_inv
    gf31_vector_add_64 gf31_vector_dot_64
    gf31_t2_matrix_mul_4 gf31_t2_matvec_4
    mod_pow2_1024 mod_pow2_round_1024
    gf16_composite_madd_256 gf256_composite_madd_256
    gf16_composite_solve_16 gf16_composite_triangular_16 gf16_composite_matvec_16
    gf16_composite_batch_16 gf256_composite_quad_16
    u256_mod_add u256_mul_div_u32 bytes_be_64
    merkle_treehash_h4 merkle_root_path_h4 merkle_auth_nodes_13
    merkle_build_13 merkle_root_proof_13
    u256_binomial_4721_27 combination_unrank_128_16 u256_reduce_bytes_4
    fp32_add fp32_sub fp32_mul fp32_mont_reduce fp32_inv
    fp16_mont_reduce fp16_barrett_reduce fp16_caddq2
    poly16_add_512 poly16_barrett_512 poly16_tomont_512
    cbd_eta5_512 cbd_eta6_512
    poly16_encode_512 poly16_compress_512 poly16_rej_uniform12_512
    poly16_repeat_encode_512 poly16_repeat_decode_512
    ntt16_128 ntt16_blocks_512 ntt16_product_512
    q4171_montgomery q4171_barrett q4171_ntt_512 q4171_invntt_512
    q4171_pointwise_512 q4171_chknorm_512
    ct_equal_32 ct_select_32
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
    pqc_record_benchmark(algorithm_dilithium_performance _)
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
    pqc_record_benchmark(algorithm_kyber_performance _)
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

if(PQC_WITH_NGCC_VDOO)
    pqc_add_test(ngcc_vdoo_kat
        SOURCES ngcc/test_vdoo_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_vdoo}
        DEFINITIONS
            ${pqc_vdoo_namespace}
            PQC_VDOO_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_SIG_VDOO-128.txt"
        TIMEOUT 900 LABELS ngcc)
    if(NOT PQC_NGCC_VDOO_PARAM STREQUAL "128")
        set_tests_properties(ngcc_vdoo_kat PROPERTIES DISABLED TRUE)
    endif()

    # VDOO hashes the salt into the target, so the whole signature field is
    # authenticated and any byte flip must be rejected.
    pqc_add_test(signature_ngcc_vdoo_correctness
        SOURCES sig/test_signature_api.c
        LIBS ${PQC_LINK_pqc_ngcc_vdoo}
        DEFINITIONS
            PQC_SIG_GETTER=pqc_ngcc_vdoo_algorithm
            PQC_SIG_TAMPER_TRAILER=0
        TIMEOUT 600 LABELS signatures)
endif()

if(PQC_WITH_NGCC_FLEXTREE)
    pqc_add_test(ngcc_flextree_kat
        SOURCES ngcc/test_flextree_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_flextree}
        DEFINITIONS
            sig_keygen=pqc_ngcc_flextree_keygen
            sig_sign=pqc_ngcc_flextree_sign
            sig_verify=pqc_ngcc_flextree_verify
            sig_get_pk_len_bytes=pqc_ngcc_flextree_pk_len
            sig_get_sk_len_bytes=pqc_ngcc_flextree_sk_len
            sig_get_sn_len_bytes=pqc_ngcc_flextree_sig_len
            PQC_FLEXTREE_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_SIG_Flextree-160f.txt"
        TIMEOUT 900 LABELS ngcc)
    # The strict build changes the OTS signer predicate (sign-11-7), so its
    # signatures intentionally differ from the published vectors.
    if(PQC_NGCC_FLEXTREE_STRICT)
        set_tests_properties(ngcc_flextree_kat PROPERTIES DISABLED TRUE)
    endif()

    # FlexTree binds every signature field into the root, including the PORS
    # padding that the reference signer zeroes, so any byte flip is rejected.
    pqc_add_test(signature_ngcc_flextree_correctness
        SOURCES sig/test_signature_api.c
        LIBS ${PQC_LINK_pqc_ngcc_flextree}
        DEFINITIONS
            PQC_SIG_GETTER=pqc_ngcc_flextree_algorithm
            PQC_SIG_TAMPER_TRAILER=0
        TIMEOUT 600 LABELS signatures)
endif()

if(PQC_WITH_NGCC_AIGIS)
    pqc_add_test(ngcc_aigis_kat
        SOURCES ngcc/test_aigis_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_aigis}
        DEFINITIONS
            PARAMS=1
            PQC_AIGIS_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_KEM_Aigis-enc1.txt"
        TIMEOUT 600 LABELS ngcc)
    # The reference entry points are force-prefixed with the same namespace
    # header that the library target uses.
    target_compile_options(ngcc_aigis_kat PRIVATE
        -include "${PROJECT_SOURCE_DIR}/alg/ngcc_aigis/namespace.h")

    pqc_add_test(kem_ngcc_aigis_correctness
        SOURCES kem/test_kem_api.c
        LIBS ${PQC_LINK_pqc_ngcc_aigis}
        DEFINITIONS PQC_KEM_GETTER=pqc_ngcc_aigis_algorithm
        TIMEOUT 600 LABELS kem)

    add_executable(kem_ngcc_aigis_performance perf/test_kem_perf.c)
    target_include_directories(kem_ngcc_aigis_performance
        PRIVATE "${PROJECT_SOURCE_DIR}/Test")
    target_link_libraries(kem_ngcc_aigis_performance
        PRIVATE ${PQC_LINK_pqc_ngcc_aigis})
    target_compile_definitions(kem_ngcc_aigis_performance
        PRIVATE PQC_KEM_GETTER=pqc_ngcc_aigis_algorithm)
    pqc_apply_compile_flags(kem_ngcc_aigis_performance)
    pqc_apply_link_flags(kem_ngcc_aigis_performance)
    foreach(case keypair encaps decaps)
        pqc_add_perf_case(kem_ngcc_aigis_${case}_performance
            kem_ngcc_aigis_performance ${case})
    endforeach()
endif()

if(PQC_WITH_NGCC_AIGIS_SIG)
    pqc_add_test(ngcc_aigis_sig_kat
        SOURCES ngcc/test_aigis_sig_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_aigis_sig}
        DEFINITIONS
            PARAMS=1
            PQC_AIGIS_SIG_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_SIG_Aigis-sig1.txt"
        TIMEOUT 900 LABELS ngcc)
    target_compile_options(ngcc_aigis_sig_kat PRIVATE
        -include "${PROJECT_SOURCE_DIR}/alg/ngcc_aigis_sig/namespace.h")
    # The strict build fixes the challenge sign schedule (sign-01-5), so its
    # signatures intentionally differ from the published vectors.
    if(PQC_NGCC_AIGIS_SIG_STRICT)
        set_tests_properties(ngcc_aigis_sig_kat PROPERTIES DISABLED TRUE)
    endif()

    # The adapter zero-pads the variable-length signature and authenticates the
    # whole fixed field, so any byte flip must be rejected.
    pqc_add_test(signature_ngcc_aigis_sig_correctness
        SOURCES sig/test_signature_api.c
        LIBS ${PQC_LINK_pqc_ngcc_aigis_sig}
        DEFINITIONS
            PQC_SIG_GETTER=pqc_ngcc_aigis_sig_algorithm
            PQC_SIG_TAMPER_TRAILER=0
        TIMEOUT 600 LABELS signatures)
endif()

if(PQC_WITH_NGCC_SCLOUD)
    pqc_add_test(ngcc_scloud_kat
        SOURCES ngcc/test_scloud_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_scloud}
        DEFINITIONS
            SCLOUDPLUS_FAMILY_AES
            SCLOUDPLUS_TIER_REFERENCE
            PQC_SCLOUD_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_KEM_Scloudplus-128-AES-packed10.txt"
        TIMEOUT 900 LABELS ngcc)
    target_compile_options(ngcc_scloud_kat PRIVATE
        -include "${PROJECT_SOURCE_DIR}/alg/ngcc_scloud/namespace.h")

    pqc_add_test(kem_ngcc_scloud_correctness
        SOURCES kem/test_kem_api.c
        LIBS ${PQC_LINK_pqc_ngcc_scloud}
        DEFINITIONS PQC_KEM_GETTER=pqc_ngcc_scloud_algorithm
        TIMEOUT 600 LABELS kem)

    add_executable(kem_ngcc_scloud_performance perf/test_kem_perf.c)
    target_include_directories(kem_ngcc_scloud_performance
        PRIVATE "${PROJECT_SOURCE_DIR}/Test")
    target_link_libraries(kem_ngcc_scloud_performance
        PRIVATE ${PQC_LINK_pqc_ngcc_scloud})
    target_compile_definitions(kem_ngcc_scloud_performance
        PRIVATE PQC_KEM_GETTER=pqc_ngcc_scloud_algorithm)
    pqc_apply_compile_flags(kem_ngcc_scloud_performance)
    pqc_apply_link_flags(kem_ngcc_scloud_performance)
    foreach(case keypair encaps decaps)
        pqc_add_perf_case(kem_ngcc_scloud_${case}_performance
            kem_ngcc_scloud_performance ${case})
    endforeach()
endif()

if(PQC_WITH_NGCC_TSUOV)
    pqc_add_test(ngcc_tsuov_kat
        SOURCES ngcc/test_tsuov_kat.c
        LIBS ${PQC_LINK_pqc_ngcc_tsuov}
        DEFINITIONS
            PQC_TSUOV_KAT_FILE="${PROJECT_SOURCE_DIR}/Test/ngcc/vectors/KAT_SIG_TSUOV-128.txt"
        TIMEOUT 900 LABELS ngcc)
    target_compile_options(ngcc_tsuov_kat PRIVATE
        -include "${PROJECT_SOURCE_DIR}/alg/ngcc_tsuov/namespace.h")
    # Only the 128-bit submitted instance ships a KAT file.
    if(NOT PQC_NGCC_TSUOV_PARAM STREQUAL "128")
        set_tests_properties(ngcc_tsuov_kat PROPERTIES DISABLED TRUE)
    endif()

    # The construction hashes the salt into the target, so the whole signature
    # field is authenticated and any byte flip must be rejected.
    pqc_add_test(signature_ngcc_tsuov_correctness
        SOURCES sig/test_signature_api.c
        LIBS ${PQC_LINK_pqc_ngcc_tsuov}
        DEFINITIONS
            PQC_SIG_GETTER=pqc_ngcc_tsuov_algorithm
            PQC_SIG_TAMPER_TRAILER=0
        TIMEOUT 600 LABELS signatures)
endif()

foreach(algorithm dove vdoo flextree aigis_sig tsuov)
    string(TOUPPER "${algorithm}" algorithm_upper)
    if(PQC_WITH_NGCC_${algorithm_upper})
        set(target signature_ngcc_${algorithm}_performance)
        add_executable(${target} perf/test_signature_perf.c)
        target_include_directories(${target} PRIVATE "${PROJECT_SOURCE_DIR}/Test")
        target_link_libraries(${target} PRIVATE ${PQC_LINK_pqc_ngcc_${algorithm}})
        target_compile_definitions(${target} PRIVATE
            PQC_SIG_GETTER=pqc_ngcc_${algorithm}_algorithm)
        pqc_apply_compile_flags(${target})
        pqc_apply_link_flags(${target})
        foreach(case keygen sign open)
            pqc_add_perf_case(signature_ngcc_${algorithm}_${case}_performance
                ${target} ${case})
        endforeach()
    endif()
endforeach()

if(PQC_WITH_NGCC_DOVE AND PQC_WITH_NGCC_VDOO AND
   PQC_WITH_NGCC_FLEXTREE AND PQC_WITH_NGCC_AIGIS AND
   PQC_WITH_NGCC_AIGIS_SIG AND PQC_WITH_NGCC_SCLOUD AND
   PQC_WITH_NGCC_TSUOV)
    set(pqc_ngcc_link_defs "")
    set(pqc_ngcc_link_libs
        ${PQC_LINK_pqc_ngcc_dove} ${PQC_LINK_pqc_ngcc_vdoo}
        ${PQC_LINK_pqc_ngcc_flextree} ${PQC_LINK_pqc_ngcc_aigis}
        ${PQC_LINK_pqc_ngcc_aigis_sig} ${PQC_LINK_pqc_ngcc_scloud}
        ${PQC_LINK_pqc_ngcc_tsuov})
    if(PQC_WITH_SNOVA_R3)
        list(APPEND pqc_ngcc_link_defs PQC_NGCC_LINK_TEST_SNOVA)
        list(APPEND pqc_ngcc_link_libs ${PQC_LINK_pqc_snova_r3})
    endif()
    pqc_add_test(ngcc_coexistence
        SOURCES ngcc/test_link.c LIBS ${pqc_ngcc_link_libs}
        DEFINITIONS ${pqc_ngcc_link_defs}
        TIMEOUT 600 LABELS ngcc)
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
    # The benchmark list and the target dependencies both come from the
    # registrations above, so a newly added algorithm cannot be silently left
    # out of bench.sh or left unbuilt by the run_performance target.
    get_property(pqc_perf_cases GLOBAL PROPERTY PQC_PERF_CASES)
    get_property(pqc_perf_executables GLOBAL PROPERTY PQC_PERF_EXECUTABLES)
    list(REMOVE_DUPLICATES pqc_perf_executables)
    set(PQC_PERF_CASES "${pqc_perf_cases}")
    configure_file(${PROJECT_SOURCE_DIR}/Test/perf/run_performance.cmake
        "${PROJECT_BINARY_DIR}/run_performance.cmake" @ONLY)

    add_custom_target(run_performance
        COMMAND ${CMAKE_COMMAND} -DBUILD_DIR="${PQC_RUNTIME_OUTPUT_DIRECTORY}"
                -P "${PROJECT_BINARY_DIR}/run_performance.cmake"
        DEPENDS ${pqc_perf_executables}
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}")

    if(UNIX)
        configure_file(${PROJECT_SOURCE_DIR}/Test/perf/bench.sh.in
            "${PROJECT_BINARY_DIR}/bench.sh" @ONLY)
        execute_process(COMMAND chmod +x "${PROJECT_BINARY_DIR}/bench.sh"
            RESULT_VARIABLE bench_chmod_result)
        if(NOT bench_chmod_result EQUAL 0)
            message(FATAL_ERROR
                "Could not make ${PROJECT_BINARY_DIR}/bench.sh executable")
        endif()
    endif()
endif()
