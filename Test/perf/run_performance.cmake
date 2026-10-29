if(NOT DEFINED BUILD_DIR)
    message(FATAL_ERROR "Set BUILD_DIR to the configured build directory")
endif()

get_filename_component(_build_dir "${BUILD_DIR}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_LIST_DIR}/../..")
set(_benchmarks
    "operator_performance|gf16_add"
    "operator_performance|gf16_sub"
    "operator_performance|gf16_mul"
    "operator_performance|gf16_inv"
    "operator_performance|gf16_matrix_mul_4"
    "operator_performance|gf16_matrix_mul_16"
    "operator_performance|gf16_matrix_mul_add_4"
    "operator_performance|gf16_matrix_det_4"
    "operator_performance|gf16_matrix_det_5"
    "operator_performance|fp32_add"
    "operator_performance|fp32_sub"
    "operator_performance|fp32_mul"
    "operator_performance|fp32_mont_reduce"
    "operator_performance|fp32_inv"
    "operator_performance|poly_add_256"
    "operator_performance|poly_sub_256"
    "operator_performance|poly_reduce_256"
    "operator_performance|poly_ntt_256"
    "operator_performance|poly_invntt_256"
    "operator_performance|poly_mul_ntt_256"
    "operator_performance|polyvec_add_4"
    "operator_performance|polyvec_ntt_4"
    "operator_performance|polyvec_inner_4"
    "operator_performance|polyvec_matrix_4"
    "operator_performance|sample_cbd_256"
    "operator_performance|randombytes_168"
    "operator_performance|sha3_256_168"
    "algorithm_snova_performance|keygen"
    "algorithm_snova_performance|sign"
    "algorithm_snova_performance|verify"
    "algorithm_snova_performance|shake"
    "algorithm_snova_performance|genkeys"
    "algorithm_snova_performance|sk_expand"
    "algorithm_snova_performance|pk_expand"
    "algorithm_snova_performance|sign_core"
    "algorithm_snova_performance|verify_core"
    "algorithm_dilithium_performance|_"
    "algorithm_kyber_performance|_")

foreach(_benchmark IN LISTS _benchmarks)
    if(DEFINED BENCH_FILTER AND NOT _benchmark MATCHES "${BENCH_FILTER}")
        continue()
    endif()
    string(REPLACE "|" ";" _parts "${_benchmark}")
    list(GET _parts 0 _program)
    list(LENGTH _parts _part_count)
    set(_command "${_build_dir}/${_program}")
    if(_part_count GREATER 1)
        list(GET _parts 1 _argument)
    else()
        set(_argument "_")
    endif()
    if(NOT _argument STREQUAL "_")
        execute_process(COMMAND "${_command}" "${_argument}"
                        RESULT_VARIABLE _result OUTPUT_VARIABLE _output ERROR_VARIABLE _error)
    else()
        execute_process(COMMAND "${_command}"
                        RESULT_VARIABLE _result OUTPUT_VARIABLE _output ERROR_VARIABLE _error)
    endif()
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "Benchmark failed: ${_benchmark}\n${_error}${_output}")
    endif()

    string(REPLACE "\r\n" "\n" _output "${_output}")
    string(REPLACE "\n" ";" _lines "${_output}")
    foreach(_line IN LISTS _lines)
        if(_line MATCHES "cycles/op")
            message("${_line}")
        endif()
    endforeach()
endforeach()
