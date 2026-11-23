if(NOT DEFINED BUILD_DIR)
    set(BUILD_DIR "@PQC_RUNTIME_OUTPUT_DIRECTORY@")
endif()

get_filename_component(_build_dir "${BUILD_DIR}" ABSOLUTE)
# Generated from the perf cases registered in Test/CMakeLists.txt: every module
# that registers a case is listed here automatically.
set(_benchmarks @PQC_PERF_CASES@)

set(_ran 0)
foreach(_benchmark IN LISTS _benchmarks)
    if(DEFINED BENCH_FILTER AND NOT _benchmark MATCHES "${BENCH_FILTER}")
        continue()
    endif()
    string(REPLACE "|" ";" _parts "${_benchmark}")
    list(GET _parts 0 _program)
    list(LENGTH _parts _part_count)
    set(_command "${_build_dir}/${_program}")
    if(NOT EXISTS "${_command}")
        continue()
    endif()
    math(EXPR _ran "${_ran} + 1")
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

if(_ran EQUAL 0)
    message(FATAL_ERROR
        "No benchmark executable found under ${_build_dir}. "
        "Build the project with PQC_PERF_TESTS=ON first.")
endif()
