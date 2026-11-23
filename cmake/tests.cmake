# Test registration helpers shared by Test/CMakeLists.txt.

function(pqc_add_test name)
    cmake_parse_arguments(ARG "" "TIMEOUT" "SOURCES;LIBS;DEFINITIONS;LABELS" ${ARGN})
    add_executable(${name} ${ARG_SOURCES})
    target_include_directories(${name} PRIVATE "${PROJECT_SOURCE_DIR}/Test")
    if(ARG_LIBS)
        target_link_libraries(${name} PRIVATE ${ARG_LIBS})
    endif()
    if(ARG_DEFINITIONS)
        target_compile_definitions(${name} PRIVATE ${ARG_DEFINITIONS})
    endif()
    pqc_apply_compile_flags(${name})
    pqc_apply_link_flags(${name})
    if(MINGW)
        target_link_options(${name} PRIVATE -Wl,--stack,16777216)
    endif()
    add_test(NAME ${name} COMMAND ${name})
    if(NOT ARG_TIMEOUT)
        set(ARG_TIMEOUT 60)
    endif()
    set(test_properties TIMEOUT ${ARG_TIMEOUT})
    if(ARG_LABELS)
        list(APPEND test_properties LABELS ${ARG_LABELS})
    endif()
    set_tests_properties(${name} PROPERTIES ${test_properties})
    set_property(GLOBAL APPEND PROPERTY PQC_CORRECTNESS_TARGETS ${name})
endfunction()

# Records one benchmark for bench.sh/run_performance. Kept separate from CTest
# registration so algorithms benchmarked through their own test executable can
# join the same list.
function(pqc_record_benchmark executable case)
    set_property(GLOBAL APPEND PROPERTY PQC_PERF_CASES "${executable}|${case}")
    set_property(GLOBAL APPEND PROPERTY PQC_PERF_EXECUTABLES ${executable})
endfunction()

# Registers one performance case of an already declared executable.
function(pqc_add_perf_case name executable case)
    add_test(NAME ${name} COMMAND ${executable} ${case})
    set_property(GLOBAL APPEND PROPERTY PQC_PERF_TEST_NAMES ${name})
    pqc_record_benchmark(${executable} ${case})
endfunction()

# Registers a whole executable as one performance benchmark.
function(pqc_add_benchmark name executable)
    add_test(NAME ${name} COMMAND ${executable})
    set_property(GLOBAL APPEND PROPERTY PQC_PERF_TEST_NAMES ${name})
    pqc_record_benchmark(${executable} _)
endfunction()

# Applies the performance label and the enabled/disabled state in one place.
function(pqc_finalize_performance_tests)
    get_property(perf_tests GLOBAL PROPERTY PQC_PERF_TEST_NAMES)
    if(NOT perf_tests)
        return()
    endif()
    if(PQC_PERF_TESTS)
        set_tests_properties(${perf_tests} PROPERTIES
            TIMEOUT 600 LABELS performance)
    else()
        set_tests_properties(${perf_tests} PROPERTIES
            TIMEOUT 600 LABELS performance DISABLED TRUE)
    endif()
endfunction()
