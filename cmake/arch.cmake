# Architecture and backend detection.
# The selected backends require GCC or Clang, because they use -mavx2, -mgfni
# or -march=armv8-a+crypto.

if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR
        "The optimized GF and signature backends require GCC or Clang; "
        "CMAKE_C_COMPILER_ID is ${CMAKE_C_COMPILER_ID}")
endif()

if(PQC_ARCH STREQUAL "AUTO")
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" PQC_PROCESSOR)
    if(PQC_PROCESSOR MATCHES "^(x86_64|amd64)$")
        set(PQC_ARCH_EFFECTIVE "X64")
    elseif(PQC_PROCESSOR MATCHES "^(aarch64|arm64)$")
        set(PQC_ARCH_EFFECTIVE "A64")
    else()
        message(FATAL_ERROR
            "Optimized GF/signatures require x86-64 or AArch64; target: "
            "${CMAKE_SYSTEM_PROCESSOR}")
    endif()
else()
    set(PQC_ARCH_EFFECTIVE "${PQC_ARCH}")
endif()

if(PQC_ARCH_EFFECTIVE STREQUAL "X64")
    set(PQC_DEFAULT_BACKEND "AVX2")
    set(PQC_ALLOWED_BACKENDS AVX2 GFNI)
elseif(PQC_ARCH_EFFECTIVE STREQUAL "A64")
    set(PQC_DEFAULT_BACKEND "NEON")
    set(PQC_ALLOWED_BACKENDS NEON)
else()
    message(FATAL_ERROR "Unknown PQC_ARCH=${PQC_ARCH}, choose AUTO, X64 or A64")
endif()

function(pqc_resolve_backend option_name out_var)
    if(${option_name} STREQUAL "AUTO")
        set(backend "${PQC_DEFAULT_BACKEND}")
    else()
        set(backend "${${option_name}}")
    endif()
    list(FIND PQC_ALLOWED_BACKENDS "${backend}" backend_index)
    if(backend_index EQUAL -1)
        string(REPLACE ";" " " allowed_text "${PQC_ALLOWED_BACKENDS}")
        message(FATAL_ERROR
            "${option_name}=${backend} is unsupported for ${PQC_ARCH_EFFECTIVE}; "
            "choose AUTO or one of: ${allowed_text}")
    endif()
    set(${out_var} "${backend}" PARENT_SCOPE)
endfunction()

pqc_resolve_backend(PQC_BACKEND PQC_BACKEND_EFFECTIVE)
pqc_resolve_backend(PQC_GF_BACKEND PQC_GF_BACKEND_EFFECTIVE)

# SNOVA R2 has no GFNI implementation.
set(PQC_SNOVA_R2_BACKEND_EFFECTIVE "${PQC_BACKEND_EFFECTIVE}")
if(PQC_BACKEND_EFFECTIVE STREQUAL "GFNI")
    set(PQC_SNOVA_R2_BACKEND_EFFECTIVE "AVX2")
    message(STATUS "SNOVA R2 has no GFNI backend; using AVX2")
endif()

# AES backend: AUTO probes the host, C is portable, NI forces AES-NI.
set(PQC_AES_EFFECTIVE "${PQC_AES}")
if(PQC_AES_EFFECTIVE STREQUAL "NI")
    if(NOT PQC_ARCH_EFFECTIVE STREQUAL "X64")
        message(FATAL_ERROR "PQC_AES=NI requires PQC_ARCH=X64")
    endif()
endif()
if(PQC_AES_EFFECTIVE STREQUAL "AUTO")
    if(PQC_ARCH_EFFECTIVE STREQUAL "X64" AND NOT CMAKE_CROSSCOMPILING)
        include(CheckCSourceRuns)
        set(CMAKE_REQUIRED_FLAGS "-maes -mavx2")
        set(PQC_AES_NI_PROBE "#if !defined(__AES__) || !defined(__AVX2__)
#error AES-NI is not available
#endif
int main(void) {
    return __builtin_cpu_supports(\"aes\") && __builtin_cpu_supports(\"avx2\") ? 0 : 1;
}")
        check_c_source_runs("${PQC_AES_NI_PROBE}" PQC_AES_NI_WORKS)
        unset(CMAKE_REQUIRED_FLAGS)
        if(PQC_AES_NI_WORKS)
            set(PQC_AES_EFFECTIVE "NI")
        else()
            set(PQC_AES_EFFECTIVE "C")
        endif()
    else()
        set(PQC_AES_EFFECTIVE "C")
    endif()
endif()
message(STATUS "AES backend: ${PQC_AES_EFFECTIVE}")
