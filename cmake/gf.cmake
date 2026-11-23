# GF(16)/GF(256) arithmetic module.
#
# The portable scalar sources and the platform kernels use different views of
# "gf16.h": the scalar sources include the public inc/gf/gf16.h, while the
# kernels include the architecture-specific internal/gf16.h. They are
# therefore built as two OBJECT libraries that end up in the same shared and
# static variants.

if(PQC_GF_BACKEND_EFFECTIVE MATCHES "^(AVX2|GFNI)$")
    set(pqc_gf_platform x86)
    set(pqc_gf_flags -mavx2 -mpclmul)
elseif(PQC_GF_BACKEND_EFFECTIVE STREQUAL "NEON")
    set(pqc_gf_platform armv8)
    set(pqc_gf_flags -march=armv8-a+crypto)
else()
    message(FATAL_ERROR "Unknown PQC_GF_BACKEND=${PQC_GF_BACKEND_EFFECTIVE}")
endif()

file(GLOB pqc_gf_kernel_sources CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/gf/${pqc_gf_platform}/core/*.c"
    "${PROJECT_SOURCE_DIR}/src/gf/${pqc_gf_platform}/adapters/*.c"
    "${PROJECT_SOURCE_DIR}/src/gf/${pqc_gf_platform}/workloads/*.c")

set(pqc_gf_definitions "")
if(PQC_GF_BACKEND_EFFECTIVE STREQUAL "GFNI")
    list(APPEND pqc_gf_flags -mgfni)
    list(APPEND pqc_gf_definitions GF16_ENABLE_GFNI=1)
else()
    # The non-GFNI packed source contains upstream placeholders.
    list(FILTER pqc_gf_kernel_sources EXCLUDE REGEX "gf16_gfni_pack\\.c$")
endif()

set(pqc_gf_kernel_includes
    "${PROJECT_SOURCE_DIR}/src/gf/${pqc_gf_platform}/internal")

add_library(pqc_gf_kernels_obj OBJECT
    ${pqc_gf_kernel_sources}
    src/gf/platform_bridge.c)
set_target_properties(pqc_gf_kernels_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(pqc_gf_kernels_obj PRIVATE ${pqc_gf_kernel_includes})
target_compile_options(pqc_gf_kernels_obj PRIVATE ${pqc_gf_flags})
if(pqc_gf_definitions)
    target_compile_definitions(pqc_gf_kernels_obj PRIVATE ${pqc_gf_definitions})
endif()
pqc_apply_compile_flags(pqc_gf_kernels_obj)

pqc_add_variants(pqc_gf
    SOURCES
        src/gf/gf16.c
        src/gf/gf16_matrix.c
        src/gf/common/vector.c
        src/gf/common/matrix.c
        src/gf/common/gf256_solve.c
        src/gf/common/byte_matrix.c
        src/gf/tower/coordinates.c
        src/gf/tower/composite.c
        src/gf/tower/composite_matrix.c
    PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
    EXTRA_OBJECTS pqc_gf_kernels_obj)
