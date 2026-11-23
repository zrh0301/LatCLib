# Installation and CMake package export.

include(CMakePackageConfigHelpers)

get_property(pqc_library_targets GLOBAL PROPERTY PQC_LIBRARY_TARGETS)

install(TARGETS ${pqc_library_targets}
    EXPORT pqcTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/inc/"
    DESTINATION include/pqc
    FILES_MATCHING PATTERN "*.h")
install(FILES "${PROJECT_BINARY_DIR}/include/core/pqc_conf.h"
    DESTINATION include/pqc/core)

if(PQC_WITH_DILITHIUM)
    install(FILES alg/dilithium/dili.h DESTINATION include/pqc)
endif()
if(PQC_WITH_KYBER)
    install(FILES alg/kyber/kyber.h DESTINATION include/pqc)
endif()

install(EXPORT pqcTargets
    NAMESPACE pqc::
    DESTINATION lib/cmake/pqc)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/pqc-config.cmake.in"
    "${PROJECT_BINARY_DIR}/pqc-config.cmake"
    INSTALL_DESTINATION lib/cmake/pqc)

write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/pqc-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES
    "${PROJECT_BINARY_DIR}/pqc-config.cmake"
    "${PROJECT_BINARY_DIR}/pqc-config-version.cmake"
    DESTINATION lib/cmake/pqc)
