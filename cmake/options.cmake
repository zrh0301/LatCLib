# Build options, module selection and parameter validation.
# All public switches use the PQC_ prefix, following the relic-toolkit style.

set(PQC_ARCH "AUTO" CACHE STRING "Target architecture: AUTO, X64 or A64")
set_property(CACHE PQC_ARCH PROPERTY STRINGS AUTO X64 A64)

set(PQC_BACKEND "AUTO" CACHE STRING "Signature backend: AUTO, AVX2, GFNI or NEON")
set_property(CACHE PQC_BACKEND PROPERTY STRINGS AUTO AVX2 GFNI NEON)

set(PQC_GF_BACKEND "AUTO" CACHE STRING "GF backend: AUTO, AVX2, GFNI or NEON")
set_property(CACHE PQC_GF_BACKEND PROPERTY STRINGS AUTO AVX2 GFNI NEON)

set(PQC_AES "AUTO" CACHE STRING "AES backend: AUTO, C or NI")
set_property(CACHE PQC_AES PROPERTY STRINGS AUTO C NI)

set(PQC_WITH "ALL" CACHE STRING
    "Modules to build: ALL or a list of DILITHIUM, KYBER, SNOVA_R2, SNOVA_R3, MAYO, MQOM, NGCC_DOVE")

option(PQC_SHLIB "Build shared libraries" ON)
option(PQC_STLIB "Build static libraries" ON)
option(PQC_TESTS "Build automated tests" ON)
option(PQC_PERF_TESTS "Enable performance tests in CTest" OFF)
option(PQC_INSTALL "Generate installation rules" ON)
option(PQC_NATIVE "Optimize for the build host CPU" OFF)
option(PQC_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)
option(PQC_STRIP "Remove unused code and data" OFF)
set(PQC_CFLAGS "" CACHE STRING "Extra compiler flags appended to every target")

if(NOT PQC_SHLIB AND NOT PQC_STLIB)
    message(FATAL_ERROR "At least one of PQC_SHLIB or PQC_STLIB must be enabled")
endif()

# Expand PQC_WITH into one boolean per module.
set(PQC_ALL_MODULES DILITHIUM KYBER SNOVA_R2 SNOVA_R3 MAYO MQOM NGCC_DOVE)
string(TOUPPER "${PQC_WITH}" PQC_WITH_UPPER)
string(REPLACE "," ";" PQC_WITH_UPPER "${PQC_WITH_UPPER}")

if(PQC_WITH_UPPER STREQUAL "ALL")
    set(PQC_WITH_EFFECTIVE ${PQC_ALL_MODULES})
else()
    set(PQC_WITH_EFFECTIVE ${PQC_WITH_UPPER})
endif()

foreach(module ${PQC_ALL_MODULES})
    list(FIND PQC_WITH_EFFECTIVE ${module} module_index)
    if(module_index EQUAL -1)
        set(PQC_WITH_${module} OFF)
    else()
        set(PQC_WITH_${module} ON)
    endif()
endforeach()

foreach(module ${PQC_WITH_EFFECTIVE})
    list(FIND PQC_ALL_MODULES ${module} module_index)
    if(module_index EQUAL -1)
        string(REPLACE ";" " " PQC_ALL_MODULES_TEXT "${PQC_ALL_MODULES}")
        message(FATAL_ERROR
            "Unknown module PQC_WITH=${module}, choose from: ${PQC_ALL_MODULES_TEXT}")
    endif()
endforeach()
list(REMOVE_DUPLICATES PQC_WITH_EFFECTIVE)

if(NOT PQC_WITH_DILITHIUM AND NOT PQC_WITH_KYBER AND
   NOT PQC_WITH_SNOVA_R2 AND NOT PQC_WITH_SNOVA_R3 AND
   NOT PQC_WITH_MAYO AND NOT PQC_WITH_MQOM AND NOT PQC_WITH_NGCC_DOVE)
    message(FATAL_ERROR "PQC_WITH selects no module")
endif()

# Algorithm parameter sets, kept under their traditional names.
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

set(PQC_MAYO_VARIANT "1" CACHE STRING "MAYO parameter: 1, 2, 3, 5")
set_property(CACHE PQC_MAYO_VARIANT PROPERTY STRINGS 1 2 3 5)
if(NOT PQC_MAYO_VARIANT MATCHES "^(1|2|3|5)$")
    message(FATAL_ERROR "PQC_MAYO_VARIANT must be 1, 2, 3, or 5")
endif()

set(PQC_MQOM_SECURITY "128" CACHE STRING "MQOM security: 128, 192, 256")
set(PQC_MQOM_FIELD "4" CACHE STRING "MQOM log2 field size: 1, 2, 4, 8")
set(PQC_MQOM_ROUNDS "3" CACHE STRING "MQOM rounds: 3 or 5")
set(PQC_MQOM_TRADEOFF "0" CACHE STRING "MQOM tradeoff: 0=fast, 1=short")
set_property(CACHE PQC_MQOM_SECURITY PROPERTY STRINGS 128 192 256)
set_property(CACHE PQC_MQOM_FIELD PROPERTY STRINGS 1 2 4 8)
set_property(CACHE PQC_MQOM_ROUNDS PROPERTY STRINGS 3 5)
set_property(CACHE PQC_MQOM_TRADEOFF PROPERTY STRINGS 0 1)
if(NOT PQC_MQOM_SECURITY MATCHES "^(128|192|256)$" OR
   NOT PQC_MQOM_FIELD MATCHES "^(1|2|4|8)$" OR
   NOT PQC_MQOM_ROUNDS MATCHES "^(3|5)$" OR
   NOT PQC_MQOM_TRADEOFF MATCHES "^[01]$")
    message(FATAL_ERROR "Invalid MQOM parameters")
endif()

set(PQC_NGCC_DOVE_PARAM "128" CACHE STRING "NGCC DOVE security level: 128, 256, 512")
set_property(CACHE PQC_NGCC_DOVE_PARAM PROPERTY STRINGS 128 256 512)
option(PQC_NGCC_DOVE_BIND_SALT
    "Bind the DOVE signature salt into the target (sign-09-1 hardening)" OFF)
if(NOT PQC_NGCC_DOVE_PARAM MATCHES "^(128|256|512)$")
    message(FATAL_ERROR "PQC_NGCC_DOVE_PARAM must be 128, 256, or 512")
endif()

# Values consumed by inc/pqc_conf.h.in.
set(PQC_CONF_WITH_DILITHIUM ${PQC_WITH_DILITHIUM})
set(PQC_CONF_WITH_KYBER ${PQC_WITH_KYBER})
set(PQC_CONF_WITH_SNOVA_R2 ${PQC_WITH_SNOVA_R2})
set(PQC_CONF_WITH_SNOVA_R3 ${PQC_WITH_SNOVA_R3})
set(PQC_CONF_WITH_MAYO ${PQC_WITH_MAYO})
set(PQC_CONF_WITH_MQOM ${PQC_WITH_MQOM})
set(PQC_CONF_WITH_NGCC_DOVE ${PQC_WITH_NGCC_DOVE})
set(PQC_CONF_BUILD_SHARED ${PQC_SHLIB})
set(PQC_CONF_BUILD_STATIC ${PQC_STLIB})
