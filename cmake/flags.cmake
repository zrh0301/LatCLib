# Common compiler and linker flags.

set(PQC_WARNING_FLAGS -Wall -Wextra)
set(PQC_NATIVE_FLAGS "")
set(PQC_STRIP_COMPILE_FLAGS "")
set(PQC_STRIP_LINK_FLAGS "")
set(PQC_SANITIZER_FLAGS "")

if(PQC_NATIVE)
    set(PQC_NATIVE_FLAGS -march=native -mtune=native)
endif()

if(PQC_STRIP)
    set(PQC_STRIP_COMPILE_FLAGS -ffunction-sections -fdata-sections)
    if(APPLE)
        set(PQC_STRIP_LINK_FLAGS -Wl,-dead_strip)
    else()
        set(PQC_STRIP_LINK_FLAGS -Wl,--gc-sections)
    endif()
endif()

if(PQC_SANITIZERS)
    set(PQC_SANITIZER_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer)
endif()

set(PQC_CFLAGS_LIST "")
if(PQC_CFLAGS)
    separate_arguments(PQC_CFLAGS_LIST UNIX_COMMAND "${PQC_CFLAGS}")
endif()

set(PQC_COMMON_COMPILE_FLAGS ${PQC_WARNING_FLAGS} ${PQC_CFLAGS_LIST}
    ${PQC_NATIVE_FLAGS} ${PQC_STRIP_COMPILE_FLAGS} ${PQC_SANITIZER_FLAGS})
set(PQC_COMMON_LINK_FLAGS ${PQC_STRIP_LINK_FLAGS} ${PQC_SANITIZER_FLAGS})

function(pqc_apply_compile_flags target)
    target_compile_options(${target} PRIVATE ${PQC_COMMON_COMPILE_FLAGS})
endfunction()

function(pqc_apply_link_flags target)
    target_link_options(${target} PRIVATE ${PQC_COMMON_LINK_FLAGS})
endfunction()
