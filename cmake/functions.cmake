# Helpers that build one OBJECT library and derive the shared and static
# variants from the same objects, so that sources compile only once.

function(pqc_add_variants base)
    cmake_parse_arguments(ARG ""
        ""
        "SOURCES;INCLUDES;PUBLIC_INCLUDES;DEFINITIONS;PUBLIC_DEFINITIONS;COMPILE_OPTIONS;DEPENDS;LINK_LIBRARIES;EXTRA_OBJECTS"
        ${ARGN})

    if(NOT ARG_SOURCES)
        message(FATAL_ERROR "pqc_add_variants(${base}) requires SOURCES")
    endif()

    add_library(${base}_obj OBJECT ${ARG_SOURCES})
    set(pqc_extra_objects "")
    foreach(extra_object IN LISTS ARG_EXTRA_OBJECTS)
        list(APPEND pqc_extra_objects $<TARGET_OBJECTS:${extra_object}>)
    endforeach()
    set_target_properties(${base}_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
    if(ARG_PUBLIC_INCLUDES)
        target_include_directories(${base}_obj PUBLIC ${ARG_PUBLIC_INCLUDES})
    endif()
    if(ARG_INCLUDES)
        target_include_directories(${base}_obj PRIVATE ${ARG_INCLUDES})
    endif()
    if(ARG_DEFINITIONS)
        target_compile_definitions(${base}_obj PRIVATE ${ARG_DEFINITIONS})
    endif()
    if(ARG_PUBLIC_DEFINITIONS)
        target_compile_definitions(${base}_obj PUBLIC ${ARG_PUBLIC_DEFINITIONS})
    endif()
    if(ARG_COMPILE_OPTIONS)
        target_compile_options(${base}_obj PRIVATE ${ARG_COMPILE_OPTIONS})
    endif()
    pqc_apply_compile_flags(${base}_obj)
    foreach(dep IN LISTS ARG_DEPENDS)
        target_link_libraries(${base}_obj PRIVATE ${dep}_obj)
    endforeach()

    set(pqc_variant_targets "")

    if(PQC_SHLIB)
        add_library(${base} SHARED
            $<TARGET_OBJECTS:${base}_obj> ${pqc_extra_objects})
        set_target_properties(${base} PROPERTIES
            SOVERSION ${PROJECT_VERSION_MAJOR}
            VERSION ${PROJECT_VERSION})
        if(ARG_PUBLIC_INCLUDES)
            target_include_directories(${base} PUBLIC ${ARG_PUBLIC_INCLUDES})
        endif()
        if(ARG_PUBLIC_DEFINITIONS)
            target_compile_definitions(${base} PUBLIC ${ARG_PUBLIC_DEFINITIONS})
        endif()
        target_link_libraries(${base} PUBLIC ${ARG_DEPENDS} ${ARG_LINK_LIBRARIES})
        if(UNIX AND NOT APPLE)
            set_target_properties(${base} PROPERTIES INSTALL_RPATH "$ORIGIN")
        endif()
        pqc_apply_link_flags(${base})
        list(APPEND pqc_variant_targets ${base})
    endif()

    if(PQC_STLIB)
        add_library(${base}_static STATIC
            $<TARGET_OBJECTS:${base}_obj> ${pqc_extra_objects})
        set_target_properties(${base}_static PROPERTIES OUTPUT_NAME ${base})
        if(WIN32 AND PQC_SHLIB)
            # Keep the import library of the shared variant apart.
            set_target_properties(${base}_static PROPERTIES
                ARCHIVE_OUTPUT_NAME ${base}_static)
        endif()
        if(ARG_PUBLIC_INCLUDES)
            target_include_directories(${base}_static PUBLIC ${ARG_PUBLIC_INCLUDES})
        endif()
        if(ARG_PUBLIC_DEFINITIONS)
            target_compile_definitions(${base}_static PUBLIC ${ARG_PUBLIC_DEFINITIONS})
        endif()
        set(pqc_static_deps "")
        foreach(dep IN LISTS ARG_DEPENDS)
            list(APPEND pqc_static_deps ${dep}_static)
        endforeach()
        target_link_libraries(${base}_static PUBLIC
            ${pqc_static_deps} ${ARG_LINK_LIBRARIES})
        pqc_apply_link_flags(${base}_static)
        list(APPEND pqc_variant_targets ${base}_static)
    endif()

    set_property(GLOBAL APPEND PROPERTY PQC_LIBRARY_TARGETS ${pqc_variant_targets})

    if(PQC_STLIB)
        set(PQC_LINK_${base} ${base}_static PARENT_SCOPE)
    else()
        set(PQC_LINK_${base} ${base} PARENT_SCOPE)
    endif()
endfunction()

# Public include directories for targets installed under include/pqc.
set(PQC_PUBLIC_INCLUDE_DIRS
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>"
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/inc>"
    "$<INSTALL_INTERFACE:include/pqc>")
