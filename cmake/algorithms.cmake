# Core library and signature algorithm modules.

# ---- Core: shared operators, hash/symmetric layers, Dilithium and Kyber ----

set(pqc_core_sources
    src/poly/fp32.c
    src/poly/fp16.c
    src/poly/poly16.c
    src/poly/poly16_pack.c
    src/poly/ntt16.c
    src/poly/ring.c
    src/poly/poly.c
    src/poly/poly_ntt.c
    src/poly/polyvec.c
    src/poly/ntt_table.c
    src/poly/ntt_compat.c
    src/poly/ring4171.c
    src/poly/poly_sample.c
    src/util/rand.c
    src/util/secure_clear.c
    src/util/ct.c
    src/util/combinatorics.c
    src/tree/merkle.c
    src/sym/keccak.c
    src/sym/aes.c
    src/sym/sm3.c
    src/ngcc/auxfunc.c
    src/ngcc/drng.c
    src/core/pqc_conf.c)

set(pqc_core_private_includes
    "${PROJECT_SOURCE_DIR}/alg/dilithium"
    "${PROJECT_SOURCE_DIR}/alg/kyber")
set(pqc_core_public_includes ${PQC_PUBLIC_INCLUDE_DIRS}
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/alg/dilithium>"
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/alg/kyber>")

if(PQC_WITH_DILITHIUM)
    list(APPEND pqc_core_sources
        alg/dilithium/dili.c
        alg/dilithium/dili_sample.c
        alg/dilithium/pack.c)
endif()

if(PQC_WITH_KYBER)
    list(APPEND pqc_core_sources
        alg/kyber/kyber.c
        alg/kyber/kyber_sample.c
        alg/kyber/kyber_pack.c)
endif()

if(PQC_AES_EFFECTIVE STREQUAL "NI")
    set_source_files_properties(src/sym/aes.c PROPERTIES
        COMPILE_OPTIONS "-maes;-mavx2")
endif()

set(pqc_core_link_libraries "")
if(WIN32)
    list(APPEND pqc_core_link_libraries bcrypt)
else()
    list(APPEND pqc_core_link_libraries m)
endif()

pqc_add_variants(pqc
    SOURCES ${pqc_core_sources}
    INCLUDES ${pqc_core_private_includes}
    PUBLIC_INCLUDES ${pqc_core_public_includes}
    DEPENDS pqc_gf
    LINK_LIBRARIES ${pqc_core_link_libraries})

# ---- SNOVA R3 ----

if(PQC_WITH_SNOVA_R3)
    set(pqc_r3_opt 10)
    if(PQC_BACKEND_EFFECTIVE MATCHES "^(AVX2|GFNI)$")
        set(pqc_r3_dir "${PROJECT_SOURCE_DIR}/alg/snova_r3/optimized/native")
        set(pqc_r3_core "${pqc_r3_dir}/snova_avx2_16.c")
        set(pqc_r3_flags -mavx2)
        set(pqc_r3_opt 20)
        if(PQC_BACKEND_EFFECTIVE STREQUAL "GFNI")
            list(APPEND pqc_r3_flags -mgfni)
        endif()
    elseif(PQC_BACKEND_EFFECTIVE STREQUAL "NEON")
        set(pqc_r3_dir "${PROJECT_SOURCE_DIR}/alg/snova_r3/optimized/armv8")
        set(pqc_r3_core "${pqc_r3_dir}/snova_opt_16.c")
        set(pqc_r3_flags -march=armv8-a+crypto)
    else()
        message(FATAL_ERROR "Unknown PQC_BACKEND=${PQC_BACKEND_EFFECTIVE}")
    endif()

    pqc_add_variants(pqc_snova_r3
        SOURCES
            ${pqc_r3_core}
            "${pqc_r3_dir}/symmetric.c"
            "${pqc_r3_dir}/sign.c"
            alg/snova_r3/adapter.c
        INCLUDES "${pqc_r3_dir}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_r3_dir}>"
        PUBLIC_DEFINITIONS SNOVA_OPT=${pqc_r3_opt}
        COMPILE_OPTIONS ${pqc_r3_flags}
        DEPENDS pqc)
endif()

# ---- SNOVA R2 ----

if(PQC_WITH_SNOVA_R2)
    if(PQC_SNOVA_R2_BACKEND_EFFECTIVE STREQUAL "AVX2")
        set(pqc_r2_dir "${PROJECT_SOURCE_DIR}/alg/snova_r2/optimized/x86")
        set(pqc_r2_flags -mavx2)
    elseif(PQC_SNOVA_R2_BACKEND_EFFECTIVE STREQUAL "NEON")
        set(pqc_r2_dir "${PROJECT_SOURCE_DIR}/alg/snova_r2/optimized/armv8")
        set(pqc_r2_flags -march=armv8-a+crypto)
    else()
        message(FATAL_ERROR
            "Unknown SNOVA R2 backend ${PQC_SNOVA_R2_BACKEND_EFFECTIVE}")
    endif()

    set(pqc_r2_sources
        "${pqc_r2_dir}/snova.c"
        "${pqc_r2_dir}/snova_vecops.c"
        "${pqc_r2_dir}/snova_transform.c"
        "${pqc_r2_dir}/snova_public_layout.c"
        "${pqc_r2_dir}/aes_armv8_crypto.c"
        "${pqc_r2_dir}/symmetric.c"
        alg/snova_r2/adapter.c
        alg/snova_r2/randombytes.c)

    pqc_add_variants(pqc_snova_r2
        SOURCES ${pqc_r2_sources}
        INCLUDES "${pqc_r2_dir}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
        DEFINITIONS
            SNOVA_v=24 SNOVA_o=5 SNOVA_q=16 SNOVA_l=4 AESCTR
            SNOVA_OPT=_avx2 SNOVA_OFFICIAL_API_ULL SNOVA_BENCH_OWN_RANDOMBYTES
            SNOVA_VERIFY_DIRECT_MTK2=1 SNOVA_VERIFY_VTL=0
            SNOVA_VERIFY_DIRECT_OFFICIAL=0 SNOVA_VTL_MULTTAB_LOOKUP=1
            randombytes=pqc_snova_r2_randombytes
            crypto_sign=pqc_r2_crypto_sign
            crypto_sign_keypair=pqc_r2_crypto_sign_keypair
            crypto_sign_open=pqc_r2_crypto_sign_open
        COMPILE_OPTIONS ${pqc_r2_flags}
            -include "${PROJECT_SOURCE_DIR}/alg/snova_r2/namespace.h"
        DEPENDS pqc)
endif()

# ---- MAYO ----

if(PQC_WITH_MAYO)
    set(pqc_mayo_dir "${PROJECT_SOURCE_DIR}/alg/mayo/common")
    set(pqc_mayo_sources
        "${pqc_mayo_dir}/mayo.c"
        "${pqc_mayo_dir}/params.c"
        "${pqc_mayo_dir}/arithmetic.c"
        "${pqc_mayo_dir}/common/aes_c.c"
        "${pqc_mayo_dir}/common/aes128ctr.c"
        "${pqc_mayo_dir}/common/aes_neon.c"
        "${pqc_mayo_dir}/common/fips202.c"
        "${pqc_mayo_dir}/common/mem.c"
        alg/mayo/adapter.c
        alg/mayo/randombytes.c)

    set(pqc_mayo_includes
        "${PROJECT_SOURCE_DIR}/alg/mayo/include"
        "${pqc_mayo_dir}"
        "${pqc_mayo_dir}/common"
        "${pqc_mayo_dir}/generic")
    set(pqc_mayo_definitions
        MAYO_VARIANT=MAYO_${PQC_MAYO_VARIANT}
        randombytes=pqc_mayo_randombytes)
    set(pqc_mayo_flags "")

    if(PQC_BACKEND_EFFECTIVE MATCHES "^(AVX2|GFNI)$")
        list(APPEND pqc_mayo_includes "${pqc_mayo_dir}/AVX2")
        list(APPEND pqc_mayo_definitions MAYO_AVX MAYO_BUILD_TYPE_AVX2)
        list(APPEND pqc_mayo_flags -mavx2)
        if(PQC_BACKEND_EFFECTIVE STREQUAL "GFNI")
            list(APPEND pqc_mayo_flags -mgfni)
            list(APPEND pqc_mayo_definitions
                MAYO_GF16_MATVEC MAYO_GF16_GFNI MAYO_GF16_GFNI_MULTAB
                MAYO_GF16_GFNI_CROSS_ROW MAYO_FUSED_SIGN_M MAYO_PACKED_SIGN_P
                MAYO_INPLACE_UNPACK MAYO_VECTOR_UNPACK)
        endif()
    elseif(PQC_BACKEND_EFFECTIVE STREQUAL "NEON")
        list(APPEND pqc_mayo_includes "${pqc_mayo_dir}/neon")
        list(APPEND pqc_mayo_definitions MAYO_NEON MAYO_BUILD_TYPE_NEON)
        list(APPEND pqc_mayo_flags -march=armv8-a+crypto)
    else()
        message(FATAL_ERROR "Unknown PQC_BACKEND=${PQC_BACKEND_EFFECTIVE}")
    endif()

    pqc_add_variants(pqc_mayo
        SOURCES ${pqc_mayo_sources}
        INCLUDES ${pqc_mayo_includes}
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
        DEFINITIONS ${pqc_mayo_definitions}
        COMPILE_OPTIONS ${pqc_mayo_flags}
            -include "${PROJECT_SOURCE_DIR}/alg/mayo/namespace.h"
        DEPENDS pqc)
endif()

# ---- MQOM ----

if(PQC_WITH_MQOM)
    if(PQC_BACKEND_EFFECTIVE MATCHES "^(AVX2|GFNI)$")
        set(pqc_mqom_dir "${PROJECT_SOURCE_DIR}/alg/mqom/optimized/x86")
        set(pqc_mqom_flags -mavx2)
        if(PQC_BACKEND_EFFECTIVE STREQUAL "GFNI")
            list(APPEND pqc_mqom_flags -mgfni)
        endif()
    elseif(PQC_BACKEND_EFFECTIVE STREQUAL "NEON")
        set(pqc_mqom_dir "${PROJECT_SOURCE_DIR}/alg/mqom/optimized/armv8")
        set(pqc_mqom_flags -march=armv8-a+crypto)
    else()
        message(FATAL_ERROR "Unknown PQC_BACKEND=${PQC_BACKEND_EFFECTIVE}")
    endif()

    set(pqc_mqom_sources "")
    foreach(file xof prg expand_mq keygen sign sign_memopt crypto_sign)
        list(APPEND pqc_mqom_sources "${pqc_mqom_dir}/${file}.c")
    endforeach()
    foreach(dir rijndael ggm_tree blc piop)
        file(GLOB pqc_mqom_part CONFIGURE_DEPENDS "${pqc_mqom_dir}/${dir}/*.c")
        list(APPEND pqc_mqom_sources ${pqc_mqom_part})
    endforeach()
    foreach(file KeccakHash KeccakSponge KeccakSpongetimes4 KeccakHashtimes4)
        list(APPEND pqc_mqom_sources "${pqc_mqom_dir}/sha3/${file}.c")
    endforeach()
    list(APPEND pqc_mqom_sources
        "${pqc_mqom_dir}/sha3/opt64/KeccakP-1600-opt64.c"
        "${pqc_mqom_dir}/sha3/opt64/KeccakP-1600-times4-on1.c"
        alg/mqom/adapter.c
        alg/mqom/randombytes.c)

    set(pqc_mqom_definitions
        MQOM2_PARAM_SECURITY=${PQC_MQOM_SECURITY}
        MQOM2_PARAM_BASE_FIELD=${PQC_MQOM_FIELD}
        MQOM2_PARAM_NBROUNDS=${PQC_MQOM_ROUNDS}
        MQOM2_PARAM_TRADEOFF=${PQC_MQOM_TRADEOFF}
        randombytes=pqc_mqom_randombytes
        APPLY_NAMESPACE=latclib_mqom_
        APPLY_PUBLIC_API_NAMESPACE=latclib_mqom_
        RIJNDAEL_CONSTANT_TIME_REF
        NO_MISALIGNED_ACCESSES)
    if(PQC_BACKEND_EFFECTIVE MATCHES "^(AVX2|GFNI)$")
        list(APPEND pqc_mqom_definitions FIELDS_AVX2)
    elseif(PQC_BACKEND_EFFECTIVE STREQUAL "NEON")
        list(APPEND pqc_mqom_definitions FIELDS_NEON)
    endif()

    pqc_add_variants(pqc_mqom
        SOURCES ${pqc_mqom_sources}
        INCLUDES
            "${pqc_mqom_dir}"
            "${pqc_mqom_dir}/sha3"
            "${pqc_mqom_dir}/sha3/opt64"
            "${pqc_mqom_dir}/rijndael"
            "${pqc_mqom_dir}/ggm_tree"
            "${pqc_mqom_dir}/blc"
            "${pqc_mqom_dir}/piop"
            "${pqc_mqom_dir}/fields"
            "${pqc_mqom_dir}/fields_bitsliced"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
        DEFINITIONS ${pqc_mqom_definitions}
        COMPILE_OPTIONS ${pqc_mqom_flags}
        DEPENDS pqc)
endif()

# ---- NGCC: DOVE (round-1 signature candidate, UOV over GF(256)) ----

if(PQC_WITH_NGCC_DOVE)
    set(pqc_dove_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_dove/reference")

    set(pqc_dove_definitions DOVE_PARAM_SET=DOVE${PQC_NGCC_DOVE_PARAM})
    if(PQC_NGCC_DOVE_BIND_SALT)
        list(APPEND pqc_dove_definitions DOVE_BIND_SALT=1)
    endif()

    pqc_add_variants(pqc_ngcc_dove
        SOURCES
            "${pqc_dove_ref}/SIG_AlgorithmInstance.c"
            alg/ngcc_dove/adapter.c
        INCLUDES "${pqc_dove_ref}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_dove_ref}>"
        PUBLIC_DEFINITIONS ${pqc_dove_definitions}
        COMPILE_OPTIONS -Wno-unused-parameter -Wno-unused-variable
        DEPENDS pqc)
endif()

# ---- NGCC: VDOO (round-1 signature candidate, layered Oil-and-Vinegar) ----

if(PQC_WITH_NGCC_VDOO)
    set(pqc_vdoo_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_vdoo/reference")

    set(pqc_vdoo_namespace
        sig_keygen=pqc_ngcc_vdoo_keygen
        sig_sign=pqc_ngcc_vdoo_sign
        sig_verify=pqc_ngcc_vdoo_verify
        sig_get_pk_len_bytes=pqc_ngcc_vdoo_pk_len
        sig_get_sk_len_bytes=pqc_ngcc_vdoo_sk_len
        sig_get_sn_len_bytes=pqc_ngcc_vdoo_sig_len)
    pqc_add_variants(pqc_ngcc_vdoo
        SOURCES
            "${pqc_vdoo_ref}/SIG_AlgorithmInstance.c"
            "${pqc_vdoo_ref}/api.c"
            "${pqc_vdoo_ref}/vdoo_keypair.c"
            "${pqc_vdoo_ref}/vdoo_sign.c"
            "${pqc_vdoo_ref}/vdoo_verif.c"
            alg/ngcc_vdoo/adapter.c
        INCLUDES
            "${pqc_vdoo_ref}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_vdoo_ref}>"
        DEFINITIONS
            ${pqc_vdoo_namespace}
            crypto_sign=pqc_vdoo_crypto_sign
            crypto_sign_keypair=pqc_vdoo_crypto_sign_keypair
            crypto_sign_open=pqc_vdoo_crypto_sign_open
        PUBLIC_DEFINITIONS VDOO_${PQC_NGCC_VDOO_PARAM}
        COMPILE_OPTIONS -Wno-unused-parameter -Wno-unused-variable
        DEPENDS pqc)
endif()

# ---- NGCC: FlexTree (round-1 signature candidate, PORS+FP and WOTS) ----

if(PQC_WITH_NGCC_FLEXTREE)
    set(pqc_flextree_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_flextree/reference")

    set(pqc_flextree_definitions
        crypto_sign=pqc_flextree_crypto_sign
        crypto_sign_keypair=pqc_flextree_crypto_sign_keypair
        crypto_sign_open=pqc_flextree_crypto_sign_open
        drng_algorithm=pqc_ngcc_flextree_drng
        sig_keygen=pqc_ngcc_flextree_keygen
        sig_sign=pqc_ngcc_flextree_sign
        sig_verify=pqc_ngcc_flextree_verify
        sig_get_pk_len_bytes=pqc_ngcc_flextree_pk_len
        sig_get_sk_len_bytes=pqc_ngcc_flextree_sk_len
        sig_get_sn_len_bytes=pqc_ngcc_flextree_sig_len)
    if(PQC_NGCC_FLEXTREE_STRICT)
        list(APPEND pqc_flextree_definitions PQC_NGCC_FLEXTREE_STRICT)
    endif()

    pqc_add_variants(pqc_ngcc_flextree
        SOURCES
            "${pqc_flextree_ref}/address.c"
            "${pqc_flextree_ref}/hash.c"
            "${pqc_flextree_ref}/merkle.c"
            "${pqc_flextree_ref}/pors_fp.c"
            "${pqc_flextree_ref}/sign.c"
            "${pqc_flextree_ref}/wots.c"
            alg/ngcc_flextree/adapter.c
        INCLUDES
            "${pqc_flextree_ref}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_flextree_ref}>"
        DEFINITIONS ${pqc_flextree_definitions}
        COMPILE_OPTIONS -Wno-unused-parameter -Wno-unused-variable
        DEPENDS pqc)
endif()

# ---- NGCC: Aigis-Enc+ (round-1 KEM candidate, module-LWE) ----

if(PQC_WITH_NGCC_AIGIS)
    set(pqc_aigis_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_aigis/reference")
    set(pqc_aigis_namespace "${PROJECT_SOURCE_DIR}/alg/ngcc_aigis/namespace.h")

    pqc_add_variants(pqc_ngcc_aigis
        SOURCES
            "${pqc_aigis_ref}/kem.c"
            "${pqc_aigis_ref}/owcpa.c"
            alg/ngcc_aigis/adapter.c
        INCLUDES
            "${pqc_aigis_ref}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_aigis_ref}>"
        DEFINITIONS PARAMS=1
        COMPILE_OPTIONS -include "${pqc_aigis_namespace}"
        DEPENDS pqc)
endif()

# ---- NGCC: Aigis-Sig+ (round-1 signature candidate, module-LWE/SIS) ----

if(PQC_WITH_NGCC_AIGIS_SIG)
    set(pqc_aigis_sig_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_aigis_sig/reference")
    set(pqc_aigis_sig_namespace
        "${PROJECT_SOURCE_DIR}/alg/ngcc_aigis_sig/namespace.h")

    set(pqc_aigis_sig_definitions PARAMS=1)
    if(PQC_NGCC_AIGIS_SIG_STRICT)
        list(APPEND pqc_aigis_sig_definitions PQC_NGCC_AIGIS_SIG_STRICT)
    endif()

    pqc_add_variants(pqc_ngcc_aigis_sig
        SOURCES
            "${pqc_aigis_sig_ref}/sign.c"
            "${pqc_aigis_sig_ref}/polyvec.c"
            "${pqc_aigis_sig_ref}/packing.c"
            "${pqc_aigis_sig_ref}/poly.c"
            "${pqc_aigis_sig_ref}/rounding.c"
            "${pqc_aigis_sig_ref}/hashkdf.c"
            "${pqc_aigis_sig_ref}/pspm.c"
            "${pqc_aigis_sig_ref}/SIG_AlgorithmInstance.c"
            alg/ngcc_aigis_sig/adapter.c
        INCLUDES
            "${pqc_aigis_sig_ref}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_aigis_sig_ref}>"
        DEFINITIONS ${pqc_aigis_sig_definitions}
        COMPILE_OPTIONS -Wno-unused-parameter -Wno-unused-variable
            -include "${pqc_aigis_sig_namespace}"
        DEPENDS pqc)
endif()

# ---- NGCC: Scloud+ (round-1 KEM candidate, unstructured LWE) ----

if(PQC_WITH_NGCC_SCLOUD)
    set(pqc_scloud_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_scloud/reference")
    set(pqc_scloud_namespace
        "${PROJECT_SOURCE_DIR}/alg/ngcc_scloud/namespace.h")

    # 提交默认实例：128 位安全级 + AES 族。矩阵 A 用 AES-128-CTR 展开；
    # F/G/H/K 使用带标签的 SHAKE256/SHA3-512。
    set(pqc_scloud_definitions
        SCLOUDPLUS_FAMILY_AES
        SCLOUDPLUS_TIER_REFERENCE
        SCLOUDPLUS_REF_FAMILY_AES)

    pqc_add_variants(pqc_ngcc_scloud
        SOURCES
            "${pqc_scloud_ref}/common/kem.c"
            "${pqc_scloud_ref}/common/pke.c"
            "${pqc_scloud_ref}/common/sample.c"
            "${pqc_scloud_ref}/common/encode.c"
            "${pqc_scloud_ref}/common/hash_aes_shake.c"
            "${pqc_scloud_ref}/common/util.c"
            "${pqc_scloud_ref}/ref/matrix_reference.c"
            "${pqc_scloud_ref}/ref/pack_reference.c"
            alg/ngcc_scloud/aes_shim.c
            alg/ngcc_scloud/random.c
            alg/ngcc_scloud/adapter.c
        INCLUDES
            "${pqc_scloud_ref}"
            "${pqc_scloud_ref}/include"
            "${pqc_scloud_ref}/common"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_scloud_ref}>"
            "$<BUILD_INTERFACE:${pqc_scloud_ref}/include>"
        DEFINITIONS ${pqc_scloud_definitions}
        COMPILE_OPTIONS -Wno-unused-parameter -Wno-unused-variable
            -include "${pqc_scloud_namespace}"
        DEPENDS pqc)
endif()

# ---- NGCC: TSUOV (round-1 signature candidate, multivariate over GF(31)) ----

if(PQC_WITH_NGCC_TSUOV)
    set(pqc_tsuov_ref "${PROJECT_SOURCE_DIR}/alg/ngcc_tsuov/reference")
    set(pqc_tsuov_namespace
        "${PROJECT_SOURCE_DIR}/alg/ngcc_tsuov/namespace.h")

    pqc_add_variants(pqc_ngcc_tsuov
        SOURCES
            "${pqc_tsuov_ref}/tsuov_core.c"
            "${pqc_tsuov_ref}/mgf.c"
            "${pqc_tsuov_ref}/SIG_AlgorithmInstance.c"
            alg/ngcc_tsuov/adapter.c
        INCLUDES
            "${pqc_tsuov_ref}"
        PUBLIC_INCLUDES ${PQC_PUBLIC_INCLUDE_DIRS}
            "$<BUILD_INTERFACE:${pqc_tsuov_ref}>"
        DEFINITIONS TSUOV_VARIANT=${PQC_NGCC_TSUOV_PARAM}
        COMPILE_OPTIONS -Wno-unused-parameter -Wno-unused-variable
            -include "${pqc_tsuov_namespace}"
        DEPENDS pqc)
endif()

# ---- Compatibility target: pqc_snova points at SNOVA R3 ----

if(PQC_WITH_SNOVA_R3)
    add_library(pqc_snova INTERFACE)
    target_link_libraries(pqc_snova INTERFACE ${PQC_LINK_pqc_snova_r3})
    set_property(GLOBAL APPEND PROPERTY PQC_LIBRARY_TARGETS pqc_snova)
endif()
