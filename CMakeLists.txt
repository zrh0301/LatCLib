cmake_minimum_required(VERSION 3.10)

# 项目名称
project(DilithiumProject C)

# ==========================================
# 1. 基础配置
# ==========================================

# 设置 C 标准
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 【修改点】设置输出目录为当前构建目录 (即 build 文件夹)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# 开启 AddressSanitizer (ASan) 和调试信息 (-g)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -g")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")

# 添加头文件路径

# 宏定义
add_compile_definitions(USE_DILITHIUM_NTT USE_KYBER_NTT)

# ==========================================
# 2. 编译选项 (CFLAGS)
# ==========================================

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    # Linux GCC/Clang 或 Windows MinGW
    add_compile_options(
        -O3 
        -march=native 
        -mtune=native 
        -ffast-math 
        -funroll-loops 
        -g 
        -Wall 
        -Wextra
    )
elseif(MSVC)
    # Windows Visual Studio
    add_compile_options(/O2 /W4 /arch:AVX2 /Zi)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

# ==========================================
# 3. 核心库构建
# ==========================================

file(GLOB COMMON_SRCS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
)
file(GLOB_RECURSE ALG_SRCS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/alg/*/*.c"
)

add_library(core STATIC ${COMMON_SRCS} ${ALG_SRCS})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc")

file(GLOB ALG_DIRS LIST_DIRECTORIES true "${CMAKE_CURRENT_SOURCE_DIR}/alg/*")
foreach(ALG_DIR ${ALG_DIRS})
    if(IS_DIRECTORY "${ALG_DIR}")
        target_include_directories(core PUBLIC "${ALG_DIR}")
    endif()
endforeach()

target_link_libraries(core m)

# ==========================================
# 4. 生成可执行文件
# ==========================================

macro(add_crypto_test TEST_NAME SOURCE_FILE)
    add_executable(${TEST_NAME} ${SOURCE_FILE})
    target_include_directories(${TEST_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Test")
    target_link_libraries(${TEST_NAME} core)
endmacro()

# 添加测试目标 (已移除 test_sha3)
add_crypto_test(test_poly     Test/poly/test_poly.c)
add_crypto_test(test_ntt      Test/ntt/test_ntt.c)
add_crypto_test(test_sampling Test/sampling/test_sampling.c)
add_crypto_test(test_pack     Test/pack/test_pack.c)
add_crypto_test(test_hash     Test/hash/test_hash.c)
add_crypto_test(test_sig      Test/sig/test_sig.c)
add_crypto_test(test_enc      Test/enc/test_enc.c)
add_crypto_test(print_zetas   Test/tools/print_zeta.c)

target_sources(test_sig PRIVATE Test/sig/cases.c)
target_sources(test_enc PRIVATE Test/enc/cases.c)
target_sources(test_hash PRIVATE Test/hash/cases.c)
target_sources(test_pack PRIVATE Test/pack/cases.c Test/pack/dilithium.c Test/pack/kyber.c)

# ==========================================
# 5. 自定义 Test 命令
# ==========================================
# 【修改点】使用 generator expression 自动定位文件位置
# 无论是在 build/ 下还是 build/Release/ 下都能找到
add_custom_target(run_test
    COMMAND $<TARGET_FILE:test_sig>
    COMMENT "Running Signature Tests..."
    DEPENDS test_sig
)
