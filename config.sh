#!/usr/bin/env bash
#
# config.sh — OpenSSL ./config 风格的交互式配置脚本。
# 逐项菜单选择构建选项, 任何一步回车即接受默认值, 最后调用 CMake 配置。
#
# 用法:
#   ./config.sh [构建目录]     # 构建目录默认 ./build
#
# 构建系统参考 relic-toolkit: 选项为 PQC_ARCH / PQC_WITH / PQC_SHLIB 等。
# 运行结束会打印本次选择等价的手工 CMake 命令。

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR" || exit 1

BUILD_DIR="${1:-build}"

if ! command -v cmake >/dev/null 2>&1; then
    echo "错误: 未找到 cmake 命令 (本库需要 CMake 3.16+)。" >&2
    exit 1
fi

# ---- 交互辅助: 菜单与提示写到 stderr, 选中值写到 stdout (供 $(...) 捕获) ----
# choose <默认值> <提示语> "值:菜单文字" ...
choose() {
    local default="$1" prompt="$2"; shift 2
    local items=("$@") i pick val disp mark
    echo "==== $prompt ====" >&2
    for i in "${!items[@]}"; do
        val="${items[$i]%%:*}"
        disp="${items[$i]#*:}"
        mark=""
        [ "$val" = "$default" ] && mark="  <- 回车默认"
        printf "  [%d] %s%s\n" "$((i + 1))" "$disp" "$mark" >&2
    done
    while :; do
        printf "选择 [1-%d, 回车=%s]: " "${#items[@]}" "$default" >&2
        if ! IFS= read -r pick; then
            echo >&2
            pick=""
        fi
        if [ -z "$pick" ]; then
            echo "$default"
            return
        fi
        case "$pick" in
            *[!0-9]*)
                echo "  无效输入, 请输入 1-${#items[@]} (回车取默认)。" >&2
                ;;
            *)
                if [ "$pick" -ge 1 ] && [ "$pick" -le "${#items[@]}" ]; then
                    echo "${items[$((pick - 1))]%%:*}"
                    return
                fi
                echo "  无效输入, 请输入 1-${#items[@]} (回车取默认)。" >&2
                ;;
        esac
    done
}

# ask_bool <说明> <默认 ON|OFF>
ask_bool() {
    local prompt="$1" default="$2" ans hint
    case "$default" in
        ON) hint="[Y/n]" ;;
        *)  hint="[y/N]" ;;
    esac
    while :; do
        printf "%s %s: " "$prompt" "$hint" >&2
        if ! IFS= read -r ans; then
            echo >&2
            ans=""
        fi
        case "$ans" in
            '')  echo "$default"; return ;;
            y|Y|yes|YES) echo "ON";  return ;;
            n|N|no|NO)   echo "OFF"; return ;;
            *)  echo "  请输入 y 或 n (回车=$default)。" >&2 ;;
        esac
    done
}

echo "本库交互式构建配置 (每项回车接受默认值)"
echo "源码目录: $SCRIPT_DIR"
echo

BUILD_TYPE=$(choose Release "构建类型" \
    "Release:Release  (发布构建, 默认)" \
    "Debug:Debug  (调试, 无优化, 含断言)" \
    "RelWithDebInfo:RelWithDebInfo  (优化 + 调试信息)" \
    "MinSizeRel:MinSizeRel  (最小体积)")

ARCH=$(choose AUTO "目标架构" \
    "AUTO:自动探测主机架构 (默认)" \
    "X64:x86-64" \
    "A64:AArch64")

BACKEND=$(choose AUTO "签名算法后端" \
    "AUTO:跟随架构 (x86-64=AVX2, AArch64=NEON, 默认)" \
    "AVX2:AVX2 实现" \
    "GFNI:AVX2 + GFNI 实现" \
    "NEON:AArch64 NEON 实现")

GF_BACKEND=$(choose AUTO "GF 运算后端" \
    "AUTO:跟随架构 (默认)" \
    "AVX2:AVX2 实现" \
    "GFNI:AVX2 + GFNI 实现" \
    "NEON:AArch64 NEON 实现")

AES_BACKEND=$(choose AUTO "AES 后端" \
    "AUTO:自动探测主机 AES-NI (默认)" \
    "C:纯 C 实现 (任何平台可用)" \
    "NI:强制 AES-NI (需 x86-64)")

echo "==== 选择要构建的算法模块 (y/N, 回车=Y) ====" >&2
MODULES=""
for module in DILITHIUM KYBER SNOVA_R2 SNOVA_R3 MAYO MQOM; do
    enabled=$(ask_bool "构建 ${module}" ON)
    if [ "$enabled" = ON ]; then
        if [ -z "$MODULES" ]; then
            MODULES="$module"
        else
            MODULES="$MODULES;$module"
        fi
    fi
done
if [ -z "$MODULES" ]; then
    echo "错误: 至少选择一个算法模块。" >&2
    exit 1
fi

SHLIB=$(ask_bool "构建共享库" ON)
STLIB=$(ask_bool "构建静态库" ON)
if [ "$SHLIB" = OFF ] && [ "$STLIB" = OFF ]; then
    echo "错误: 共享库与静态库至少构建一种。" >&2
    exit 1
fi

TESTS=$(ask_bool "构建并注册 CTest 正确性测试" ON)
PERF_TESTS=$(ask_bool "把性能测试注册进 CTest (开启后另生成 bench.sh)" OFF)
NATIVE_OPT=$(ask_bool "针对本机 CPU 优化 (-march=native, 产物不可跨机移植)" OFF)
SANITIZERS=$(ask_bool "启用 ASan + UBSan (建议搭配 Debug 构建)" OFF)

DILITHIUM_MODE=$(choose 0 "Dilithium 参数集" \
    "0:全部参数集 (默认)" \
    "2:Dilithium2" \
    "3:Dilithium3" \
    "5:Dilithium5")

KYBER_MODE=$(choose 0 "Kyber 参数集" \
    "0:全部参数集 (默认)" \
    "2:Kyber2 (PKE)" \
    "3:Kyber3 (PKE)" \
    "5:Kyber5 (PKE)")

MAYO_VARIANT=$(choose 1 "MAYO 参数集" \
    "1:MAYO-1 (默认)" \
    "2:MAYO-2" \
    "3:MAYO-3" \
    "5:MAYO-5")

if [ "$SANITIZERS" = ON ] && [ "$BUILD_TYPE" != Debug ]; then
    echo
    echo "提示: sanitizer 通常配合 Debug 构建; Release 下结果不代表常规性能。" >&2
fi

CMAKE_ARGS=(-S . -B "$BUILD_DIR"
             -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
             -DPQC_ARCH="$ARCH"
             -DPQC_BACKEND="$BACKEND"
             -DPQC_GF_BACKEND="$GF_BACKEND"
             -DPQC_AES="$AES_BACKEND"
             -DPQC_WITH="$MODULES"
             -DPQC_SHLIB="$SHLIB"
             -DPQC_STLIB="$STLIB"
             -DPQC_TESTS="$TESTS"
             -DPQC_PERF_TESTS="$PERF_TESTS"
             -DPQC_NATIVE="$NATIVE_OPT"
             -DPQC_SANITIZERS="$SANITIZERS"
             -DPQC_DILITHIUM_MODE="$DILITHIUM_MODE"
             -DPQC_KYBER_MODE="$KYBER_MODE"
             -DPQC_MAYO_VARIANT="$MAYO_VARIANT")

echo
echo "==================== 选择汇总 ===================="
printf '构建目录     : %s\n' "$BUILD_DIR"
printf '构建类型     : %s\n' "$BUILD_TYPE"
printf '架构         : %s\n' "$ARCH"
printf '签名后端     : %s\n' "$BACKEND"
printf 'GF 后端      : %s\n' "$GF_BACKEND"
printf 'AES 后端     : %s\n' "$AES_BACKEND"
printf '算法模块     : %s\n' "$MODULES"
printf '共享/静态库  : %s / %s\n' "$SHLIB" "$STLIB"
printf '正确性测试   : %s\n' "$TESTS"
printf '性能测试注册 : %s\n' "$PERF_TESTS"
printf '本机优化     : %s\n' "$NATIVE_OPT"
printf 'ASan+UBSan   : %s\n' "$SANITIZERS"
printf 'Dilithium    : %s\n' "$DILITHIUM_MODE"
printf 'Kyber        : %s\n' "$KYBER_MODE"
printf 'MAYO         : %s\n' "$MAYO_VARIANT"
echo "=================================================="
echo "执行: cmake ${CMAKE_ARGS[*]}"
echo

cmake "${CMAKE_ARGS[@]}" || exit $?

echo
echo "配置完成。接下来执行:"
echo "  cmake --build $BUILD_DIR -j"
if [ "$TESTS" = ON ]; then
    echo "  ctest --test-dir $BUILD_DIR --output-on-failure -LE performance"
fi
if [ "$PERF_TESTS" = ON ]; then
    echo "  ctest --test-dir $BUILD_DIR -L performance    # 性能基准"
fi
echo "  cmake --install $BUILD_DIR --prefix <安装前缀>   # 可选安装"
