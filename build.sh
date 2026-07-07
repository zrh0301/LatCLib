#!/bin/bash
# Build script for DilithiumProject
# Usage: ./build.sh [clean|rebuild]

set -e  # 遇到错误立即退出

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"

echo "==> Project dir : ${PROJECT_DIR}"
echo "==> Build dir   : ${BUILD_DIR}"

# 默认行为：先清理再构建（rebuild）
echo "==> Cleaning previous build..."
rm -rf "${BUILD_DIR}"

# 创建 build 目录并进入
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# 运行 CMake 配置
echo "==> Running cmake .."
cmake ..

# 编译（使用所有 CPU 核心）
echo "==> Running make -j$(nproc)"
make -j$(nproc)

echo ""
echo "========================================"
echo "  Build completed successfully!"
echo "  Executables are in: ${BUILD_DIR}"
echo "========================================"
