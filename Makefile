# ==========================================
# My_PQC Full Suite Makefile
# ==========================================

CC = gcc

# 编译选项
CFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops -Wall -Wextra -I./inc -I./test/inc

# Windows 需要链接 advapi32 库
LDFLAGS = -ladvapi32

# 输出目录
OUTPUT_DIR = output

# 测试程序目标
TARGET_POLY = test_poly.exe
TARGET_NTT = test_ntt.exe

# 核心库源文件
CORE_SRCS = src/poly.c src/poly_ntt.c src/rand.c
CORE_OBJS = $(addprefix $(OUTPUT_DIR)/, $(notdir $(CORE_SRCS:.c=.o)))

# 测试源文件（新位置：test/src/）
TEST_POLY_SRCS = test/src/test_poly.c
TEST_NTT_SRCS = test/src/test_ntt.c
TEST_POLY_OBJS = $(addprefix $(OUTPUT_DIR)/, $(notdir $(TEST_POLY_SRCS:.c=.o)))
TEST_NTT_OBJS = $(addprefix $(OUTPUT_DIR)/, $(notdir $(TEST_NTT_SRCS:.c=.o)))

# 头文件依赖
DEPS = inc/poly.h inc/types.h inc/rand.h test/inc/test_common.h

# 默认目标：编译所有测试
all: $(OUTPUT_DIR) $(TARGET_POLY) $(TARGET_NTT)

# 创建输出目录
$(OUTPUT_DIR):
	@if not exist $(OUTPUT_DIR) mkdir $(OUTPUT_DIR)

# 链接 test_poly
$(TARGET_POLY): $(CORE_OBJS) $(TEST_POLY_OBJS)
	@echo "Linking $(TARGET_POLY)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_ntt
$(TARGET_NTT): $(CORE_OBJS) $(TEST_NTT_OBJS)
	@echo "Linking $(TARGET_NTT)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 编译规则 - src 目录
$(OUTPUT_DIR)/%.o: src/%.c $(DEPS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 编译规则 - test/src 目录
$(OUTPUT_DIR)/%.o: test/src/%.c $(DEPS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 运行所有测试
test: all
	@echo.
	@echo =============================================
	@echo    Running Polynomial Tests
	@echo =============================================
	@echo.
	./$(TARGET_POLY)
	@echo.
	@echo.
	@echo =============================================
	@echo    Running NTT Tests
	@echo =============================================
	@echo.
	./$(TARGET_NTT)

# 清理
clean:
	-del /Q $(OUTPUT_DIR)\*.o $(TARGET_POLY) $(TARGET_NTT) 2>nul
	-rmdir $(OUTPUT_DIR) 2>nul

.PHONY: all test clean