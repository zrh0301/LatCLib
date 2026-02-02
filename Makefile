CC = gcc
# 编译选项：
# -O3: 最高级优化
# -march=native: 针对本机CPU指令集优化 (对NTT性能至关重要)
# -Wall -Wextra: 开启警告
# -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT: 开启代码中的宏，使得两套表都能被加载
CFLAGS = -O3 -march=native -Wall -Wextra -fomit-frame-pointer \
         -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT

<<<<<<< HEAD
# 包含路径
INCLUDES = -Iinc

# 目录定义
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# 源文件和目标文件
SRCS = src/test.c src/poly.c
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
TARGET = $(BIN_DIR)/test_poly

# 默认目标
all: prepare $(TARGET)

# 创建目录
prepare:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

# 链接
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# 编译 .c -> .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# 清理
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
=======
# 编译选项
CFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops -Wall -Wextra -I./inc -I./test/inc
# 启用 Dilithium 和 Kyber NTT
CFLAGS += -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT

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
>>>>>>> new_structure

# 运行测试
run: all
	./$(TARGET)

.PHONY: all clean run prepare