CC = gcc

# 编译选项
# -O3: 最高级优化
# -march=native: 针对本机CPU指令集优化
# -ffast-math -funroll-loops: 激进的数学和循环优化
# -Wall -Wextra: 开启警告
CFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops -Wall -Wextra
CFLAGS += -Iinc -Itest/inc

# 宏定义：启用 Dilithium 和 Kyber NTT
CFLAGS += -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT

# Linux 下通常不需要链接 advapi32，如果你的 rand.c 需要 math 库，可以解开下面注释
# LDFLAGS = -lm

# 目录定义
SRC_DIR = src
TEST_SRC_DIR = test/src
OBJ_DIR = obj
BIN_DIR = bin

# 核心库源文件
CORE_SRCS = $(SRC_DIR)/poly.c $(SRC_DIR)/poly_ntt.c $(SRC_DIR)/rand.c
# 将 src/xxx.c 映射为 obj/xxx.o
CORE_OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(CORE_SRCS))

# 测试源文件
TEST_POLY_SRCS = $(TEST_SRC_DIR)/test_poly.c
TEST_NTT_SRCS = $(TEST_SRC_DIR)/test_ntt.c
TEST_SAMPLING_SRCS = $(TEST_SRC_DIR)/test_sampling.c

# 测试对象文件 (映射 test/src/xxx.c -> obj/xxx.o)
TEST_POLY_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_POLY_SRCS))
TEST_NTT_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_NTT_SRCS))
TEST_SAMPLING_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_SAMPLING_SRCS))

# 目标文件
TARGET_POLY = $(BIN_DIR)/test_poly
TARGET_NTT = $(BIN_DIR)/test_ntt
TARGET_SAMPLING = $(BIN_DIR)/test_sampling

# 头文件依赖 (用于检测更新)
DEPS = inc/poly.h inc/types.h inc/rand.h test/inc/test_common.h

# 默认目标
all: prepare $(TARGET_POLY) $(TARGET_NTT) $(TARGET_SAMPLING)

# 创建目录
prepare:
	@if not exist $(OBJ_DIR) mkdir $(OBJ_DIR)
	@if not exist $(BIN_DIR) mkdir $(BIN_DIR)

# 链接 test_poly
$(TARGET_POLY): $(CORE_OBJS) $(TEST_POLY_OBJS)
	@echo "Linking $(TARGET_POLY)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_ntt
$(TARGET_NTT): $(CORE_OBJS) $(TEST_NTT_OBJS)
	@echo "Linking $(TARGET_NTT)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_sampling
$(TARGET_SAMPLING): $(CORE_OBJS) $(TEST_SAMPLING_OBJS)
	@echo "Linking $(TARGET_SAMPLING)..."
	$(CC) $(CFLAGS) -o $@ $^ -lm $(LDFLAGS)

# 编译 src 目录下的源文件
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(DEPS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 编译 test/src 目录下的源文件
$(OBJ_DIR)/%.o: $(TEST_SRC_DIR)/%.c $(DEPS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 运行所有测试
test: all
	@echo ""
	@echo "============================================="
	@echo "   Running Polynomial Tests"
	@echo "============================================="
	@echo ""
	./$(TARGET_POLY)
	@echo ""
	@echo "============================================="
	@echo "   Running NTT Tests"
	@echo "============================================="
	@echo ""
	./$(TARGET_NTT)
	@echo ""
	@echo "============================================="
	@echo "   Running Sampling Tests"
	@echo "============================================="
	@echo ""
	./$(TARGET_SAMPLING)

# 清理
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# 伪目标
.PHONY: all clean test prepare