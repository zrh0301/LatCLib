CC = gcc

# 编译选项
# -O3: 最高级优化
# -march=native: 针对本机CPU指令集优化
# -ffast-math -funroll-loops: 激进的数学和循环优化
# -Wall -Wextra: 开启警告
CFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops -g -Wall -Wextra
# 添加头文件搜索路径
CFLAGS += -Iinc -Itest/inc 

# 宏定义：启用 Dilithium 和 Kyber NTT
CFLAGS += -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT

# 链接选项 (Linux 下使用 math 库必须显式链接 -lm)
LDFLAGS = -lm -g

# 目录定义
SRC_DIR = src
TEST_SRC_DIR = test/src
OBJ_DIR = obj
BIN_DIR = bin

# ==========================================
# 1. 核心源文件
# ==========================================
CORE_SRCS = $(SRC_DIR)/poly.c $(SRC_DIR)/poly_ntt.c $(SRC_DIR)/rand.c $(SRC_DIR)/keccak.c $(SRC_DIR)/dili.c $(SRC_DIR)/polyvec.c $(SRC_DIR)/pack.c
# 将 src/xxx.c 映射为 obj/xxx.o
CORE_OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(CORE_SRCS))

# 测试源文件
TEST_POLY_SRCS = $(TEST_SRC_DIR)/test_poly.c
TEST_NTT_SRCS = $(TEST_SRC_DIR)/test_ntt.c
TEST_SHA_SRCS = $(TEST_SRC_DIR)/test_sha3.c
TEST_SAMPLING_SRCS = $(TEST_SRC_DIR)/test_sampling.c
TEST_DILI_SRCS = $(TEST_SRC_DIR)/test_dili.c

# 测试对象文件 (映射 test/src/xxx.c -> obj/xxx.o)
TEST_POLY_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_POLY_SRCS))
TEST_NTT_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_NTT_SRCS))
TEST_SHA_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_SHA_SRCS))
TEST_SAMPLING_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_SAMPLING_SRCS))
TEST_DILI_OBJS = $(patsubst $(TEST_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_DILI_SRCS))

# 目标文件
TARGET_POLY = $(BIN_DIR)/test_poly
TARGET_NTT = $(BIN_DIR)/test_ntt
TARGET_SHA = $(BIN_DIR)/test_sha3
TARGET_SAMPLING = $(BIN_DIR)/test_sampling
TARGET_DILI = $(BIN_DIR)/test_dili

# 头文件依赖 (用于检测更新)
DEPS = inc/poly.h inc/types.h inc/rand.h inc/keccak.h inc/dili.h test/inc/test_common.h 

# 默认目标
all: prepare $(TARGET_POLY) $(TARGET_NTT) $(TARGET_SAMPLING) $(TARGET_DILI)

# 创建目录 (Linux 格式修改点)
# mkdir -p: 如果目录存在不报错，不存在则创建，且支持创建多级目录
prepare:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

# 链接 test_poly
$(TARGET_POLY): $(CORE_OBJS) $(TEST_POLY_OBJS)
	@echo "Linking $(TARGET_POLY)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_ntt
$(TARGET_NTT): $(CORE_OBJS) $(TEST_NTT_OBJS)
	@echo "Linking $(TARGET_NTT)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_sha3
$(TARGET_SHA): $(CORE_OBJS) $(TEST_SHA_OBJS)
	@echo "Linking $(TARGET_SHA)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_sampling
$(TARGET_SAMPLING): $(CORE_OBJS) $(TEST_SAMPLING_OBJS)
	@echo "Linking $(TARGET_SAMPLING)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 链接 test_dili
$(TARGET_DILI): $(CORE_OBJS) $(TEST_DILI_OBJS)
	@echo "Linking $(TARGET_DILI)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 编译 src 目录下的源文件
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(DEPS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 编译 test/src 目录下的源文件
$(OBJ_DIR)/%.o: $(TEST_SRC_DIR)/%.c $(DEPS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 测试目标
test: all
	@echo ""
	@echo "============================================="
	@echo "   Running Dilithium Tests"
	@echo "============================================="
	./$(TARGET_DILI)

# 清理 (Linux 格式修改点)
# rm -rf: 强制递归删除，不询问
clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	@echo "Cleaned build directories."

# 伪目标
.PHONY: all clean test prepare