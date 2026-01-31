CC = gcc
# 编译选项：
# -O3: 最高级优化
# -march=native: 针对本机CPU指令集优化 (对NTT性能至关重要)
# -Wall -Wextra: 开启警告
# -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT: 开启代码中的宏，使得两套表都能被加载
CFLAGS = -O3 -march=native -Wall -Wextra -fomit-frame-pointer \
         -DUSE_DILITHIUM_NTT -DUSE_KYBER_NTT

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

# 运行测试
run: all
	./$(TARGET)

.PHONY: all clean run prepare