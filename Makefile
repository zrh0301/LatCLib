# ==========================================
# Dilithium NTT Test Makefile
# ==========================================

# 编译器设置
CC = gcc

# 编译选项
# -O3: 开启最高级别优化 (性能测试必需)
# -I./inc: 指定头文件搜索路径
# -Wall -Wextra: 开启警告
CFLAGS = -O3 -I./inc -Wall -Wextra

# 目标可执行文件名
TARGET = test_ntt

# 源文件列表
# 注意：poly_ntt.c 不需要在此处列出，因为通过宏技巧
# 它已经被包含在 src/test.c 中了。如果单独列出会被编译两次导致符号冲突。
SRCS = src/test.c src/poly.c

# 自动生成对象文件名 (.c -> .o)
OBJS = $(SRCS:.c=.o)

# ==========================================
# 规则定义
# ==========================================

# 默认目标：构建并运行
all: $(TARGET) run

# 链接步骤
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

# 编译步骤 (通用规则: .c -> .o)
%.o: %.c
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# 运行测试
run: $(TARGET)
	@echo "Running test..."
	@./$(TARGET)

# 清理生成的文件
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)

# 伪目标声明 (防止目录下有名为 all/clean/run 的文件干扰)
.PHONY: all clean run