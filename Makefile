# ==========================================
# My_PQC Full Suite Makefile
# ==========================================

CC = gcc

# 编译选项
# -O3    : 最高级优化 (性能测试必需)
# -Wall  : 开启警告
# -I./inc: 头文件目录
CFLAGS = -O3 -Wall -Wextra -I./inc

# 目标
TARGET = run_full_test

# 源文件 (包含核心库和新的测试文件)
SRCS = src/poly.c src/poly_ntt.c src/test_poly.c

# 自动生成 .o 文件名
OBJS = $(SRCS:.c=.o)

# 头文件依赖
DEPS = inc/poly.h inc/types.h

# 默认目标
all: $(TARGET)

# 链接
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CC) $(CFLAGS) -o $@ $^

# 编译 (.c -> .o)
%.o: %.c $(DEPS)
	@echo "Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@

# 运行测试
test: $(TARGET)
	@echo "\n============================================="
	@echo ">>> Running Complete Correctness & Benchmark <<<"
	@echo "=============================================\n"
	@./$(TARGET)

# 清理
clean:
	rm -f src/*.o $(TARGET)

.PHONY: all test clean