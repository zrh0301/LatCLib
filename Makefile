# ==========================================
# My_PQC Ultimate Benchmark Makefile
# ==========================================

CC = gcc
# 关键：-O3 优化，-g 用于调试（如果需要），-I./inc 头文件
CFLAGS = -O3 -Wall -Wextra -I./inc

DEPS = inc/poly.h inc/types.h
COMMON_SRCS = src/poly.c src/poly_ntt.c
COMMON_OBJS = $(COMMON_SRCS:.c=.o)

# 目标: 全函数性能测试
TARGET_BENCH = bench_all
SRC_BENCH    = src/test_perf_all.c
OBJ_BENCH    = $(SRC_BENCH:.c=.o)

# ------------------------------------------

all: $(TARGET_BENCH)

# 链接
$(TARGET_BENCH): $(OBJ_BENCH) $(COMMON_OBJS)
	@echo "Linking $(TARGET_BENCH)..."
	$(CC) $(CFLAGS) -o $@ $^

# 编译
%.o: %.c $(DEPS)
	@echo "Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@

# 运行
test: $(TARGET_BENCH)
	@echo "\n>>> Running All Functions Performance Test <<<\n"
	@./$(TARGET_BENCH)

clean:
	rm -f src/*.o $(TARGET_BENCH)

.PHONY: all test clean