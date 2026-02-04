CC = gcc
CFLAGS = -O3 -Wall -Wextra

# 目标文件
TARGET = test_sha3
OBJS = test.o keccak.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

keccak.o: keccak.c keccak.h
	$(CC) $(CFLAGS) -c keccak.c

test.o: test.c keccak.h
	$(CC) $(CFLAGS) -c test.c

clean:
	rm -f $(OBJS) $(TARGET)