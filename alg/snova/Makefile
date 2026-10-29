ifeq ($(ARCH),)
ARCH = native
endif

ifeq ($(OPT),)
OPT = OPT
endif

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -fstack-usage -O3 -march=$(ARCH) -mtune=$(ARCH) -I../../inc
# CFLAGS += -g
# CFLAGS += -fsanitize=address -fsanitize=undefined -fsanitize-address-use-after-scope

# Select appropriate source file
ifeq ($(OPT), REF)
OBJS = ./build/snova_ref.o ./build/symmetric_ref.o
CFLAGS += -DSNOVA_OPT=0
else ifeq ($(OPT), OPT)
OBJS = ./build/snova_opt_16.o ./build/symmetric.o
CFLAGS += -DSNOVA_OPT=10
else ifeq ($(OPT), AVX2)
OBJS = ./build/snova_avx2_16.o ./build/symmetric.o
CFLAGS += -DSNOVA_OPT=20
else ifeq ($(OPT), MEM)
OBJS = ./build/snova_memopt_16.o ./build/symmetric.o
CFLAGS += -DSNOVA_OPT=5
else ifeq ($(OPT), TINYMEM)
OBJS = ./build/snova_tinymem_16.o ./build/symmetric_tiny.o
CFLAGS += -DSNOVA_OPT=5
else
$(error OPT=$(OPT) is not supported)
endif

OBJS += ./build/rng.o ./build/sign.o ./build/aes.o
OBJS += ./build/pqc_keccak.o ./build/pqc_rand.o ./build/pqc_gf16.o ./build/pqc_gf16_matrix.o
HDRS = api.h rng.h symmetric.h snova.h snova_params.h
# LIBS = -lcrypto


all: build kat speed

clean:
	rm -rf build *.req *.rsp *.log

build:
	mkdir -p build

./build/%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(P) -o $@ -c $<

./build/pqc_%.o: ../../src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

./build/PQCgenKAT_sign: $(OBJS) ./build/PQCgenKAT_sign.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

kat: build ./build/PQCgenKAT_sign
	./build/PQCgenKAT_sign

./build/speed: $(OBJS) ./build/speed.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

speed: build ./build/speed
	@./build/speed
