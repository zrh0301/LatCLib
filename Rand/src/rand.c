
#include "rand.h"
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define ROTL32(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

#define QR(a, b, c, d) do { \
	a += b; d ^= a; d = ROTL32(d, 16); \
	c += d; b ^= c; b = ROTL32(b, 12); \
	a += b; d ^= a; d = ROTL32(d, 8);  \
	c += d; b ^= c; b = ROTL32(b, 7);  \
} while(0)

typedef struct {
	uint32_t state[16];
	uint8_t keystream[64];
	int available;
	int initialized;
} stream_cipher_ctx;

static stream_cipher_ctx g_rng = {0};

static void chacha20_block(uint32_t out[16], const uint32_t in[16]) {
	int i;
	for (i = 0; i < 16; i++) {
		out[i] = in[i];
	}
	for (i = 0; i < 10; i++) {
		QR(out[0], out[4], out[8],  out[12]);
		QR(out[1], out[5], out[9],  out[13]);
		QR(out[2], out[6], out[10], out[14]);
		QR(out[3], out[7], out[11], out[15]);
		QR(out[0], out[5], out[10], out[15]);
		QR(out[1], out[6], out[11], out[12]);
		QR(out[2], out[7], out[8],  out[13]);
		QR(out[3], out[4], out[9],  out[14]);
	}
	for (i = 0; i < 16; i++) {
		out[i] += in[i];
	}
}

void rand_seed(uint8_t* buf, int size)
{
	unsigned int seed;
	int i;
	uint8_t key_material[32];
	seed = (unsigned int)time(NULL);
	for (i = 0; i < size && i < 32; i++)
	{
		seed += buf[i] << ((i * 3) % 0x1F);
		key_material[i] = (seed >> ((i % 4) * 8)) & 0xFF;
	}
	for (; i < 32; i++) {
		seed = seed * 1103515245 + 12345;
		key_material[i] = (seed >> 16) & 0xFF;
	}
	g_rng.state[0] = 0x61707865;
	g_rng.state[1] = 0x3320646e;
	g_rng.state[2] = 0x79622d32;
	g_rng.state[3] = 0x6b206574;
	for (i = 0; i < 8; i++) {
		g_rng.state[4 + i] = 
			((uint32_t)key_material[i*4+0]) |
			((uint32_t)key_material[i*4+1] << 8) |
			((uint32_t)key_material[i*4+2] << 16) |
			((uint32_t)key_material[i*4+3] << 24);
	}
	
	g_rng.state[12] = 0;
	
	g_rng.state[13] = (uint32_t)time(NULL);
	g_rng.state[14] = seed;
	g_rng.state[15] = seed ^ 0xDEADBEEF;
	
	g_rng.available = 0;
	g_rng.initialized = 1;
}

void rand_bytes(uint8_t* buf, int size)
{
	int i, j;
	uint32_t block[16];
	
	if (!g_rng.initialized) {
		uint8_t default_seed[4] = {0x12, 0x34, 0x56, 0x78};
		rand_seed(default_seed, 4);
	}

	for (i = 0; i < size; i++)
	{
		if (g_rng.available == 0) {
			chacha20_block(block, g_rng.state);
			g_rng.state[12]++;
			if (g_rng.state[12] == 0) {
				g_rng.state[13]++;
			}
			for (j = 0; j < 16; j++) {
				g_rng.keystream[j*4+0] = block[j] & 0xFF;
				g_rng.keystream[j*4+1] = (block[j] >> 8) & 0xFF;
				g_rng.keystream[j*4+2] = (block[j] >> 16) & 0xFF;
				g_rng.keystream[j*4+3] = (block[j] >> 24) & 0xFF;
			}
			g_rng.available = 64;
		}
		
		buf[i] = g_rng.keystream[64 - g_rng.available];
		g_rng.available--;
	}
}

void rand_bits(uint8_t* buf, int bits)
{
	int i;
	int bytes = bits >> 3;
	uint8_t resmask = (1 << (bits % 8)) - 1;
	if (bytes > 0) {
		rand_bytes(buf, bytes);
	}
	if (bits % 8) {
		uint8_t temp;
		rand_bytes(&temp, 1);
		buf[bytes] = temp & resmask;
	}
}