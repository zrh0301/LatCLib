#include "kyber_sample.h"
#include "kyber_pack.h"
#include "keccak.h"
#include <string.h>

void kyber_g(uint8_t out[64], const uint8_t *in, size_t inlen) {
    sha3_512(out, in, inlen);
}

void kyber_prf(uint8_t *out, size_t outlen, const uint8_t key[KYBER_SYMBYTES],
               uint8_t nonce) {
    uint8_t extkey[KYBER_SYMBYTES + 1];

    memcpy(extkey, key, KYBER_SYMBYTES);
    extkey[KYBER_SYMBYTES] = nonce;
    shake256(out, outlen, extkey, sizeof(extkey));
}

void kyber_sample_ntt(poly a, const uint8_t seed[34], lp_ctx *ctx) {
    KECCAK1600_CTX kc;
    uint8_t buf[3];
    int ctr = 0;
    (void)ctx;

    shake128_init(&kc);
    sha3_update(&kc, seed, 34);
    shake_xof_final(&kc);

    while(ctr < KYBER_N) {
        uint32_t d1;
        uint32_t d2;

        shake_squeeze(&kc, buf, sizeof(buf));
        d1 = (uint32_t)buf[0] | (((uint32_t)buf[1] & 0x0Fu) << 8);
        d2 = ((uint32_t)buf[1] >> 4) | ((uint32_t)buf[2] << 4);

        if(d1 < KYBER_Q) {
            a[ctr++] = (lssi_t)d1;
        }
        if(ctr < KYBER_N && d2 < KYBER_Q) {
            a[ctr++] = (lssi_t)d2;
        }
    }
}

static int kyber_get_bit(const uint8_t *buf, int bitidx) {
    return (buf[bitidx >> 3] >> (bitidx & 7)) & 1;
}

void kyber_sample_cbd(poly a, const uint8_t *buf, int eta) {
    for(int i = 0; i < KYBER_N; i++) {
        int x = 0;
        int y = 0;
        int base = 2 * eta * i;

        for(int j = 0; j < eta; j++) {
            x += kyber_get_bit(buf, base + j);
            y += kyber_get_bit(buf, base + eta + j);
        }
        a[i] = kyber_reduce((lssi_t)(x - y));
    }
}

void kyber_sample_noise(poly a, const uint8_t sigma[KYBER_SYMBYTES],
                        uint8_t nonce, int eta) {
    uint8_t buf[64 * 3];
    size_t outlen = (size_t)64 * (size_t)eta;

    kyber_prf(buf, outlen, sigma, nonce);
    kyber_sample_cbd(a, buf, eta);
}
