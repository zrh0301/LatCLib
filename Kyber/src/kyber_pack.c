#include "kyber_pack.h"
#include <string.h>

int kyber_param_is_valid(kyber_param param) {
    return (param.k >= 2 && param.k <= KYBER_MAX_K &&
            (param.mode == KYBER_MODE_2 ||
             param.mode == KYBER_MODE_3 ||
             param.mode == KYBER_MODE_5));
}

lssi_t kyber_reduce(lssi_t a) {
    a %= KYBER_Q;
    if(a < 0) {
        a += KYBER_Q;
    }
    return a;
}

lssi_t kyber_compress_one(lssi_t x, int d) {
    int32_t mask = (1 << d) - 1;
    x = kyber_reduce(x);
    return (lssi_t)(((((int64_t)x << d) + KYBER_Q / 2) / KYBER_Q) & mask);
}

lssi_t kyber_decompress_one(lssi_t y, int d) {
    return (lssi_t)(((int64_t)y * KYBER_Q + (1 << (d - 1))) >> d);
}

void kyber_poly_compress(poly r, const poly a, int d) {
    for(int i = 0; i < KYBER_N; i++) {
        r[i] = kyber_compress_one(a[i], d);
    }
}

void kyber_poly_decompress(poly r, const poly a, int d) {
    int32_t mask = (1 << d) - 1;
    for(int i = 0; i < KYBER_N; i++) {
        r[i] = kyber_decompress_one(a[i] & mask, d);
    }
}

void kyber_byte_encode(uint8_t *out, const poly a, int d) {
    size_t outlen = (size_t)KYBER_N * (size_t)d / 8;
    uint32_t mask = (d == 32) ? 0xFFFFFFFFu : ((1u << d) - 1u);
    size_t bitpos = 0;

    memset(out, 0, outlen);
    for(int i = 0; i < KYBER_N; i++) {
        uint32_t t = (d == KYBER_QBITS) ? (uint32_t)kyber_reduce(a[i])
                                        : ((uint32_t)a[i] & mask);
        for(int j = 0; j < d; j++, bitpos++) {
            if((t >> j) & 1u) {
                out[bitpos >> 3] |= (uint8_t)(1u << (bitpos & 7));
            }
        }
    }
}

void kyber_byte_decode(poly a, const uint8_t *in, int d) {
    uint32_t mask = (d == 32) ? 0xFFFFFFFFu : ((1u << d) - 1u);
    size_t bitpos = 0;

    for(int i = 0; i < KYBER_N; i++) {
        uint32_t t = 0;
        for(int j = 0; j < d; j++, bitpos++) {
            t |= ((uint32_t)((in[bitpos >> 3] >> (bitpos & 7)) & 1u)) << j;
        }
        a[i] = (lssi_t)(t & mask);
    }
}

void kyber_pack_ek(uint8_t *ek, const poly *t_hat, const uint8_t rho[KYBER_SYMBYTES],
                   kyber_param param) {
    uint8_t *p = ek;

    for(int i = 0; i < param.k; i++) {
        kyber_byte_encode(p, t_hat[i], KYBER_QBITS);
        p += KYBER_POLYBYTES;
    }
    memcpy(p, rho, KYBER_SYMBYTES);
}

void kyber_unpack_ek(poly *t_hat, uint8_t rho[KYBER_SYMBYTES], const uint8_t *ek,
                     kyber_param param) {
    const uint8_t *p = ek;

    for(int i = 0; i < param.k; i++) {
        kyber_byte_decode(t_hat[i], p, KYBER_QBITS);
        p += KYBER_POLYBYTES;
    }
    memcpy(rho, p, KYBER_SYMBYTES);
}

void kyber_pack_dk_pke(uint8_t *dk_pke, const poly *s_hat, kyber_param param) {
    uint8_t *p = dk_pke;

    for(int i = 0; i < param.k; i++) {
        kyber_byte_encode(p, s_hat[i], KYBER_QBITS);
        p += KYBER_POLYBYTES;
    }
}

void kyber_unpack_dk_pke(poly *s_hat, const uint8_t *dk_pke, kyber_param param) {
    const uint8_t *p = dk_pke;

    for(int i = 0; i < param.k; i++) {
        kyber_byte_decode(s_hat[i], p, KYBER_QBITS);
        p += KYBER_POLYBYTES;
    }
}

void kyber_pack_ct(uint8_t *ct, const poly *u, const poly v, kyber_param param) {
    uint8_t *p = ct;
    poly tmp;

    for(int i = 0; i < param.k; i++) {
        kyber_poly_compress(tmp, u[i], param.du);
        kyber_byte_encode(p, tmp, param.du);
        p += KYBER_N * param.du / 8;
    }

    kyber_poly_compress(tmp, v, param.dv);
    kyber_byte_encode(p, tmp, param.dv);
}

void kyber_unpack_ct(poly *u, poly v, const uint8_t *ct, kyber_param param) {
    const uint8_t *p = ct;
    poly tmp;

    for(int i = 0; i < param.k; i++) {
        kyber_byte_decode(tmp, p, param.du);
        kyber_poly_decompress(u[i], tmp, param.du);
        p += KYBER_N * param.du / 8;
    }

    kyber_byte_decode(tmp, p, param.dv);
    kyber_poly_decompress(v, tmp, param.dv);
}

int kyber_check_ek_modulus(const uint8_t *ek, kyber_param param) {
    const uint8_t *p = ek;
    poly t;

    for(int i = 0; i < param.k; i++) {
        kyber_byte_decode(t, p, KYBER_QBITS);
        for(int j = 0; j < KYBER_N; j++) {
            if(t[j] >= KYBER_Q) {
                return 0;
            }
        }
        p += KYBER_POLYBYTES;
    }
    return 1;
}
