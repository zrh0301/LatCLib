#include "test.h"
#include "kyber.h"
#include "kyber_pack.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void print_case(const char *name, int ok) {
    printf("  %-28s %s\n", name, ok ? "[PASS]" : "[FAIL]");
}

static void fill_poly(poly a, int seed) {
    for(int i = 0; i < KYBER_N; i++) {
        a[i] = (lssi_t)((seed * 97 + i * 31 + (i >> 1)) % KYBER_Q);
    }
}

static void fill_bytes(uint8_t *buf, int len, uint8_t seed) {
    for(int i = 0; i < len; i++) {
        buf[i] = (uint8_t)(seed + 13u * (uint8_t)i + (uint8_t)(i >> 2));
    }
}

static int poly_eq(const poly a, const poly b) {
    for(int i = 0; i < KYBER_N; i++) {
        if(a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static int test_byte_encode_decode_one(int d) {
    poly a, decoded;
    uint8_t buf[KYBER_POLYBYTES];
    int mask = (1 << d) - 1;

    memset(buf, 0, sizeof(buf));
    for(int i = 0; i < KYBER_N; i++) {
        a[i] = (d == KYBER_QBITS) ? (lssi_t)((i * 41 + 9) % KYBER_Q) : (lssi_t)((i * 41 + 9) & mask);
    }

    kyber_byte_encode(buf, a, d);
    kyber_byte_decode(decoded, buf, d);
    return poly_eq(a, decoded);
}

static int test_byte_encode_decode(void) {
    const int widths[] = {1, 4, 5, 10, 11, 12};
    int ok = 1;

    for(size_t i = 0; i < sizeof(widths) / sizeof(widths[0]); i++) {
        ok &= test_byte_encode_decode_one(widths[i]);
    }
    return ok;
}

static int test_ek_dk_pack_mode(int mode) {
    kyber_param param;
    uint8_t ek[KYBER_MAX_EK_BYTES];
    uint8_t dk[KYBER_MAX_DK_PKE_BYTES];
    uint8_t rho[KYBER_SYMBYTES];
    uint8_t decoded_rho[KYBER_SYMBYTES];
    poly t_hat[KYBER_MAX_K];
    poly decoded_t_hat[KYBER_MAX_K];
    poly s_hat[KYBER_MAX_K];
    poly decoded_s_hat[KYBER_MAX_K];
    int ok = 1;

    kyber_init_param(&param, mode);
    fill_bytes(rho, KYBER_SYMBYTES, (uint8_t)(mode & 0xFF));
    for(int i = 0; i < param.k; i++) {
        fill_poly(t_hat[i], mode + i);
        fill_poly(s_hat[i], mode + 17 + i);
    }

    kyber_pack_ek(ek, t_hat, rho, param);
    kyber_unpack_ek(decoded_t_hat, decoded_rho, ek, param);
    kyber_pack_dk_pke(dk, s_hat, param);
    kyber_unpack_dk_pke(decoded_s_hat, dk, param);

    ok &= memcmp(rho, decoded_rho, KYBER_SYMBYTES) == 0;
    ok &= kyber_check_ek_modulus(ek, param);
    for(int i = 0; i < param.k; i++) {
        ok &= poly_eq(t_hat[i], decoded_t_hat[i]);
        ok &= poly_eq(s_hat[i], decoded_s_hat[i]);
    }

    ek[0] = 0xFF;
    ek[1] = 0x0F;
    ok &= !kyber_check_ek_modulus(ek, param);
    return ok;
}

static int test_ct_pack_mode(int mode) {
    kyber_param param;
    uint8_t ct[KYBER_MAX_CT_BYTES];
    poly u[KYBER_MAX_K], decoded_u[KYBER_MAX_K];
    poly v, decoded_v;
    poly compressed;
    poly expected;

    kyber_init_param(&param, mode);
    for(int i = 0; i < param.k; i++) {
        fill_poly(u[i], mode + 31 + i);
    }
    fill_poly(v, mode + 63);

    kyber_pack_ct(ct, u, v, param);
    kyber_unpack_ct(decoded_u, decoded_v, ct, param);

    for(int i = 0; i < param.k; i++) {
        kyber_poly_compress(compressed, u[i], param.du);
        kyber_poly_decompress(expected, compressed, param.du);
        if(!poly_eq(expected, decoded_u[i])) {
            return 0;
        }
    }

    kyber_poly_compress(compressed, v, param.dv);
    kyber_poly_decompress(expected, compressed, param.dv);
    return poly_eq(expected, decoded_v);
}

static int test_mode_pack(int mode) {
    return test_ek_dk_pack_mode(mode) && test_ct_pack_mode(mode);
}

static int run_kyber_pack_suite(void) {
    int ok = 1;
    int case_ok;

    case_ok = test_byte_encode_decode();
    print_case("byte encode/decode", case_ok);
    ok &= case_ok;

    case_ok = test_mode_pack(KYBER_MODE_2);
    print_case("level 2 ek/dk/ct", case_ok);
    ok &= case_ok;

    case_ok = test_mode_pack(KYBER_MODE_3);
    print_case("level 3 ek/dk/ct", case_ok);
    ok &= case_ok;

    case_ok = test_mode_pack(KYBER_MODE_5);
    print_case("level 5 ek/dk/ct", case_ok);
    ok &= case_ok;

    return ok;
}

const test_pack_suite test_pack_kyber_suite = {"Kyber pack", run_kyber_pack_suite};
