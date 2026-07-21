#include "test.h"
#include "dili.h"
#include "pack.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void print_case(const char *name, int ok) {
    printf("  %-28s %s\n", name, ok ? "[PASS]" : "[FAIL]");
}

static int poly_eq(const poly a, const poly b) {
    for(int i = 0; i < DILI_N; i++) {
        if(a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static int hint_eq(uint8_t a[][DILI_N], uint8_t b[][DILI_N], int k) {
    for(int i = 0; i < k; i++) {
        for(int j = 0; j < DILI_N; j++) {
            if(a[i][j] != b[i][j]) {
                return 0;
            }
        }
    }
    return 1;
}

static void fill_z(poly z, int gamma1) {
    for(int i = 0; i < DILI_N; i++) {
        z[i] = (i * 137 - 50000) % (2 * gamma1);
        if(z[i] > gamma1) {
            z[i] -= 2 * gamma1;
        }
    }
}

static int test_pack_z_case(int gamma1) {
    poly z, decoded;
    uint8_t packed[3 * DILI_N];

    memset(packed, 0, sizeof(packed));
    fill_z(z, gamma1);
    pack_z(packed, z, gamma1);
    unpack_z(decoded, packed, gamma1);
    return poly_eq(z, decoded);
}

static void put_hint(uint8_t h[][DILI_N], int row, int col) {
    if(row < DILI_K && col < DILI_N) {
        h[row][col] = 1;
    }
}

static int test_pack_h_case(int k, int omega) {
    uint8_t h[DILI_K][DILI_N];
    uint8_t decoded[DILI_K][DILI_N];
    uint8_t packed[DILI_OMEGA + DILI_K];

    memset(h, 0, sizeof(h));
    memset(decoded, 0, sizeof(decoded));
    memset(packed, 0, sizeof(packed));

    put_hint(h, 0, 10);
    put_hint(h, 0, 25);
    put_hint(h, 0, 100);
    put_hint(h, 1, 5);
    put_hint(h, 1, 50);
    put_hint(h, 2, 15);
    put_hint(h, 3, 200);
    if(k >= 6) {
        put_hint(h, 4, 30);
        put_hint(h, 5, 120);
    }
    if(k >= 8) {
        put_hint(h, 6, 60);
        put_hint(h, 7, 180);
    }

    if(pack_h(packed, h, k, omega) != 0) {
        return 0;
    }
    if(unpack_h(decoded, packed, k, omega) != 0) {
        return 0;
    }
    return hint_eq(h, decoded, k);
}

static int test_pack_w1_case(int gamma2) {
    poly w1;
    uint8_t packed1[DILI_N];
    uint8_t packed2[DILI_N];
    int w1_range = (DILI_Q - 1) / (2 * gamma2);
    size_t packed_len = (gamma2 == ((DILI_Q - 1) / 88)) ? (DILI_N * 3u / 4u) : (DILI_N / 2u);

    memset(packed1, 0, sizeof(packed1));
    memset(packed2, 0xA5, sizeof(packed2));
    for(int i = 0; i < DILI_N; i++) {
        w1[i] = i % w1_range;
    }

    pack_w1(packed1, w1, gamma2);
    pack_w1(packed2, w1, gamma2);
    return memcmp(packed1, packed2, packed_len) == 0;
}

static int test_polyt1_pack_case(void) {
    poly t1, decoded;
    uint8_t packed[DILI_POLYT1_PACKEDBYTES];

    for(int i = 0; i < DILI_N; i++) {
        t1[i] = (i * 17) & 0x3FF;
    }

    dili_polyt1_pack(packed, t1, DILI_N, DILI_D);
    dili_polyt1_unpack(decoded, packed, DILI_N, DILI_D);
    return poly_eq(t1, decoded);
}

static int test_pk_pack_case(void) {
    uint8_t rho[32];
    uint8_t decoded_rho[32];
    poly t1[DILI_K];
    poly decoded_t1[DILI_K];
    uint8_t pk[32 + DILI_K * DILI_POLYT1_PACKEDBYTES];

    for(int i = 0; i < 32; i++) {
        rho[i] = (uint8_t)(i * 7 + 123);
    }
    for(int k = 0; k < DILI_K; k++) {
        for(int i = 0; i < DILI_N; i++) {
            t1[k][i] = ((k * 256 + i) * 19) & 0x3FF;
        }
    }

    dili_pack_pk(pk, rho, t1, DILI_K);
    dili_unpack_pk(decoded_rho, decoded_t1, pk, DILI_K);
    if(memcmp(rho, decoded_rho, sizeof(rho)) != 0) {
        return 0;
    }
    for(int k = 0; k < DILI_K; k++) {
        if(!poly_eq(t1[k], decoded_t1[k])) {
            return 0;
        }
    }
    return 1;
}

static int run_dilithium_pack_suite(void) {
    int ok = 1;
    int case_ok;

    case_ok = test_pack_z_case(1 << 17);
    print_case("pack_z level 2", case_ok);
    ok &= case_ok;

    case_ok = test_pack_z_case(1 << 19);
    print_case("pack_z level 3/5", case_ok);
    ok &= case_ok;

    case_ok = test_pack_h_case(4, 80);
    print_case("pack_h level 2", case_ok);
    ok &= case_ok;

    case_ok = test_pack_h_case(6, 55);
    print_case("pack_h level 3", case_ok);
    ok &= case_ok;

    case_ok = test_pack_h_case(8, 75);
    print_case("pack_h level 5", case_ok);
    ok &= case_ok;

    case_ok = test_pack_w1_case((DILI_Q - 1) / 88);
    print_case("pack_w1 level 2", case_ok);
    ok &= case_ok;

    case_ok = test_pack_w1_case((DILI_Q - 1) / 32);
    print_case("pack_w1 level 3/5", case_ok);
    ok &= case_ok;

    case_ok = test_polyt1_pack_case();
    print_case("polyt1 pack", case_ok);
    ok &= case_ok;

    case_ok = test_pk_pack_case();
    print_case("public key pack", case_ok);
    ok &= case_ok;

    return ok;
}

const test_pack_suite test_pack_dilithium_suite = {"Dilithium pack", run_dilithium_pack_suite};
