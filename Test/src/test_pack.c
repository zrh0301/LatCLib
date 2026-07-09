#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "poly.h"
#include "dili.h"
#include "pack.h"
#include "kyber.h"
#include "kyber_pack.h"
#include "types.h"

void print_bytes(const char* label, const uint8_t* data, int len) {
    printf("%s: ", label);
    for(int i = 0; i < len && i < 32; i++) {
        printf("%02x ", data[i]);
    }
    if(len > 32) printf("...");
    printf("\n");
}

void print_poly(const char* label, const poly p, int n) {
    printf("%s: ", label);
    for(int i = 0; i < n && i < 8; i++) {
        printf("%d ", p[i]);
    }
    if(n > 8) printf("...");
    printf("\n");
}

void test_pack_unpack_z_with_params(int gamma1, const char* mode_name) {
    printf("\n========================================\n");
    printf("Test: pack_z / unpack_z (%s, gamma1=%d)\n", mode_name, gamma1);
    printf("========================================\n");
    
    poly z_original, z_unpacked;
    uint8_t packed[3 * DILI_N];
    
    // Initialize test data
    for(int i = 0; i < DILI_N; i++) {
        z_original[i] = (i * 137 - 50000) % (2 * gamma1);
        if(z_original[i] > gamma1) z_original[i] -= 2 * gamma1;
    }
    
    printf("\nOriginal z data:\n");
    print_poly("  z_original", z_original, 8);
    printf("  z_original[0] = %d\n", z_original[0]);
    printf("  z_original[1] = %d\n", z_original[1]);
    printf("  z_original[255] = %d\n", z_original[255]);
    
    // Pack
    pack_z(packed, z_original, gamma1);
    printf("\nPacked bytes:\n");
    print_bytes("  packed", packed, 3 * DILI_N);
    printf("  packed[0] = 0x%02x\n", packed[0]);
    printf("  packed[1] = 0x%02x\n", packed[1]);
    printf("  packed[2] = 0x%02x\n", packed[2]);
    
    // Unpack
    unpack_z(z_unpacked, packed, gamma1);
    printf("\nUnpacked z data:\n");
    print_poly("  z_unpacked", z_unpacked, 8);
    printf("  z_unpacked[0] = %d\n", z_unpacked[0]);
    printf("  z_unpacked[1] = %d\n", z_unpacked[1]);
    printf("  z_unpacked[255] = %d\n", z_unpacked[255]);
    
    // Verify
    int errors = 0;
    for(int i = 0; i < DILI_N; i++) {
        if(z_original[i] != z_unpacked[i]) {
            errors++;
            if(errors <= 5) {
                printf("  Error at index %d: original=%d, unpacked=%d\n", 
                       i, z_original[i], z_unpacked[i]);
            }
        }
    }
    
    if(errors == 0) {
        printf("\n[PASS] %s pack_z/unpack_z test passed!\n", mode_name);
    } else {
        printf("\n[FAIL] %s pack_z/unpack_z has %d errors!\n", mode_name, errors);
    }
}

void test_pack_unpack_z() {
    // Test Dilithium2
    test_pack_unpack_z_with_params(1 << 17, "Dilithium2");
    
    // Test Dilithium3/5
    test_pack_unpack_z_with_params(1 << 19, "Dilithium3/5");
}

void test_pack_unpack_h_with_params(int k, int omega, const char* mode_name) {
    printf("\n========================================\n");
    printf("Test: pack_h / unpack_h (%s, k=%d, omega=%d)\n", mode_name, k, omega);
    printf("========================================\n");
    
    uint8_t h_original[DILI_K][DILI_N];
    uint8_t h_unpacked[DILI_K][DILI_N];
    uint8_t packed[DILI_OMEGA + DILI_K];
    
    // Initialize test data - sparse hint vectors
    for(int i = 0; i < DILI_K; i++) {
        for(int j = 0; j < DILI_N; j++) {
            h_original[i][j] = 0;
        }
    }
    
    // Set some hint bits
    h_original[0][10] = 1;
    h_original[0][25] = 1;
    h_original[0][100] = 1;
    h_original[1][5] = 1;
    h_original[1][50] = 1;
    h_original[2][15] = 1;
    h_original[3][200] = 1;
    if(k >= 6) {
        h_original[4][30] = 1;
        h_original[5][120] = 1;
    }
    if(k >= 8) {
        h_original[6][60] = 1;
        h_original[7][180] = 1;
    }
    
    printf("\nOriginal hint data:\n");
    int total_ones = 0;
    for(int i = 0; i < k; i++) {
        int ones = 0;
        for(int j = 0; j < DILI_N; j++) {
            if(h_original[i][j] == 1) {
                ones++;
                if(ones <= 3) {
                    printf("  h[%d][%d] = 1\n", i, j);
                }
            }
        }
        total_ones += ones;
        printf("  h[%d] has %d ones\n", i, ones);
    }
    printf("  Total: %d ones (omega=%d)\n", total_ones, omega);
    
    // Pack
    int pack_result = pack_h(packed, h_original, k, omega);
    printf("\nPack result: %d\n", pack_result);
    print_bytes("  packed", packed, omega + k);
    
    // Unpack
    int unpack_result = unpack_h(h_unpacked, packed, k, omega);
    printf("\nUnpacked hint data:\n");
    for(int i = 0; i < k; i++) {
        int ones = 0;
        for(int j = 0; j < DILI_N; j++) {
            if(h_unpacked[i][j] == 1) {
                ones++;
                if(ones <= 3) {
                    printf("  h[%d][%d] = 1\n", i, j);
                }
            }
        }
        printf("  h[%d] has %d ones\n", i, ones);
    }
    
    // Verify
    int errors = 0;
    for(int i = 0; i < k; i++) {
        for(int j = 0; j < DILI_N; j++) {
            if(h_original[i][j] != h_unpacked[i][j]) {
                errors++;
                if(errors <= 5) {
                    printf("  Error at h[%d][%d]: original=%d, unpacked=%d\n", 
                           i, j, h_original[i][j], h_unpacked[i][j]);
                }
            }
        }
    }
    
    if(errors == 0 && pack_result == 0 && unpack_result == 0) {
        printf("\n[PASS] %s pack_h/unpack_h test passed!\n", mode_name);
    } else {
        printf("\n[FAIL] %s pack_h/unpack_h has %d errors! (pack=%d, unpack=%d)\n", 
               mode_name, errors, pack_result, unpack_result);
    }
}

void test_pack_unpack_h() {
    // Test Dilithium2
    test_pack_unpack_h_with_params(4, 80, "Dilithium2");
    
    // Test Dilithium3
    test_pack_unpack_h_with_params(6, 55, "Dilithium3");
    
    // Test Dilithium5
    test_pack_unpack_h_with_params(8, 75, "Dilithium5");
}

void test_pack_w1_with_params(int gamma2, const char* mode_name) {
    printf("\n========================================\n");
    printf("Test: pack_w1 (%s, gamma2=%d)\n", mode_name, gamma2);
    printf("========================================\n");
    
    poly w1;
    uint8_t packed[DILI_N];
    
    // Calculate range based on gamma2
    int w1_range = (DILI_Q - 1) / (2 * gamma2);
    
    // Initialize test data
    for(int i = 0; i < DILI_N; i++) {
        w1[i] = i % w1_range;
    }
    
    printf("\nOriginal w1 data (range [0, %d)):\n", w1_range);
    print_poly("  w1", w1, 16);
    printf("  w1[0] = %d\n", w1[0]);
    printf("  w1[1] = %d\n", w1[1]);
    printf("  w1[2] = %d\n", w1[2]);
    printf("  w1[3] = %d\n", w1[3]);
    
    // Pack
    pack_w1(packed, w1, gamma2);
    printf("\nPacked bytes:\n");
    print_bytes("  packed", packed, 32);
    printf("  packed[0] = 0x%02x\n", packed[0]);
    printf("  packed[1] = 0x%02x\n", packed[1]);
    printf("  packed[2] = 0x%02x\n", packed[2]);
    
    printf("\n[INFO] %s pack_w1 completed (w1_range=%d)\n", mode_name, w1_range);
}

void test_pack_w1() {
    // Test Dilithium2 (gamma2 = (q-1)/88 = 95232, w1 in [0, 44))
    test_pack_w1_with_params((DILI_Q - 1) / 88, "Dilithium2");
    
    // Test Dilithium3/5 (gamma2 = (q-1)/32 = 261888, w1 in [0, 16))
    test_pack_w1_with_params((DILI_Q - 1) / 32, "Dilithium3/5");
}

void test_pack_unpack_t1() {
    printf("\n========================================\n");
    printf("Test: dili_polyt1_pack / dili_polyt1_unpack\n");
    printf("========================================\n");
    
    poly t1_original, t1_unpacked;
    // t1 uses 10-bit standard packing: 256 coefficients -> 320 bytes.
    uint8_t packed[DILI_POLYT1_PACKEDBYTES];

    // Initialize test data in the 10-bit t1 range [0, 1023].
    for(int i = 0; i < DILI_N; i++) {
        t1_original[i] = (i * 17) & 0x3FF;
    }
    
    printf("\nOriginal t1 data:\n");
    print_poly("  t1_original", t1_original, 8);
    printf("  t1_original[0] = %d\n", t1_original[0]);
    printf("  t1_original[1] = %d\n", t1_original[1]);
    printf("  t1_original[255] = %d\n", t1_original[255]);
    
    // Pack
    dili_polyt1_pack(packed, t1_original, DILI_N, DILI_D);
    printf("\nPacked bytes:\n");
    print_bytes("  packed", packed, DILI_POLYT1_PACKEDBYTES);
    
    // Unpack
    dili_polyt1_unpack(t1_unpacked, packed, DILI_N, DILI_D);
    printf("\nUnpacked t1 data:\n");
    print_poly("  t1_unpacked", t1_unpacked, 8);
    printf("  t1_unpacked[0] = %d\n", t1_unpacked[0]);
    printf("  t1_unpacked[1] = %d\n", t1_unpacked[1]);
    printf("  t1_unpacked[255] = %d\n", t1_unpacked[255]);
    
    // Verify
    int errors = 0;
    for(int i = 0; i < DILI_N; i++) {
        int32_t expected = (i * 17) & 0x3FF;
        
        if(t1_unpacked[i] != expected) {
            errors++;
            if(errors <= 5) {
                printf("  Error at index %d: expected=%d, unpacked=%d\n", 
                       i, expected, t1_unpacked[i]);
            }
        }
    }
    
    if(errors == 0) {
        printf("\n[PASS] dili_polyt1_pack/unpack test passed!\n");
    } else {
        printf("\n[FAIL] dili_polyt1_pack/unpack has %d errors!\n", errors);
    }
}

void test_pack_unpack_pk() {
    printf("\n========================================\n");
    printf("Test: dili_pack_pk / dili_unpack_pk\n");
    printf("========================================\n");
    
    const int K = DILI_K;
    uint8_t rho_original[32];
    uint8_t rho_unpacked[32];
    poly t1_original[K];
    poly t1_unpacked[K];
    
    // Public key size: rho + K * packed t1.
    uint8_t pk_bytes[32 + DILI_K * DILI_POLYT1_PACKEDBYTES];
    
    // Initialize test data
    printf("\nInitializing test data...\n");
    
    // Initialize rho (32-byte seed)
    for(int i = 0; i < 32; i++) {
        rho_original[i] = (uint8_t)(i * 7 + 123);
    }
    printf("  rho initialized: ");
    for(int i = 0; i < 8; i++) {
        printf("%02x ", rho_original[i]);
    }
    printf("...\n");
    
    // Initialize t1 vectors
    for(int k = 0; k < K; k++) {
        for(int i = 0; i < DILI_N; i++) {
            t1_original[k][i] = ((k * 256 + i) * 19) & 0x3FF;
        }
        printf("  t1[%d][0..3] = %d %d %d %d\n", 
               k, t1_original[k][0], t1_original[k][1], 
               t1_original[k][2], t1_original[k][3]);
    }
    
    // Pack public key
    printf("\nPacking public key...\n");
    dili_pack_pk(pk_bytes, rho_original, t1_original, K);
    printf("  pk_bytes[0..7]: ");
    for(int i = 0; i < 8; i++) {
        printf("%02x ", pk_bytes[i]);
    }
    printf("...\n");
    printf("  pk_bytes[32..39]: ");
    for(int i = 32; i < 40; i++) {
        printf("%02x ", pk_bytes[i]);
    }
    printf("...\n");
    
    // Unpack public key
    printf("\nUnpacking public key...\n");
    dili_unpack_pk(rho_unpacked, t1_unpacked, pk_bytes, K);
    
    printf("  rho unpacked: ");
    for(int i = 0; i < 8; i++) {
        printf("%02x ", rho_unpacked[i]);
    }
    printf("...\n");
    
    for(int k = 0; k < K; k++) {
        printf("  t1[%d][0..3] = %d %d %d %d\n", 
               k, t1_unpacked[k][0], t1_unpacked[k][1], 
               t1_unpacked[k][2], t1_unpacked[k][3]);
    }
    
    // Verify rho
    printf("\nVerifying rho...\n");
    int rho_errors = 0;
    for(int i = 0; i < 32; i++) {
        if(rho_original[i] != rho_unpacked[i]) {
            rho_errors++;
            if(rho_errors <= 3) {
                printf("  Error at rho[%d]: original=0x%02x, unpacked=0x%02x\n",
                       i, rho_original[i], rho_unpacked[i]);
            }
        }
    }
    
    // Verify t1
    printf("Verifying t1 vectors...\n");
    int t1_errors = 0;
    for(int k = 0; k < K; k++) {
        for(int i = 0; i < DILI_N; i++) {
            if(t1_original[k][i] != t1_unpacked[k][i]) {
                t1_errors++;
                if(t1_errors <= 5) {
                    printf("  Error at t1[%d][%d]: original=%d, unpacked=%d\n",
                           k, i, t1_original[k][i], t1_unpacked[k][i]);
                }
            }
        }
    }
    
    if(rho_errors == 0 && t1_errors == 0) {
        printf("\n[PASS] dili_pack_pk/unpack_pk test passed!\n");
    } else {
        printf("\n[FAIL] dili_pack_pk/unpack_pk has errors: rho=%d, t1=%d\n", 
               rho_errors, t1_errors);
    }
}

static void fill_kyber_poly(poly a, int seed) {
    for(int i = 0; i < KYBER_N; i++) {
        a[i] = (lssi_t)((seed * 97 + i * 31 + (i >> 1)) % KYBER_Q);
    }
}

static void fill_kyber_bytes(uint8_t *buf, int len, uint8_t seed) {
    for(int i = 0; i < len; i++) {
        buf[i] = (uint8_t)(seed + 13u * (uint8_t)i + (uint8_t)(i >> 2));
    }
}

static int kyber_poly_match(const poly a, const poly b) {
    for(int i = 0; i < KYBER_N; i++) {
        if(a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static int test_kyber_byte_encode_decode_one(int d) {
    poly a, b;
    uint8_t buf[KYBER_POLYBYTES];
    int mask = (1 << d) - 1;

    for(int i = 0; i < KYBER_N; i++) {
        a[i] = (d == KYBER_QBITS) ? (lssi_t)((i * 41 + 9) % KYBER_Q)
                                  : (lssi_t)((i * 41 + 9) & mask);
    }

    kyber_byte_encode(buf, a, d);
    kyber_byte_decode(b, buf, d);

    for(int i = 0; i < KYBER_N; i++) {
        if(a[i] != b[i]) {
            printf("  [FAIL] Kyber ByteEncode%d mismatch at %d: %d != %d\n",
                   d, i, a[i], b[i]);
            return 0;
        }
    }
    return 1;
}

static int test_kyber_byte_encode_decode(void) {
    int ds[] = {1, 4, 5, 10, 11, 12};
    int ok = 1;

    printf("\n========================================\n");
    printf("Test: Kyber ByteEncode / ByteDecode\n");
    printf("========================================\n");

    for(int i = 0; i < (int)(sizeof(ds) / sizeof(ds[0])); i++) {
        ok &= test_kyber_byte_encode_decode_one(ds[i]);
    }

    printf("%s Kyber ByteEncode / ByteDecode tests\n", ok ? "[PASS]" : "[FAIL]");
    return ok;
}

static int test_kyber_ek_dk_pack_mode(int mode) {
    kyber_param param;
    uint8_t ek[KYBER_MAX_EK_BYTES];
    uint8_t dk[KYBER_MAX_DK_PKE_BYTES];
    uint8_t rho[KYBER_SYMBYTES];
    uint8_t rho2[KYBER_SYMBYTES];
    poly t_hat[KYBER_MAX_K], t_hat2[KYBER_MAX_K];
    poly s_hat[KYBER_MAX_K], s_hat2[KYBER_MAX_K];
    int ok = 1;

    kyber_init_param(&param, mode);
    fill_kyber_bytes(rho, KYBER_SYMBYTES, (uint8_t)(mode & 0xFF));
    for(int i = 0; i < param.k; i++) {
        fill_kyber_poly(t_hat[i], mode + i);
        fill_kyber_poly(s_hat[i], mode + 17 + i);
    }

    kyber_pack_ek(ek, t_hat, rho, param);
    kyber_unpack_ek(t_hat2, rho2, ek, param);
    kyber_pack_dk_pke(dk, s_hat, param);
    kyber_unpack_dk_pke(s_hat2, dk, param);

    if(memcmp(rho, rho2, KYBER_SYMBYTES) != 0) {
        printf("  [FAIL] %s rho mismatch after ek unpack\n", kyber_mode_name(mode));
        ok = 0;
    }
    for(int i = 0; i < param.k; i++) {
        if(!kyber_poly_match(t_hat[i], t_hat2[i])) {
            printf("  [FAIL] %s t_hat[%d] mismatch after ek unpack\n",
                   kyber_mode_name(mode), i);
            ok = 0;
            break;
        }
        if(!kyber_poly_match(s_hat[i], s_hat2[i])) {
            printf("  [FAIL] %s s_hat[%d] mismatch after dk unpack\n",
                   kyber_mode_name(mode), i);
            ok = 0;
            break;
        }
    }
    if(!kyber_check_ek_modulus(ek, param)) {
        printf("  [FAIL] %s valid ek rejected by modulus check\n", kyber_mode_name(mode));
        ok = 0;
    }

    ek[0] = 0xFF;
    ek[1] = 0x0F;
    if(kyber_check_ek_modulus(ek, param)) {
        printf("  [FAIL] %s invalid ek accepted by modulus check\n", kyber_mode_name(mode));
        ok = 0;
    }

    return ok;
}

static int test_kyber_ct_pack_mode(int mode) {
    kyber_param param;
    uint8_t ct[KYBER_MAX_CT_BYTES];
    poly u[KYBER_MAX_K], u2[KYBER_MAX_K];
    poly v, v2;
    poly compressed;
    poly expected;

    kyber_init_param(&param, mode);
    for(int i = 0; i < param.k; i++) {
        fill_kyber_poly(u[i], mode + 31 + i);
    }
    fill_kyber_poly(v, mode + 63);

    kyber_pack_ct(ct, u, v, param);
    kyber_unpack_ct(u2, v2, ct, param);

    for(int i = 0; i < param.k; i++) {
        kyber_poly_compress(compressed, u[i], param.du);
        kyber_poly_decompress(expected, compressed, param.du);
        if(!kyber_poly_match(expected, u2[i])) {
            printf("  [FAIL] %s u[%d] mismatch after ct pack/unpack\n",
                   kyber_mode_name(mode), i);
            return 0;
        }
    }

    kyber_poly_compress(compressed, v, param.dv);
    kyber_poly_decompress(expected, compressed, param.dv);
    if(!kyber_poly_match(expected, v2)) {
        printf("  [FAIL] %s v mismatch after ct pack/unpack\n", kyber_mode_name(mode));
        return 0;
    }
    return 1;
}

static int test_kyber_pack_all_modes(void) {
    int modes[] = {KYBER_MODE_2, KYBER_MODE_3, KYBER_MODE_5};
    int ok = 1;

    printf("\n========================================\n");
    printf("Test: Kyber ek/dk/ct Pack / Unpack\n");
    printf("========================================\n");

    for(int i = 0; i < (int)(sizeof(modes) / sizeof(modes[0])); i++) {
        int mode_ok = test_kyber_ek_dk_pack_mode(modes[i]) &&
                      test_kyber_ct_pack_mode(modes[i]);
        printf("  %-16s %s\n", kyber_mode_name(modes[i]), mode_ok ? "[PASS]" : "[FAIL]");
        ok &= mode_ok;
    }

    printf("%s Kyber pack/unpack mode tests\n", ok ? "[PASS]" : "[FAIL]");
    return ok;
}

int test_kyber_pack_suite(void) {
    int ok = 1;

    printf("\n========================================\n");
    printf("  Kyber Pack Functions Test Suite\n");
    printf("========================================\n");

    ok &= test_kyber_byte_encode_decode();
    ok &= test_kyber_pack_all_modes();
    return ok;
}

int main() {
    int kyber_ok;

    printf("\n");
    printf("========================================\n");
    printf("  Dilithium Pack Functions Test Suite\n");
    printf("========================================\n");
    
    printf("\nParameters:\n");
    printf("  DILI_N = %d\n", DILI_N);
    printf("  DILI_K = %d\n", DILI_K);
    printf("  DILI_L = %d\n", DILI_L);
    printf("  DILI_D = %d\n", DILI_D);
    printf("  DILI_GAMMA1 = %d\n", DILI_GAMMA1);
    printf("  DILI_OMEGA = %d\n", DILI_OMEGA);
    
    test_pack_unpack_z();
    test_pack_unpack_h();
    test_pack_w1();
    test_pack_unpack_t1();
    test_pack_unpack_pk();
    kyber_ok = test_kyber_pack_suite();

    printf("\n========================================\n");
    printf("  All Pack Tests Complete\n");
    printf("========================================\n\n");

    return kyber_ok ? 0 : 1;
}
