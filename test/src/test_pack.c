#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "poly.h"
#include "dili.h"
#include "pack.h"
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

void test_pack_unpack_z() {
    printf("\n========================================\n");
    printf("Test: pack_z / unpack_z\n");
    printf("========================================\n");
    
    poly z_original, z_unpacked;
    uint8_t packed[3 * DILI_N];
    
    // Initialize test data
    for(int i = 0; i < DILI_N; i++) {
        z_original[i] = (i * 137 - 50000) % (2 * DILI_GAMMA1);
        if(z_original[i] > DILI_GAMMA1) z_original[i] -= 2 * DILI_GAMMA1;
    }
    
    printf("\nOriginal z data:\n");
    print_poly("  z_original", z_original, 8);
    printf("  z_original[0] = %d\n", z_original[0]);
    printf("  z_original[1] = %d\n", z_original[1]);
    printf("  z_original[255] = %d\n", z_original[255]);
    
    // Pack
    pack_z(packed, z_original);
    printf("\nPacked bytes:\n");
    print_bytes("  packed", packed, 3 * DILI_N);
    printf("  packed[0] = 0x%02x\n", packed[0]);
    printf("  packed[1] = 0x%02x\n", packed[1]);
    printf("  packed[2] = 0x%02x\n", packed[2]);
    
    // Unpack
    unpack_z(z_unpacked, packed);
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
        printf("\n[PASS] pack_z/unpack_z test passed!\n");
    } else {
        printf("\n[FAIL] pack_z/unpack_z has %d errors!\n", errors);
    }
}

void test_pack_unpack_h() {
    printf("\n========================================\n");
    printf("Test: pack_h / unpack_h\n");
    printf("========================================\n");
    
    const int K = DILI_K;
    uint8_t h_original[K][DILI_N];
    uint8_t h_unpacked[K][DILI_N];
    uint8_t packed[DILI_OMEGA + K];
    
    // Initialize test data - sparse hint vectors
    for(int i = 0; i < K; i++) {
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
    
    printf("\nOriginal hint data:\n");
    int total_ones = 0;
    for(int i = 0; i < K; i++) {
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
    printf("  Total: %d ones (OMEGA=%d)\n", total_ones, DILI_OMEGA);
    
    // Pack
    int pack_result = pack_h(packed, h_original, K);
    printf("\nPack result: %d\n", pack_result);
    print_bytes("  packed", packed, DILI_OMEGA + K);
    
    // Unpack
    unpack_h(h_unpacked, packed, K);
    printf("\nUnpacked hint data:\n");
    for(int i = 0; i < K; i++) {
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
    for(int i = 0; i < K; i++) {
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
    
    if(errors == 0) {
        printf("\n[PASS] pack_h/unpack_h test passed!\n");
    } else {
        printf("\n[FAIL] pack_h/unpack_h has %d errors!\n", errors);
    }
}

void test_pack_w1() {
    printf("\n========================================\n");
    printf("Test: pack_w1\n");
    printf("========================================\n");
    
    poly w1;
    uint8_t packed[DILI_N];
    
    // Initialize test data - w1 values should be in [0, 43]
    for(int i = 0; i < DILI_N; i++) {
        w1[i] = i % 44;
    }
    
    printf("\nOriginal w1 data:\n");
    print_poly("  w1", w1, 16);
    printf("  w1[0] = %d\n", w1[0]);
    printf("  w1[1] = %d\n", w1[1]);
    printf("  w1[2] = %d\n", w1[2]);
    printf("  w1[3] = %d\n", w1[3]);
    
    // Pack
    pack_w1(packed, w1);
    printf("\nPacked bytes:\n");
    print_bytes("  packed", packed, DILI_N);
    printf("  packed[0] = 0x%02x (should contain w1[0] and w1[1] low bits)\n", packed[0]);
    printf("  packed[1] = 0x%02x\n", packed[1]);
    printf("  packed[2] = 0x%02x\n", packed[2]);
    
    // Manual verification of the first few packed values
    // pack_w1 uses 6 bits to store each coefficient
    uint32_t t0 = w1[0] & 0x3F;
    uint32_t t1 = w1[1] & 0x3F;
    uint32_t expected_byte0 = t0 | (t1 << 6);
    printf("\nManual verification:\n");
    printf("  w1[0] & 0x3F = 0x%02x\n", t0);
    printf("  w1[1] & 0x3F = 0x%02x\n", t1);
    printf("  expected packed[0] = 0x%02x\n", (uint8_t)expected_byte0);
    
    if((uint8_t)expected_byte0 == packed[0]) {
        printf("\n[PASS] pack_w1 first bytes correct!\n");
    } else {
        printf("\n[FAIL] pack_w1 packing incorrect!\n");
    }
}

void test_pack_unpack_t1() {
    printf("\n========================================\n");
    printf("Test: dili_polyt1_pack / dili_polyt1_unpack\n");
    printf("========================================\n");
    
    poly t1_original, t1_unpacked;
    // Each t1 coefficient uses (D-1) bits, simplified implementation uses 2 bytes/coeff
    // 256 coefficients * 2 bytes = 512 bytes
    uint8_t packed[512];
    
    // Initialize test data - t1 values should be in [0, 2^(D-1)), for D=13 that's [0, 4096)
    for(int i = 0; i < DILI_N; i++) {
        t1_original[i] = (i * 17) % (1 << (DILI_D - 1)); // Range [0, 4096)
    }
    
    printf("\nOriginal t1 data:\n");
    print_poly("  t1_original", t1_original, 8);
    printf("  t1_original[0] = %d\n", t1_original[0]);
    printf("  t1_original[1] = %d\n", t1_original[1]);
    printf("  t1_original[255] = %d\n", t1_original[255]);
    
    // Pack
    dili_polyt1_pack(packed, t1_original, DILI_N, DILI_D);
    printf("\nPacked bytes:\n");
    print_bytes("  packed", packed, 512);
    
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
        int32_t expected = (i * 17) % (1 << (DILI_D - 1)); // [0, 4096)
        
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
    
    // Public key size: 32 bytes (rho) + K * (D-1) * N / 8 bytes (t1)
    // For simplified implementation: 32 + K * 512 = 32 + 4*512 = 2080 bytes
    uint8_t pk_bytes[2080];
    
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
            t1_original[k][i] = ((k * 256 + i) * 19) % (1 << (DILI_D - 1));
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

int main() {
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
    
    printf("\n========================================\n");
    printf("  All Pack Tests Complete\n");
    printf("========================================\n\n");
    
    return 0;
}
