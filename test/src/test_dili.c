#include "../inc/test_common.h"
#include "../../inc/dili.h"
#include "../../inc/poly.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>

// Dilithium2 constants for size checks
#define SIG_SIZE_BYTES 2420 // Approx for Dili2 standard
// Our implementation size: 
// 32 (c_tilde) + 4*768 (z) + (80 + 4) (h + lengths) = 32 + 3072 + 84 = 3188 bytes approx
// Let's use a safe buffer size.

// Test key generation
void test_keygen() {
    printf("\n============================================\n");
    printf("  Test: Dilithium Key Generation\n");
    printf("============================================\n");
    
    // Param struct is just for legacy compat if keygen still uses it
    dili_param param;
    param.mode = 0;  // Dilithium2
    param.k = 4;
    param.l = 4;
    param.eta = 2;
    
    dili_keypair_t keypair;
    
    printf("Generating keypair...\n");
    
    dili_keygen(&keypair, param);
    
    // Basic validity check
    int zero_rho = 1;
    for(int i=0; i<32; i++) if(keypair.pk.rho[i] != 0) zero_rho = 0;
    
    if(!zero_rho) {
        printf("PK Rho: ");
        for(int i=0; i<8; i++) printf("%02X", keypair.pk.rho[i]);
        printf("...\n");
        printf("\n[PASS] Key generation completed successfully (Rho is non-zero)\n");
    } else {
        printf("\n[FAIL] Key generation produced zero Rho\n");
    }
}

// Test sign and verify
void test_sign_verify() {
    printf("\n============================================\n");
    printf("  Test: Dilithium Sign & Verify\n");
    printf("============================================\n");
    
    dili_param param;
    param.mode = 0;
    param.k = 4;
    param.l = 4;
    param.eta = 2;
    
    dili_keypair_t keypair;
    printf("Generating keys...\n");
    dili_keygen(&keypair, param);
    
    uint8_t message[] = "Hello, Dilithium! This is a test message for signing.";
    int message_len = (int)strlen((char*)message);
    
    uint8_t signature[4096]; 
    memset(signature, 0, sizeof(signature));
    
    printf("Signing message: \"%s\" (%d bytes)\n", message, message_len);
    printf("Computing signature...\n");
    
    int sign_result = dili_sign(signature, message, message_len, &keypair.sk);
    
    if(sign_result != 0) {
        printf("[FAIL] Signature generation failed!\n");
        return;
    }

    printf("[INFO] Signature generated successfully.\n");
    
    // Check if signature was written
    int non_zero = 0;
    for(int i=0; i<100; i++) {
        if(signature[i] != 0) non_zero = 1;
    }
    
    if(non_zero) {
        printf("Signature (first 32 bytes): ");
        for(int i=0; i<32; i++) printf("%02X", signature[i]);
        printf("\n");
        
        // Test verification
        printf("\nVerifying signature...\n");
        int verify_result = dili_verify(signature, message, message_len, &keypair.pk);
        
        if(verify_result == 0) {
            printf("[PASS] Signature verified successfully!\n");
        } else {
            printf("[FAIL] Signature verification failed (result=%d)\n", verify_result);
        } 
    } else {
        printf("\n[FAIL] Signature is all zeros\n");
    }
}

// Test packing/unpacking (Keep existing logic)
void test_packing() {
    printf("\n============================================\n");
    printf("  Test: Pack/Unpack Public Key\n");
    printf("============================================\n");
    
    dili_param param;
    param.mode = 0;
    param.k = 4;
    param.l = 4;
    param.eta = 2;
    
    dili_keypair_t keypair1, keypair2;
    dili_keygen(&keypair1, param);
    
    uint8_t pk_bytes[32 + 4 * 416];
    dili_pack_pk(pk_bytes, keypair1.pk.rho, keypair1.pk.t1, param.k);
    
    dili_unpack_pk(keypair2.pk.rho, keypair2.pk.t1, pk_bytes, param.k);
    
    int rho_match = memcmp(keypair1.pk.rho, keypair2.pk.rho, 32) == 0;
    int t1_match = 1;
    for(int i = 0; i < 10; i++) {
        // Just checking first poly
        for(int j=0; j<256; j++) {
             if(keypair1.pk.t1[0][j] != keypair2.pk.t1[0][j]) {
                t1_match = 0;
                break;
            }
        }
    }
    
    if(rho_match && t1_match) {
        printf("[PASS] Pack/Unpack test successful\n");
    } else {
        printf("[FAIL] Pack/Unpack mismatch\n");
        printf("Rho match: %d, T1 match: %d\n", rho_match, t1_match);
    }
}

int main() {
    printf("\n");
    printf("╔════════════════════════════════════════╗\n");
    printf("║   Dilithium Test Suite                 ║\n");
    printf("╚════════════════════════════════════════╝\n");
    
    test_keygen();
    // test_packing(); 
    test_sign_verify(); // Test both sign and verify
    
    printf("\n============================================\n");
    printf("  All Tests Complete\n");
    printf("============================================\n\n");
    
    return 0;
}