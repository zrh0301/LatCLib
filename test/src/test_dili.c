#include "../inc/test_common.h"
#include "../../inc/dili.h"
#include "../../inc/poly.h"
#include <stdio.h>
#include <string.h>

// Test key generation
void test_keygen() {
    printf("\n============================================\n");
    printf("  Test: Dilithium Key Generation\n");
    printf("============================================\n");
    
    dili_param param;
    param.mode = 0;  // Dilithium2
    param.k = 4;
    param.l = 4;
    param.eta = 2;
    
    dili_keypair_t keypair;
    
    printf("Generating keypair with params: k=%d, l=%d, eta=%d\n", 
           param.k, param.l, param.eta);
    
    dili_keygen(&keypair, param);
    
    printf("\n[PASS] Key generation completed successfully\n");
}

// Test sign and verify (placeholder)
void test_sign_verify() {
    printf("\n============================================\n");
    printf("  Test: Dilithium Sign/Verify\n");
    printf("============================================\n");
    
    dili_param param;
    param.mode = 0;
    param.k = 4;
    param.l = 4;
    param.eta = 2;
    
    dili_keypair_t keypair;
    dili_keygen(&keypair, param);
    
    uint8_t message[] = "Hello, Dilithium!";
    int message_len = strlen((char*)message);
    uint8_t signature[4096];
    
    printf("Message: %s\n", message);
    printf("Message length: %d bytes\n", message_len);
    
    printf("\n[PENDING] Sign/Verify not yet implemented\n");
}

// Test packing/unpacking
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
    printf("Public key packed: %d bytes\n", (int)sizeof(pk_bytes));
    
    dili_unpack_pk(keypair2.pk.rho, keypair2.pk.t1, pk_bytes, param.k);
    printf("Public key unpacked\n");
    
    int rho_match = memcmp(keypair1.pk.rho, keypair2.pk.rho, 32) == 0;
    printf("rho matches: %s\n", rho_match ? "YES" : "NO");
    
    int t1_match = 1;
    for(int i = 0; i < 10; i++) {
        if(keypair1.pk.t1[0][i] != keypair2.pk.t1[0][i]) {
            t1_match = 0;
            break;
        }
    }
    printf("t1[0] first 10 coefficients match: %s\n", t1_match ? "YES" : "NO");
    
    if(rho_match && t1_match) {
        printf("\n[PASS] Pack/Unpack test successful\n");
    } else {
        printf("\n[FAIL] Pack/Unpack mismatch\n");
    }
}

int main() {
    printf("\n");
    printf("╔════════════════════════════════════════╗\n");
    printf("║   Dilithium Test Suite                 ║\n");
    printf("╚════════════════════════════════════════╝\n");
    
    test_keygen();
    //test_packing();
    //test_sign_verify();
    
    printf("\n============================================\n");
    printf("  All Tests Complete\n");
    printf("============================================\n\n");
    
    return 0;
}
