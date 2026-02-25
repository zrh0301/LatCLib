#include "../inc/test_common.h"
#include "../../inc/dili.h"
#include "../../inc/poly.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define SIG_SIZE_BYTES 2420 // Approx for Dili2 standard

// 初始化不同模式的参数
void init_dili_param(dili_param *param, int mode) {
    param->mode = mode;
    
    switch(mode) {
        case 0: // Dilithium2
            param->k = 4;
            param->l = 4;
            param->eta = 2;
            param->tau = 39;
            param->beta = 78;
            param->omega = 80;
            param->gamma1 = (1 << 17);      // 2^17 = 131072
            param->gamma2 = (DILI_Q-1)/88;  // 95232
            param->w1_packed_bytes = 192;   // 6-bit encoding: (256/4)*3 = 192
            break;
        case 1: // Dilithium3
            param->k = 6;
            param->l = 5;
            param->eta = 4;
            param->tau = 49;
            param->beta = 196;
            param->omega = 55;
            param->gamma1 = (1 << 19);      // 2^19 = 524288
            param->gamma2 = (DILI_Q-1)/32;  // 261888
            param->w1_packed_bytes = 128;   // 4-bit encoding: 256/2 = 128
            break;
        case 2: // Dilithium5
            param->k = 8;
            param->l = 7;
            param->eta = 2;
            param->tau = 60;
            param->beta = 120;
            param->omega = 75;
            param->gamma1 = (1 << 19);      // 2^19 = 524288
            param->gamma2 = (DILI_Q-1)/32;  // 261888
            param->w1_packed_bytes = 128;   // 4-bit encoding: 256/2 = 128
            break;
        default:
            printf("[ERROR] Unknown mode %d, defaulting to Dilithium2\n", mode);
            param->k = 4;
            param->l = 4;
            param->eta = 2;
            param->tau = 39;
            param->beta = 78;
            param->omega = 80;
            param->gamma1 = (1 << 17);
            param->gamma2 = (DILI_Q-1)/88;
            param->mode = 0;
    }
}

// 获取模式名称
const char* get_mode_name(int mode) {
    switch(mode) {
        case 0: return "Dilithium2";
        case 1: return "Dilithium3";
        case 2: return "Dilithium5";
        default: return "Unknown";
    }
}

// Test key generation with specified mode
void test_keygen_mode(int mode) {
    printf("\n============================================\n");
    printf("  Test: %s Key Generation\n", get_mode_name(mode));
    printf("============================================\n");
    
    dili_param param;
    init_dili_param(&param, mode);
    
    printf("Parameters: k=%d, l=%d, eta=%d\n", param.k, param.l, param.eta);
    
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
        printf("\n[PASS] %s key generation completed successfully (Rho is non-zero)\n", get_mode_name(mode));
    } else {
        printf("\n[FAIL] %s key generation produced zero Rho\n", get_mode_name(mode));
    }
}

// Test key generation for all modes
void test_keygen() {
    printf("\n============================================\n");
    printf("  Test: All Dilithium Modes Key Generation\n");
    printf("============================================\n");
    
    for(int mode = 0; mode <= 2; mode++) {
        test_keygen_mode(mode);
    }
}

// Test sign and verify with specified mode
void test_sign_verify_mode(int mode) {
    printf("\n============================================\n");
    printf("  Test: %s Sign & Verify\n", get_mode_name(mode));
    printf("============================================\n");
    
    dili_param param;
    init_dili_param(&param, mode);
    
    printf("Parameters: k=%d, l=%d, eta=%d\n", param.k, param.l, param.eta);
    printf("Parameters: k=%d, l=%d, eta=%d\n", param.k, param.l, param.eta);
    
    dili_keypair_t keypair;
    printf("Generating keys...\n");
    dili_keygen(&keypair, param);
    
    uint8_t message[] = "Hello, Dilithium! This is a test message for signing.";
    int message_len = (int)strlen((char*)message);
    
    uint8_t signature[8192]; 
    memset(signature, 0, sizeof(signature));
    
    printf("Signing message: \"%s\" (%d bytes)\n", message, message_len);
    printf("Computing signature...\n");
    
    int sign_result = dili_sign(signature, message, message_len, &keypair.sk, param);
    
    if(sign_result != 0) {
        printf("[FAIL] %s signature generation failed!\n", get_mode_name(mode));
        return;
    }

    printf("[INFO] %s signature generated successfully.\n", get_mode_name(mode));
    
    // Check if signature was written
    int non_zero = 0;
    for(int i=0; i<100; i++) {
        if(signature[i] != 0) non_zero = 1;
    }
    printf("Signature (first 32 bytes): ");
    for(int i=0; i<32; i++) printf("%02X", signature[i]);
    printf("\n");

    if(non_zero) {

        printf("\nVerifying signature with correct message...\n");
        // Test verification
        int verify_result = dili_verify(signature, message, message_len, &keypair.pk, param);
        
        if(verify_result == 0) {
            printf("[PASS] %s signature verified successfully!\n", get_mode_name(mode));
        } else {
            printf("[FAIL] %s signature verification failed (result=%d)\n", get_mode_name(mode), verify_result);
        } 
    } else {
        printf("\n[FAIL] Signature is all zeros\n");
    }

    char message2[] = "hello, Dilithium! This is a test message for signing.";
    printf("\nVerifying signature with wrong message (should fail)...\n");
    if(non_zero) {

        // Test verification with wrong message
        int verify_result = dili_verify(signature, (uint8_t*)message2, (int)strlen(message2), &keypair.pk, param);
        
        if(verify_result != 0) {
            printf("[PASS] %s correctly rejected wrong message (result=%d)\n", get_mode_name(mode), verify_result);
        } else {
            printf("[FAIL] %s incorrectly verified wrong message!\n", get_mode_name(mode));
        } 
    } else {
        printf("\n[FAIL] Signature is all zeros\n");
    }
}

// Test sign and verify for all modes
void test_sign_verify() {
    printf("\n============================================\n");
    printf("  Test: All Dilithium Modes Sign & Verify\n");
    printf("============================================\n");
    
    for(int mode = 0; mode <= 2; mode++) {
        test_sign_verify_mode(mode);
    }
}

// Test packing/unpacking with specified mode
void test_packing_mode(int mode) {
    printf("\n============================================\n");
    printf("  Test: %s Pack/Unpack Public Key\n", get_mode_name(mode));
    printf("============================================\n");
    
    dili_param param;
    init_dili_param(&param, mode);
    
    printf("Parameters: k=%d, l=%d, eta=%d\n", param.k, param.l, param.eta);
    
    dili_keypair_t keypair1, keypair2;
    dili_keygen(&keypair1, param);
    
    uint8_t pk_bytes[32 + 8 * 416]; // Max size for Dilithium5 (k=8)
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
        printf("[PASS] %s Pack/Unpack test successful\n", get_mode_name(mode));
    } else {
        printf("[FAIL] %s Pack/Unpack mismatch\n", get_mode_name(mode));
        printf("Rho match: %d, T1 match: %d\n", rho_match, t1_match);
    }
}

// Test packing for all modes (Keep existing logic)
void test_packing() {
    printf("\n============================================\n");
    printf("  Test: All Dilithium Modes Pack/Unpack\n");
    printf("============================================\n");
    
    for(int mode = 0; mode <= 2; mode++) {
        test_packing_mode(mode);
    }
}

int main() {
    printf("\n");
    printf("**************************************************\n");
    printf("*  Dilithium Multi-Mode Test Suite              *\n");
    printf("*  Testing: Dilithium2, Dilithium3, Dilithium5  *\n");
    printf("**************************************************\n");
    
    test_keygen();
    test_sign_verify();
    //test_packing();
    
    printf("\n============================================\n");
    printf("  All Tests Complete\n");
    printf("============================================\n\n");
    
    return 0;
}