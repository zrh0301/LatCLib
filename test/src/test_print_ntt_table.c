#include <stdio.h>
#include "../../inc/poly.h"
#include "../../inc/dili.h"

int main() {
    lp_ctx Q;
    
    // Dilithium NTT parameters
    init_ntt(&Q, 256, 1753, NTT_FULL_MODE);
    printf("dili: q=%d, n=%d, root=1753, n_inv=%d\n\n", Q.q, Q.n, Q.n_inv);
    
    printf("static const lssi_t DILI_NTT_W[%d] = {\n", Q.n);
    for(int i = 0; i < Q.n; i++) {
        if(i % 8 == 0) printf("    ");
        printf("%d", Q.w[i]);
        if(i < Q.n - 1) printf(", ");
        if(i % 8 == 7) printf("\n");
    }
    printf("};\n\n");
    
    printf("static const lssi_t DILI_NTT_W_INV[%d] = {\n", Q.n);
    for(int i = 0; i < Q.n; i++) {
        if(i % 8 == 0) printf("    ");
        printf("%d", Q.w_inv[i]);
        if(i < Q.n - 1) printf(", ");
        if(i % 8 == 7) printf("\n");
    }
    printf("};\n\n");
    
    printf("static const lssi_t DILI_NTT_N_INV = %d;\n\n", Q.n_inv);
    
    printf("static const int DILI_NTT_REV_TABLE[%d] = {\n", Q.n);
    for(int i = 0; i < Q.n; i++) {
        if(i % 16 == 0) printf("    ");
        printf("%3d", Q.rev_table[i]);
        if(i < Q.n - 1) printf(", ");
        if(i % 16 == 15) printf("\n");
    }
    printf("};\n\n");
    
    printf("========================================\n\n");
    
    // Kyber NTT parameters
    init_ntt(&Q, 256, 17, NTT_INCOMPLETE_MODE);
    printf("kyber: q=%d, n=%d, root=17, n_inv=%d\n\n", Q.q, Q.n, Q.n_inv);
    
    printf("static const lssi_t KYBER_NTT_W[%d] = {\n", Q.n);
    for(int i = 0; i < Q.n; i++) {
        if(i % 8 == 0) printf("    ");
        printf("%d", Q.w[i]);
        if(i < Q.n - 1) printf(", ");
        if(i % 8 == 7) printf("\n");
    }
    printf("};\n\n");
    
    printf("static const lssi_t KYBER_NTT_W_INV[%d] = {\n", Q.n);
    for(int i = 0; i < Q.n; i++) {
        if(i % 8 == 0) printf("    ");
        printf("%d", Q.w_inv[i]);
        if(i < Q.n - 1) printf(", ");
        if(i % 8 == 7) printf("\n");
    }
    printf("};\n\n");
    
    printf("static const lssi_t KYBER_NTT_N_INV = %d;\n\n", Q.n_inv);
    
    printf("static const int KYBER_NTT_REV_TABLE[%d] = {\n", Q.n);
    for(int i = 0; i < Q.n; i++) {
        if(i % 16 == 0) printf("    ");
        printf("%3d", Q.rev_table[i]);
        if(i < Q.n - 1) printf(", ");
        if(i % 16 == 15) printf("\n");
    }
    printf("};\n");
    
    return 0;
}
