#include <stdio.h>

#include "poly.h"

static void print_zetas(const char *name, const lp_ctx *ctx, int used) {
    printf("%s dynamic zetas\n", name);
    printf("q=%d, n=%d, layer=%d, used=%d\n", ctx->q, ctx->n, ctx->layer, used);

    for(int i = 0; i < used; i++) {
        if(i % 8 == 0) {
            printf("  ");
        }

        printf("%8d", ctx->zetas[i]);

        if(i % 8 == 7 || i == used - 1) {
            printf("\n");
        } else {
            printf(" ");
        }
    }
    printf("\n");
}

int main() {
    lp_ctx ctx;

    INIT_NTT_DILI(&ctx);
    print_zetas("Dilithium", &ctx, 256);

    INIT_NTT_KYBER(&ctx);
    print_zetas("Kyber", &ctx, 128);

    return 0;
}
