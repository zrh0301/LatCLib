#include "../inc/dili.h"
#include "../inc/poly.h"
#include "../inc/rand.h"
#include <string.h>
#include <stdio.h>

// ============================================
//  Helper Functions
// ============================================

// Power2Round: decompose r into (r1, r0) where r = r1*2^d + r0
static void power2round(poly *r1, poly *r0, poly *r, int k, int n, int d) {
    lssi_t mask = (1 << d) - 1;
    for(int j = 0; j < k; j++) {
        for(int i = 0; i < n; i++) {
            r0[j][i] = r[j][i] & mask;
            r1[j][i] = (r[j][i] + (1 << (d - 1))) >> d;
        }
    }
}

// Decompose: split r into high and low bits
static void decompose(lssi_t *r1, lssi_t *r0, const lssi_t *r, int n, int alpha) {
    // TODO: Implement decompose
    // Used for signature generation
}

// HighBits: extract high bits
static void high_bits(lssi_t *r1, const lssi_t *r, int n, int alpha) {
    // TODO: Implement high_bits
}

// LowBits: extract low bits
static void low_bits(lssi_t *r0, const lssi_t *r, int n, int alpha) {
    // TODO: Implement low_bits
}

// MakeHint: create hint bits for rounding correction
static int make_hint(uint8_t *hint, const lssi_t *a, const lssi_t *b, int n, int alpha) {
    // TODO: Implement make_hint
    // Returns number of 1's in hint
    return 0;
}

// UseHint: use hint to recover correct high bits
static void use_hint(lssi_t *r1, const uint8_t *hint, const lssi_t *r, int n, int alpha) {
    // TODO: Implement use_hint
}

// ExpandA: expand seed into matrix A
static void expand_a(poly *matrix_A, int k, int l, const uint8_t *rho, lp_ctx *ctx) {
    uint8_t seed[34];
    uint8_t buf[1024];
    for(int i = 0; i < k; i++) {
        for(int j = 0; j < l; j++) {
            memcpy(seed, rho, 32);
            seed[32] = (uint8_t)j;
            seed[33] = (uint8_t)i;
            shake128(buf, sizeof(buf), seed, 34);
            sample_uniform(matrix_A[i * l + j], ctx->n, buf, sizeof(buf), ctx);
        }
    }
}

// ExpandS: expand seed into secret vectors
static void expand_s(poly *s1, poly *s2, const uint8_t *rho_prime, int k, int l, int eta, lp_ctx *ctx) {
    uint8_t seed[66];
    uint8_t buf[256];

    for(int i = 0; i < l; i++) {
        memcpy(seed, rho_prime, 32);
        seed[32] = (uint8_t)i;
        seed[33] = 0x00;
        shake256(buf, sizeof(buf), seed, 34);
        sample_eta(s1[i], ctx->n, buf, sizeof(buf), ctx, eta);
    }
    for(int i = 0; i < k; i++) {
        memcpy(seed, rho_prime, 32);
        seed[32] = (uint8_t)(l + i);
        seed[33] = 0x00;
        shake256(buf, sizeof(buf), seed, 34);
        sample_eta(s2[i], ctx->n, buf, sizeof(buf), ctx, eta);
    }
}

// Pack/Unpack polynomial t1 with coefficients fitting in 13 bits (DILI_D-1)
static void dili_polyt1_pack(uint8_t *r, const poly p, int n, int d) {
    unsigned int i;
    int bits = d - 1; // 13 bits for t1
    
    // Simplified: pack each coefficient as 2 bytes (can be optimized)
    for(i = 0; i < n; i++) {
        uint32_t t = (1 << bits) - p[i];
        r[i*2] = t & 0xFF;
        r[i*2+1] = (t >> 8) & 0x1F;
    }
}

static void dili_polyt1_unpack(poly p, const uint8_t *a, int n, int d) {
    unsigned int i;
    int bits = d - 1;
    
    for(i = 0; i < n; i++) {
        p[i] = a[i*2] | ((uint32_t)(a[i*2+1] & 0x1F) << 8);
        p[i] = (1 << bits) - p[i];
    }
}

static void calc_t(poly *t, poly *A, int k, int l, lp_ctx *ctx, dili_keypair_t *keypair) {
    poly t_ntt[k];
    for(int i = 0; i < k; i++) {
        POLY_ZERO(t_ntt[i], ctx->n);
        for(int j = 0; j < l; j++) {
            poly tmp;
            // Access A as flat array: A[i * l + j]
            poly_pointwise(tmp, A[i * l + j], keypair->sk.s1[j], ctx);
            POLY_ADD(t_ntt[i], t_ntt[i], tmp, ctx->n);
        }
    }
    for(int i = 0; i < k; i++) {
        poly_invntt(t_ntt[i], ctx);
        POLY_ADD(t[i], t_ntt[i], keypair->sk.s2[i], ctx->n);
    }
}

void dili_pack_pk(uint8_t *pk_bytes, const uint8_t *rho, const poly *t1, int k) {
    unsigned int i;
    
    for(i = 0; i < 32; i++)
        pk_bytes[i] = rho[i];
    pk_bytes += 32;
    
    for(i = 0; i < k; i++) {
        dili_polyt1_pack(pk_bytes, t1[i], DILI_N, DILI_D);
        pk_bytes += (DILI_D - 1) * DILI_N / 8;
    }
}


void dili_unpack_pk(uint8_t *rho, poly *t1, const uint8_t *pk_bytes, int k) {
    unsigned int i;
    
    for(i = 0; i < 32; i++)
        rho[i] = pk_bytes[i];
    pk_bytes += 32;
    
    for(i = 0; i < k; i++) {
        dili_polyt1_unpack(t1[i], pk_bytes, DILI_N, DILI_D);
        pk_bytes += (DILI_D - 1) * DILI_N / 8;
    }
}

// ExpandMask: expand seed into masking vector y
static void expand_mask(lssi_t *y, const uint8_t *seed, int len, int gamma1) {
    // TODO: Sample y with coefficients in [-gamma1+1, gamma1]
}

// Challenge: hash to challenge polynomial c
static void challenge(lssi_t *c, const uint8_t *hash, int tau) {
    // TODO: Generate challenge polynomial with tau non-zero coefficients {-1, 1}
}

// ============================================
//  Main Algorithm Functions
// ============================================



// Key Generation
void dili_keygen(dili_keypair_t* keypair, dili_param param) {
    lp_ctx ctx;
    INIT_NTT_DILI(&ctx);

    uint8_t seed[32];
    rand_bytes(seed, 32);
    
    uint8_t buf[96];
    uint8_t rho[32], rho_prime[32], K[32];

    shake256(buf, 96, seed, 32);
    memcpy(rho, buf, 32);
    memcpy(rho_prime, buf + 32, 32);
    memcpy(K, buf + 64, 32);

    poly matrix_A[param.k][param.l];
    expand_a((poly*)matrix_A, param.k, param.l, rho, &ctx);
    for(int i = 0; i < param.k; i++) {
        for(int j = 0; j < param.l; j++) {
            poly_ntt(matrix_A[i][j], &ctx);
        }
    }
    
    expand_s(keypair->sk.s1, keypair->sk.s2, rho_prime, param.k, param.l, param.eta, &ctx);
    
    for(int i = 0; i < param.l; i++) {
        poly_ntt(keypair->sk.s1[i], &ctx);
    }
    for(int i = 0; i < param.k; i++) {
        poly_ntt(keypair->sk.s2[i], &ctx);
    }

    poly t[param.k];
    calc_t(t, (poly*)matrix_A, param.k, param.l, &ctx, keypair);
    
    power2round(keypair->pk.t1, keypair->sk.t0, t, param.k, ctx.n, DILI_D);

    uint8_t pk_packed[32 + param.k * 416];
    dili_pack_pk(pk_packed, rho, keypair->pk.t1, param.k);
    
    shake256(keypair->sk.tr, 64, pk_packed, 32 + param.k * 416);

    memcpy(keypair->pk.rho, rho, 32);
    memcpy(keypair->sk.rho, rho, 32);
    memcpy(keypair->sk.K, K, 32);
    
    // Print verification
    printf("\n=== Dilithium Key Generation Complete ===\n");
    printf("Parameters: k=%d, l=%d, eta=%d\n", param.k, param.l, param.eta);
    
    printf("\nPublic Key (pk.rho, first 16 bytes): ");
    for(int i = 0; i < 16; i++) printf("%02x", keypair->pk.rho[i]);
    
    printf("\n\nSecret Key (sk.rho, first 16 bytes): ");
    for(int i = 0; i < 16; i++) printf("%02x", keypair->sk.rho[i]);
    
    printf("\n\nSecret Key (sk.K, first 16 bytes): ");
    for(int i = 0; i < 16; i++) printf("%02x", keypair->sk.K[i]);
    
    printf("\n\nSecret Key (sk.tr, first 32 bytes): ");
    for(int i = 0; i < 32; i++) printf("%02x", keypair->sk.tr[i]);
    
    printf("\n\nPublic Key t1[0], first 8 coefficients: ");
    for(int i = 0; i < 8; i++) printf("%d ", keypair->pk.t1[0][i]);
    
    printf("\n\nSecret Key s1[0], first 8 coefficients: ");
    for(int i = 0; i < 8; i++) printf("%d ", keypair->sk.s1[0][i]);
    
    printf("\n\nSecret Key s2[0], first 8 coefficients: ");
    for(int i = 0; i < 8; i++) printf("%d ", keypair->sk.s2[0][i]);
    
    printf("\n\nSecret Key t0[0], first 8 coefficients: ");
    for(int i = 0; i < 8; i++) printf("%d ", keypair->sk.t0[0][i]);
    
    printf("\n\nVerification: rho matches? %s\n", 
           memcmp(keypair->pk.rho, keypair->sk.rho, 32) == 0 ? "YES" : "NO");
    printf("=========================================\n\n");
}

// Sign
void dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk) {
    // TODO: Implement signing
    // 1. Expand A from rho
    // 2. Sample masking vector y
    // 3. Compute w = A*y
    // 4. Extract high bits w1 = HighBits(w, 2*gamma2)
    // 5. Compute challenge c = H(mu || w1) where mu = H(tr || M)
    // 6. Compute z = y + c*s1
    // 7. Check ||z||_inf < gamma1 - beta
    // 8. Compute r0 = LowBits(w - c*s2, 2*gamma2)
    // 9. Check ||r0||_inf < gamma2 - beta
    // 10. Compute hint h = MakeHint(-c*t0, w - c*s2 + c*t0, 2*gamma2)
    // 11. Check ||c*t0||_inf < gamma2 and number of 1's in h <= omega
    // 12. Output sigma = (c, z, h)
}

// Verify
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk) {
    // TODO: Implement verification
    // 1. Parse sigma = (c, z, h)
    // 2. Check ||z||_inf < gamma1 - beta
    // 3. Expand A from rho
    // 4. Compute w' = A*z - c*t1*2^d
    // 5. Use hint to recover w1' = UseHint(h, w', 2*gamma2)
    // 6. Compute c' = H(mu || w1') where mu = H(tr || M)
    // 7. Return (c == c')
    return 0;
}
