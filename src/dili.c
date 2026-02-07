#include "../inc/dili.h"
#include "../inc/poly.h"
#include "../inc/polyvec.h"
#include "../inc/pack.h"
#include "../inc/rand.h"
#include <string.h>
#include <stdio.h>

// ============================================
//  Helper Functions
// ============================================

#define CMOD_Q(a) ({ \
    int32_t _r = (int32_t)(a)%DILI_Q; \
    _r += (_r>>31)&DILI_Q; \
    _r -= (((DILI_Q>>1)-_r)>>31)&DILI_Q; \
    _r; \
})

#define SGN_MASK(x) ((int32_t)(x)>>31)

#define ABS_VAL(x) ({ \
    int32_t _m = SGN_MASK(x); \
    ((x)^_m)-_m; \
})

#define CMOV(cond, a, b) ({ \
    int32_t _c = -(int32_t)(cond); \
    ((a)&_c) | ((b)&(~c)); \
}) // return cond?a:b;

#define REDUCE32(x) ({\
    int32_t _r = (x)%DILI_Q; \
    _r += (DILI_Q & SGN_MASK(_r)); \
    _r; \
}) // int32 -> [0,q-1]

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
    int i;
    lssi_t r0_val, r1_val, a, mask;
    for(i=0;i<n;++i){
        a = REDUCE32(r[i]);
        r0_val = a%alpha;
        mask = SGN_MASK((alpha>>1) - r0_val);
        r0_val -= (alpha&mask);
        r1_val = (a-r0_val)/alpha;

        r0[i] = r0_val;
        r1[i] = r1_val; 
    }
}

// HighBits: extract high bits
static void high_bits(lssi_t *r1, const lssi_t *r, int n, int alpha) {
    int i;
    lssi_t a, r0_val, r1_val, mask;
    for(i=0;i<n;++i){
        a = REDUCE32(r[i]);
        r0_val = a % alpha;
        mask = SGN_MASK((alpha>>1) - r0_val);
        r0_val -= (alpha & mask);
        r1_val = (a - r0_val) / alpha;
        r1[i] = r1_val;
    }
}

static int check_norm(const poly p, int bound){
    int fail, i;
    lssi_t val, abs_val, diff;
    fail = 0;
    for(i=0;i<DILI_N;++i){
        val = CMOD_Q(p[i]);
        abs_val = ABS_VAL(val);
        diff = (bound-1)-abs_val;
        fail |= SGN_MASK(diff);
    }
    return (fail&1);
}

// LowBits: extract low bits
static void low_bits(lssi_t *r0, const lssi_t *r, int n, int alpha) {
    int i;
    lssi_t a, r0_val, mask;
    for(i=0;i<n;++i){
        a = REDUCE32(r[i]);
        r0_val = a % alpha;
        mask = SGN_MASK((alpha>>1) - r0_val);
        r0_val -= (alpha & mask);
        r0[i] = r0_val;
    }
}

// MakeHint: create hint bits for rounding correction
static int make_hint(uint8_t *hint, const lssi_t *a, const lssi_t *b, int n, int alpha) {
    // TODO: Implement make_hint
    // Returns number of 1's in hint
    int weight, i;
    lssi_t b_std, r0_b, r1_b, sum_std, r0_sum, r1_sum, diff, is_diff;
    weight = 0;
    for(i=0;i<n;++i){
        b_std = REDUCE32(b[i]);
        r0_b = b_std%alpha;
        r0_b -= (alpha&SGN_MASK((alpha>>1)-r0_b));
        r1_b = (b_std-r0_b)/alpha;

        sum_std = REDUCE32(a[i]+b[i]);
        r0_sum = sum_std%alpha;
        r0_sum -= (alpha&SGN_MASK((alpha>>1)-r0_sum));
        r1_sum = (sum_std-r0_sum)/alpha;
        // r1_b != r1_sum => hint = 1
        diff = r1_b^r1_sum;
        is_diff = ((diff|-diff)>>31)&1;
        hint[i] = (uint8_t)is_diff;
        weight += is_diff;
    }

    return weight;
}

// UseHint: use hint to recover correct high bits
static void use_hint(lssi_t *r1, const uint8_t *hint, const lssi_t *r, int n, int alpha) {
    int i;
    lssi_t a, r0, r1_val;
    int m = (DILI_Q - 1) / alpha;
    
    for(i = 0; i < n; ++i) {
        a = REDUCE32(r[i]);
        r0 = a % alpha;
        
        // Center r0 around 0: if r0 > alpha/2, then r0 = r0 - alpha
        if(r0 > (alpha >> 1)) {
            r0 -= alpha;
        }
        
        r1_val = (a - r0) / alpha;
        
        if(hint[i] == 1) {
            // Apply hint correction
            if(r0 > 0) {
                r1_val = (r1_val + 1) % (m + 1);
            } else {
                r1_val = (r1_val == 0) ? m : r1_val - 1;
            }
        }
        
        r1[i] = r1_val;
    }
}

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

static void calc_t(poly *t, poly *A_mat, int k, int l, lp_ctx *ctx, poly *s1_ntt, poly *s2) {
    // t = A*s1 + s2 (计算在NTT域)
    polyvec_matrix_pointwise_montgomery(t, A_mat, k, s1_ntt, l, ctx);
    polyvec_invntt_tomont(t, k, ctx);
    polyvec_add(t, t, s2, k);
}

// ExpandMask: expand seed into masking vector y
static void expand_mask(lssi_t *y, const uint8_t *seed, int seed_len, int len, int gamma1) {
    // TODO: Sample y with coefficients in [-gamma1+1, gamma1]
    uint8_t buf[1024];
    
    KECCAK1600_CTX ctx;
    shake256_init(&ctx);
    sha3_update(&ctx, seed, seed_len);
    shake_xof_final(&ctx);
    
    int ctr = 0;
    while(ctr < len){
        shake_squeeze(&ctx, buf, sizeof(buf));
        
        int i;
        for(i=0; i<sizeof(buf)-3 && ctr<len;i+=3){
            uint32_t val = buf[i];
            val |= (uint32_t)buf[i+1] << 8;
            val |= (uint32_t)buf[i+2] << 16;
            val &= 0xFFFFF;
            
            uint32_t bound = 2 * gamma1 - 1;
            
            // diff = val - bound. if val < bound, diff < 0 (high bit 1)
            int32_t diff = (int32_t)val - (int32_t)bound;
            int32_t mask = SGN_MASK(diff); // -1 if valid, 0 if invalid
            
            if(mask){
                y[ctr++] = (lssi_t)(gamma1 - 1) - (lssi_t)val;
            }
        }
    }    
}

// Challenge: hash to challenge polynomial c
static void challenge(lssi_t *c, const uint8_t *hash, int tau) {
    // TODO: Generate challenge polynomial with tau non-zero coefficients {-1, 1}
        uint8_t buf[256];
    uint64_t signs = 0;
    int i, b;
    shake256(buf, 256, hash, 32);

    for(i = 0; i < 8; ++i)
        signs |= (uint64_t)buf[i] << (8 * i);

    for(i = 0; i < DILI_N; ++i) c[i] = 0;

    int pos = 8;
    for(i = DILI_N - tau; i < DILI_N; ++i) {
        // Rejection sampling for position b <= i
        do {
            if(pos >= 256) {
                shake256(buf, 256, hash, 32);
                pos = 0;
            }
            b = buf[pos++];
        } while(b > i);

        c[i] = c[b];
        c[b] = (signs & 1) ? -1 : 1;
        signs >>= 1;
    }
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

    poly matrix_A[DILI_K * DILI_L];
    polyvec_matrix_expand(matrix_A, param.k, param.l, rho, &ctx);
    
    expand_s(keypair->sk.s1, keypair->sk.s2, rho_prime, param.k, param.l, param.eta, &ctx);
    
    poly s1_ntt[DILI_L];
    polyvec_copy(s1_ntt, keypair->sk.s1, param.l);
    polyvec_ntt(s1_ntt, param.l, &ctx);

    poly t[DILI_K];
    poly s2[DILI_K];
    polyvec_copy(s2, keypair->sk.s2, param.k);
    calc_t(t, matrix_A, param.k, param.l, &ctx, s1_ntt, s2);
    
    poly t_arr[param.k];
    polyvec_copy(t_arr, t, param.k);
    
    power2round(keypair->pk.t1, keypair->sk.t0, t_arr, param.k, ctx.n, DILI_D);

    uint8_t pk_packed[32 + param.k * 416];
    dili_pack_pk(pk_packed, rho, keypair->pk.t1, param.k);
    
    shake256(keypair->sk.tr, 64, pk_packed, 32 + param.k * 416);

    memcpy(keypair->pk.rho, rho, 32);
    memcpy(keypair->sk.rho, rho, 32);
    memcpy(keypair->sk.K, K, 32);
}

// Sign
void dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk) {
    lp_ctx ctx;
    INIT_NTT_DILI(&ctx);
    
    const int K = DILI_K;
    const int L = DILI_L;
    
    // Polys
    poly matrix_A[K * L];  // 使用一维数组
    poly y[L];
    poly w[K], w1[K];
    poly cp; // Challenge polynomial
    poly z[L];
    poly cs1[L], cs2[K], ct0[K]; // Temps for checks
    uint8_t h[K][DILI_N]; // Hints (binary array)

    // Buffers
    uint8_t mu[64];
    uint8_t c_tilde[32];
    uint8_t c_tilde_temp_buf[64 + K * DILI_N]; // Seed + w1 packed
    uint8_t rho_prime_nonce[64]; // For sampling y
    uint16_t nonce = 0;
    int retry = 1;
    int loop = 0;
#ifdef DEBUG
    printf("expanding A\n");
#endif

    //  Expand Matrix A
    polyvec_matrix_expand(matrix_A, K, L, sk->rho, &ctx);

    // Compute mu = CRH(tr || M)
    KECCAK1600_CTX shake_ctx;
    shake256_init(&shake_ctx);
    sha3_update(&shake_ctx, sk->tr, 64);
    sha3_update(&shake_ctx, message, message_len);
    shake_xof_final(&shake_ctx);
    shake_squeeze(&shake_ctx, mu, 64);
#ifdef DEBUG
    printf("entering rejection loop\n");
#endif
    // Rejection Sampling Loop
    do{
        ++loop;
#ifdef DEBUG
        printf("attempt %d\n", loop);
#endif
        // sample y
        // use sk->K and mu combined with nonce.
        memcpy(rho_prime_nonce, sk->K, 32);
        memcpy(rho_prime_nonce + 32, mu, 32); // Using mu as diversifier
        
        uint8_t y_seed_base[64];
        shake256(y_seed_base, 64, rho_prime_nonce, 64);

        for(int i=0; i<L; i++) {
            // Append nonce to seed
            uint8_t seed_with_nonce[66];
            memcpy(seed_with_nonce, y_seed_base, 64);
            seed_with_nonce[64] = nonce & 0xFF;
            seed_with_nonce[65] = (nonce >> 8) & 0xFF;
            nonce++;
            
            expand_mask(y[i], seed_with_nonce, 66, DILI_N, DILI_GAMMA1);
        }

        // w = A * y
        poly y_ntt[L];
        polyvec_copy(y_ntt, y, L);
        polyvec_ntt(y_ntt, L, &ctx);

        polyvec_matrix_pointwise_montgomery(w, matrix_A, K, y_ntt, L, &ctx);
        polyvec_invntt_tomont(w, K, &ctx);
        
        // Extract HighBits
        for(int i=0; i<K; i++) {
            high_bits(w1[i], w[i], DILI_N, DILI_ALPHA);
        }

        // c_tilde = H(mu || w1)
        memcpy(c_tilde_temp_buf, mu, 64);
        int offset = 64;
        for(int i=0; i<K; i++) {
            pack_w1(c_tilde_temp_buf + offset, w1[i]);
            offset += DILI_N; // 256 bytes per poly in simple packing
        }
        shake256(c_tilde, 32, c_tilde_temp_buf, offset);

        // Compute Challenge Polynomial c
        challenge(cp, c_tilde, DILI_TAU);
        poly cp_ntt;
        POLY_COPY(cp_ntt, cp, DILI_N);
        poly_ntt(cp_ntt, &ctx);

        // z = y + c * s1
        // 先批量计算 cs1 = c * s1
        poly s1_ntt[L];
        polyvec_copy(s1_ntt, sk->s1, L);
        polyvec_ntt(s1_ntt, L, &ctx);

        for(int i=0; i<L; i++) {
            poly_pointwise(cs1[i], cp_ntt, s1_ntt[i], &ctx);
        }

        polyvec_invntt_tomont(cs1, L, &ctx);
        
#ifdef DEBUG
        {
            lssi_t max_cs1 = 0;
            for(int k=0;k<DILI_N;++k){
                lssi_t val = CMOD_Q(cs1[0][k]);
                if(abs(val) > max_cs1) max_cs1 = abs(val);
            }
            printf("max coeff: %d\n", max_cs1);
            printf("raw cs1: %d %d %d %d\n",cs1[0][0], cs1[0][1], cs1[0][2], cs1[0][3]);
        }
#endif
        
        // z = y + cs1
        polyvec_add(z, y, cs1, L);

#ifdef DEBUG
        {
            lssi_t y_max = 0;
            for(int k=0;k<DILI_N;++k){
                if(abs(y[0][k] > y_max)) y_max = abs(y[0][k]);
            }
            printf("max y: %d\n", y_max);
        }
#endif
        
        // Check ||z|| < gamma1 - beta
        int z_check_fail = 0;
        for(int i=0; i<L; i++) {
            if(check_norm(z[i], DILI_GAMMA1 - DILI_BETA)) {
                z_check_fail = 1;
#ifdef DEBUG  

void poly_print(const char* name, const poly p, int n) {
    printf("%s: ", name);
    for(int i=0; i<n; i++) {
        printf("%d ", p[i]);
    }
    printf("\n");
}

                printf("reject because of z(%d) check failed at %d\n", z[i][0], i);
                //打印y cs1
                poly_print("y", y[i], DILI_N);
                poly_print("cs1", cs1[i], DILI_N);
                poly_print("z", z[i], DILI_N);
#endif
                break;
            }
        }
        if(z_check_fail) continue; // Restart rejection loop

        // 批量计算 cs2 = c * s2 和 ct0 = c * t0
        for(int i=0; i<K; i++) {
            poly_pointwise(cs2[i], cp_ntt, sk->s2[i], &ctx);
            poly_pointwise(ct0[i], cp_ntt, sk->t0[i], &ctx);
        }
        polyvec_invntt_tomont(cs2, K, &ctx);
        polyvec_invntt_tomont(ct0, K, &ctx);
        
        // Check LowBits(w - c*s2)
        int r0_check_fail = 0;
        for(int i=0; i<K; i++) {
            poly diff;
            POLY_SUB(diff, w[i], cs2[i], DILI_N);
            
            poly r0;
            low_bits(r0, diff, DILI_N, DILI_ALPHA);
            
            // Check ||r0|| < gamma2 - beta
            if(check_norm(r0, DILI_GAMMA2 - DILI_BETA)) {
                r0_check_fail = 1;
#ifdef DEBUG
                printf("reject because of r0 check fail at %d\n", i);
#endif                
                break;
            }
        }
        if(r0_check_fail) continue;

        // Check ||c * t0||
        int ct0_check_fail = 0;
        for(int i=0; i<K; i++) {
            // Check ||ct0|| < gamma2
            if(check_norm(ct0[i], DILI_GAMMA2)) {
                ct0_check_fail = 1;
                break;
            }
        }
        if(ct0_check_fail) continue;

        // 3.8 Make Hints
        // h = MakeHint(-ct0, w - cs2 + ct0)
        // Note: MakeHint(z, r) -> checks HighBits(r) vs HighBits(r+z)
        // Here: r = w - cs2 + ct0. z = -ct0. 
        // Then r + z = w - cs2.
        // We compare HighBits(w - cs2 + ct0) with HighBits(w - cs2)
        // Wait, spec says: h = MakeHint(-ct0, w - cs2 + ct0)
        // Our MakeHint(a, b) does HighBits(b) vs HighBits(b+a).
        // Let a = -ct0, b = w - cs2 + ct0.
        // b + a = w - cs2.
        // Correct.
        
        int hint_ones = 0;
        for(int i=0; i<K; i++) {
            poly w_cs2_ct0;
            poly diff;
            POLY_SUB(diff, w[i], cs2[i], DILI_N); // w - cs2
            POLY_ADD(w_cs2_ct0, diff, ct0[i], DILI_N); // w - cs2 + ct0
            
            poly neg_ct0;
            // Negate ct0: -ct0
            for(int j=0; j<DILI_N; j++) neg_ct0[j] = -ct0[i][j];

            hint_ones += make_hint(h[i], neg_ct0, w_cs2_ct0, DILI_N, DILI_ALPHA);
        }

        if(hint_ones > DILI_OMEGA) {
#ifdef DEBUG
            printf("hint check failed\n");
#endif            
            continue; // Too many hints, restart
        }

        // If we reached here, all checks passed!
        retry = 0;
    } while(retry && loop < 50);
    printf("Success at loop %d\n", loop);

    // pack Signature c_tilde (32) || z (L polys) || h (K polys + omega)
    
    uint8_t *sig_ptr = signature;
    
    // c_tilde
    memcpy(sig_ptr, c_tilde, 32);
    sig_ptr += 32;
    
    // z
    for(int i=0; i<L; i++) {
        pack_z(sig_ptr, z[i]);
        sig_ptr += (3 * DILI_N); // 768 bytes
    }
    
    // h
    if(pack_h(sig_ptr, h, K)) {
        // Should not happen if OMEGA check passed
        // Handle error or do nothing
#ifdef DEBUG
        printf("error!");
#endif
    }
}


// Verify
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk) {
    lp_ctx ctx;
    INIT_NTT_DILI(&ctx);
    
    const int K = DILI_K;
    const int L = DILI_L;
    
    const uint8_t *c_tilde = signature;
    const uint8_t *z_bytes = signature + 32;
    const uint8_t *h_bytes = z_bytes + (L * 3 * DILI_N);
    
    poly z[L];  
    uint8_t h[K][DILI_N]; 
    poly c; 
    
    for(int i = 0; i < L; i++) {
        unpack_z(z[i], z_bytes + i * 3 * DILI_N);
    }
    unpack_h((poly*)h, h_bytes, K);
    
    for(int i = 0; i < L; i++) {
        if(check_norm(z[i], DILI_GAMMA1 - DILI_BETA)) {
            return -1;
        }
    }
    
    uint8_t tr[64];
    uint8_t pk_packed[32 + K * 416];
    dili_pack_pk(pk_packed, pk->rho, pk->t1, K);
    shake256(tr, 64, pk_packed, sizeof(pk_packed));
    
    uint8_t mu[64];
    KECCAK1600_CTX shake_ctx;
    shake256_init(&shake_ctx);
    sha3_update(&shake_ctx, tr, 64);
    sha3_update(&shake_ctx, message, message_len);
    shake_xof_final(&shake_ctx);
    shake_squeeze(&shake_ctx, mu, 64);
    
    poly A[DILI_K * DILI_L];
    polyvec_matrix_expand(A, DILI_K, DILI_L, pk->rho, &ctx);

    challenge(c, c_tilde, DILI_TAU);
    poly c_ntt;
    POLY_COPY(c_ntt, c, DILI_N);
    poly_ntt(c_ntt, &ctx);
    
    poly z_ntt[DILI_L];
    polyvec_copy(z_ntt, z, DILI_L);
    polyvec_ntt(z_ntt, DILI_L, &ctx);
    
    poly Az_ntt_vec[DILI_K];
    polyvec_matrix_pointwise_montgomery(Az_ntt_vec, A, DILI_K, z_ntt, DILI_L, &ctx);
    
    poly Az_ntt[K];
    polyvec_copy(Az_ntt, Az_ntt_vec, DILI_K);
    
    poly t1_2d_ntt[K];
    for(int i = 0; i < K; i++) {
        poly t1_shifted;
        for(int j = 0; j < DILI_N; j++) {
            t1_shifted[j] = pk->t1[i][j] << DILI_D;
        }
        poly_ntt(t1_shifted, &ctx);
        POLY_COPY(t1_2d_ntt[i], t1_shifted, DILI_N);
    }
    
    poly w_prime_ntt[K];
    for(int i = 0; i < K; i++) {
        poly ct1_ntt;
        poly_pointwise(ct1_ntt, c_ntt, t1_2d_ntt[i], &ctx);
        POLY_SUB(w_prime_ntt[i], Az_ntt[i], ct1_ntt, DILI_N);
    }
    
    poly w_prime[K];
    for(int i = 0; i < K; i++) {
        POLY_COPY(w_prime[i], w_prime_ntt[i], DILI_N);
        poly_invntt(w_prime[i], &ctx);
        for(int j = 0; j < DILI_N; j++) {
            w_prime[i][j] = REDUCE32(w_prime[i][j]);
            if(w_prime[i][j] < 0) w_prime[i][j] += DILI_Q;
        }
    }
    
    poly w1_prime[K];
    for(int i = 0; i < K; i++) {
        use_hint(w1_prime[i], h[i], w_prime[i], DILI_N, DILI_ALPHA);
    }
    
    uint8_t c_tilde_prime[32];
    uint8_t w1_prime_packed[32 + K * DILI_N];
    memcpy(w1_prime_packed, mu, 32);
    int offset = 32;
    for(int i = 0; i < K; i++) {
        pack_w1(w1_prime_packed + offset, w1_prime[i]);
        offset += DILI_N;
    }
    shake256(c_tilde_prime, 32, w1_prime_packed, offset);
    
    for(int i = 0; i < 32; i++) {
        if(c_tilde[i] != c_tilde_prime[i]) {
            printf("challenge fail!");
            return -1;
        }
    }
    
    return 0;
}