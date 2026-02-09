#include "../inc/dili.h"
#include "../inc/poly.h"
#include "../inc/polyvec.h"
#include "../inc/pack.h"
#include "../inc/rand.h"
#include <string.h>

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

// Compare bytes: returns 0 if equal, non-zero if different
static int compare_bytes(const uint8_t *a, const uint8_t *b, int len) {
    int i;
    int diff = 0;
    for(i = 0; i < len; i++) {
        diff |= (a[i] ^ b[i]);
    }
    return diff;
}

// Power2Round: decompose r into (r1, r0) where r = r1*2^d + r0
static void power2round(poly *r1, poly *r0, poly *r, int k, int n, int d) {
    for(int j = 0; j < k; j++) {
        for(int i = 0; i < n; i++) {
            // Standard rounding: r1 = round(r / 2^d)
            r1[j][i] = (r[j][i] + (1 << (d - 1))) >> d;
            // r0 = r - r1 * 2^d
            r0[j][i] = r[j][i] - (r1[j][i] << d);
        }
    }
}

// Decompose: split r into high and low bits
static void decompose(lssi_t *r1, lssi_t *r0, const lssi_t *r, int n, int alpha) {
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

static void high_bits(lssi_t *r1, const lssi_t *r, int n, int alpha) {
    int i;
    lssi_t a, r1_val;
    
    for(i=0;i<n;++i){
        a = REDUCE32(r[i]);
        r1_val = (a + 127) >> 7;
        #if DILI_GAMMA2 == ((DILI_Q-1)/88)
        r1_val = (r1_val * 11275 + (1 << 23)) >> 24;
        r1_val ^= ((43 - r1_val) >> 31) & r1_val;
        #else
        r1_val = (r1_val * 1025 + (1 << 21)) >> 22;
        r1_val &= 15;
        #endif
        
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

static int make_hint(uint8_t *hint, const lssi_t *a, const lssi_t *b, int n, int alpha) {
    int weight = 0;
    int i;
    (void)alpha;
    
    for(i = 0; i < n; ++i) {
        lssi_t a_val = REDUCE32(a[i]);
        lssi_t b_val = REDUCE32(b[i]);
        lssi_t r1_b = (b_val + 127) >> 7;
        #if DILI_GAMMA2 == ((DILI_Q-1)/88)
        r1_b = (r1_b * 11275 + (1 << 23)) >> 24;
        r1_b ^= ((43 - r1_b) >> 31) & r1_b;
        #else
        r1_b = (r1_b * 1025 + (1 << 21)) >> 22;
        r1_b &= 15;
        #endif
        
        lssi_t sum = REDUCE32(a_val + b_val);
        lssi_t r1_sum = (sum + 127) >> 7;
        #if DILI_GAMMA2 == ((DILI_Q-1)/88)
        r1_sum = (r1_sum * 11275 + (1 << 23)) >> 24;
        r1_sum ^= ((43 - r1_sum) >> 31) & r1_sum;
        #else
        r1_sum = (r1_sum * 1025 + (1 << 21)) >> 22;
        r1_sum &= 15;
        #endif
        if(r1_b != r1_sum) {
            hint[i] = 1;
            weight++;
        } else {
            hint[i] = 0;
        }
    }
    
    return weight;
}

static void use_hint(lssi_t *r1, const uint8_t *hint, const lssi_t *r, int n, int alpha) {
    int i;
    
    for(i = 0; i < n; ++i) {
        lssi_t a = REDUCE32(r[i]);
        lssi_t r1_val = (a + 127) >> 7;
        #if DILI_GAMMA2 == ((DILI_Q-1)/88)
        r1_val = (r1_val * 11275 + (1 << 23)) >> 24;
        r1_val ^= ((43 - r1_val) >> 31) & r1_val;
        #else
        r1_val = (r1_val * 1025 + (1 << 21)) >> 22;
        r1_val &= 15;
        #endif
        
        if(hint[i] == 0) {
            r1[i] = r1_val;
            continue;
        }
        
        lssi_t r0 = a - r1_val * 2 * DILI_GAMMA2;
        r0 -= (((DILI_Q-1)/2 - r0) >> 31) & DILI_Q;
        #if DILI_GAMMA2 == ((DILI_Q-1)/88)
        if(r0 > 0) {
            r1_val = (r1_val == 43) ? 0 : (r1_val + 1);
        } else {
            r1_val = (r1_val == 0) ? 43 : (r1_val - 1);
        }
        #else
        if(r0 > 0) {
            r1_val = (r1_val + 1) & 15;
        } else {
            r1_val = (r1_val - 1) & 15;
        }
        #endif
        
        r1[i] = r1_val;
    }
}

static void expand_s(poly *s1, poly *s2, const uint8_t *rho_prime, int k, int l, int eta, lp_ctx *ctx) {
    uint8_t seed[66];
    uint8_t buf[256];
    int i;

    for(i = 0; i < l; i++) {
        memcpy(seed, rho_prime, 32);
        seed[32] = (uint8_t)i;
        seed[33] = 0x00;
        shake256(buf, sizeof(buf), seed, 34);
        sample_eta(s1[i], ctx->n, buf, sizeof(buf), ctx, eta);
    }
    for(i = 0; i < k; i++) {
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
    uint8_t buf[1024];
    KECCAK1600_CTX ctx;
    int ctr, i;
    uint32_t val, bound;
    int32_t diff, mask;
    
    shake256_init(&ctx);
    sha3_update(&ctx, seed, seed_len);
    shake_xof_final(&ctx);
    
    ctr = 0;
    bound = 2 * gamma1 - 1;
    
    while(ctr < len){
        shake_squeeze(&ctx, buf, sizeof(buf));
        
        for(i=0; i<sizeof(buf)-3 && ctr<len; i+=3){
            val = buf[i];
            val |= (uint32_t)buf[i+1] << 8;
            val |= (uint32_t)buf[i+2] << 16;
            val &= 0xFFFFF;
            
            diff = (int32_t)val - (int32_t)bound;
            mask = SGN_MASK(diff);
            
            if(mask){
                y[ctr++] = (lssi_t)(gamma1 - 1) - (lssi_t)val;
            }
        }
    }    
}

// Challenge: hash to challenge polynomial c
static void challenge(lssi_t *c, const uint8_t *hash, int tau) {
    uint8_t buf[256];
    uint64_t signs;
    int i, b, pos;
    
    signs = 0;
    shake256(buf, 256, hash, 32);

    for(i = 0; i < 8; ++i)
        signs |= (uint64_t)buf[i] << (8 * i);

    for(i = 0; i < DILI_N; ++i) 
        c[i] = 0;

    pos = 8;
    for(i = DILI_N - tau; i < DILI_N; ++i) {
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
    uint8_t seed[32];
    uint8_t buf[96];
    uint8_t rho[32], rho_prime[32], K[32];
    uint8_t pk_packed[32 + DILI_K * 416];
    poly matrix_A[DILI_K * DILI_L];
    poly s1_ntt[DILI_L];
    poly t[DILI_K];
    poly s2[DILI_K];
    poly t_arr[DILI_K];
    
    INIT_NTT_DILI(&ctx);
    rand_bytes(seed, 32);
    
    shake256(buf, 96, seed, 32);
    memcpy(rho, buf, 32);
    memcpy(rho_prime, buf + 32, 32);
    memcpy(K, buf + 64, 32);

    polyvec_matrix_expand(matrix_A, param.k, param.l, rho, &ctx);
    
    expand_s(keypair->sk.s1, keypair->sk.s2, rho_prime, param.k, param.l, param.eta, &ctx);
    
    polyvec_copy(s1_ntt, keypair->sk.s1, param.l);
    polyvec_ntt(s1_ntt, param.l, &ctx);

    polyvec_copy(s2, keypair->sk.s2, param.k);
    calc_t(t, matrix_A, param.k, param.l, &ctx, s1_ntt, s2);
    
    polyvec_copy(t_arr, t, param.k);
    
    power2round(keypair->pk.t1, keypair->sk.t0, t_arr, param.k, ctx.n, DILI_D);

    dili_pack_pk(pk_packed, rho, keypair->pk.t1, param.k);
    
    shake256(keypair->sk.tr, 64, pk_packed, 32 + param.k * 416);

    memcpy(keypair->pk.rho, rho, 32);
    memcpy(keypair->sk.rho, rho, 32);
    memcpy(keypair->sk.K, K, 32);
}

// Sign
int dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk) {
    lp_ctx ctx;
    const int K = DILI_K;
    const int L = DILI_L;
    
    /* All variable declarations at the beginning */
    poly matrix_A[K * L];
    poly y[L];
    poly w[K], w1[K];
    poly cp, cp_ntt;
    poly z[L];
    poly cs1[L], cs2[K], ct0[K];
    poly s1_ntt[L];
    poly s2_ntt[K], t0_ntt[K];
    poly y_ntt[L];
    poly w_cs2_ct0[K];
    poly diff, r0;
    poly neg_ct0;
    uint8_t h[K][DILI_N];
    uint8_t mu[64];
    uint8_t c_tilde[32];
    uint8_t c_tilde_temp_buf[64 + K * DILI_N];
    uint8_t rho_prime_nonce[64];
    uint8_t y_seed_base[64];
    uint8_t seed_with_nonce[66];
    uint8_t *sig_ptr;
    KECCAK1600_CTX shake_ctx;
    uint16_t nonce;
    int retry, loop, offset;
    int z_check_fail, r0_check_fail, ct0_check_fail;
    int hint_ones;
    int i, j;
    
    INIT_NTT_DILI(&ctx);
    
    nonce = 0;
    retry = 1;
    loop = 0;

    polyvec_matrix_expand(matrix_A, K, L, sk->rho, &ctx);

    shake256_init(&shake_ctx);
    sha3_update(&shake_ctx, sk->tr, 64);
    sha3_update(&shake_ctx, message, message_len);
    shake_xof_final(&shake_ctx);
    shake_squeeze(&shake_ctx, mu, 64);

    /* Rejection Sampling Loop */
    do{
        ++loop;

        memcpy(rho_prime_nonce, sk->K, 32);
        memcpy(rho_prime_nonce + 32, mu, 32);
        
        shake256(y_seed_base, 64, rho_prime_nonce, 64);

        for(i=0; i<L; i++) {
            memcpy(seed_with_nonce, y_seed_base, 64);
            seed_with_nonce[64] = nonce & 0xFF;
            seed_with_nonce[65] = (nonce >> 8) & 0xFF;
            nonce++;
            
            expand_mask(y[i], seed_with_nonce, 66, DILI_N, DILI_GAMMA1);
        }

        polyvec_copy(y_ntt, y, L);
        polyvec_ntt(y_ntt, L, &ctx);

        polyvec_matrix_pointwise_montgomery(w, matrix_A, K, y_ntt, L, &ctx);
        polyvec_invntt_tomont(w, K, &ctx);

        for(i=0; i<K; i++) {
            high_bits(w1[i], w[i], DILI_N, DILI_ALPHA);
        }

        memcpy(c_tilde_temp_buf, mu, 64);
        offset = 64;
        for(i=0; i<K; i++) {
            pack_w1(c_tilde_temp_buf + offset, w1[i]);
            offset += (DILI_N / 4) * 3;
        }
        shake256(c_tilde, 32, c_tilde_temp_buf, offset);

        challenge(cp, c_tilde, DILI_TAU);
        POLY_COPY(cp_ntt, cp, DILI_N);
        poly_ntt(cp_ntt, &ctx);

        polyvec_copy(s1_ntt, sk->s1, L);
        polyvec_ntt(s1_ntt, L, &ctx);

        for(i=0; i<L; i++) {
            poly_pointwise(cs1[i], cp_ntt, s1_ntt[i], &ctx);
        }

        polyvec_invntt_tomont(cs1, L, &ctx);
        polyvec_add(z, y, cs1, L);
        
        z_check_fail = 0;
        for(i=0; i<L; i++) {
            if(check_norm(z[i], DILI_GAMMA1 - DILI_BETA)) {
                z_check_fail = 1;
                break;
            }
        }
        if(z_check_fail) continue;

        for(i=0; i<K; ++i){
            POLY_COPY(s2_ntt[i], sk->s2[i], DILI_N);
            POLY_COPY(t0_ntt[i], sk->t0[i], DILI_N);
        }
        polyvec_ntt(s2_ntt, K, &ctx);
        polyvec_ntt(t0_ntt, K, &ctx);
        for(i=0; i<K; i++) {
            poly_pointwise(cs2[i], cp_ntt, s2_ntt[i], &ctx);
            poly_pointwise(ct0[i], cp_ntt, t0_ntt[i], &ctx);
        }
        polyvec_invntt_tomont(cs2, K, &ctx);
        polyvec_invntt_tomont(ct0, K, &ctx);
        
        r0_check_fail = 0;
        for(i=0; i<K; i++) {
            POLY_SUB(diff, w[i], cs2[i], DILI_N);
            low_bits(r0, diff, DILI_N, DILI_ALPHA);
            
            if(check_norm(r0, DILI_GAMMA2 - DILI_BETA)) {
                r0_check_fail = 1;
                break;
            }
        }
        if(r0_check_fail) continue;

        ct0_check_fail = 0;
        for(i=0; i<K; i++) {
            if(check_norm(ct0[i], DILI_GAMMA2)) {
                ct0_check_fail = 1;
                break;
            }
        }
        if(ct0_check_fail) continue;

        for(i=0; i<K; i++) {
            POLY_SUB(diff, w[i], cs2[i], DILI_N);
            POLY_ADD(w_cs2_ct0[i], diff, ct0[i], DILI_N);
        }
        
        polyvec_caddq(w_cs2_ct0, K);
        
        hint_ones = 0;
        for(i=0; i<K; i++) {
            for(j=0; j<DILI_N; j++) 
                neg_ct0[j] = -ct0[i][j];

            hint_ones += make_hint(h[i], neg_ct0, w_cs2_ct0[i], DILI_N, DILI_ALPHA);
        }

        if(hint_ones > DILI_OMEGA) {
            continue;
        }

        retry = 0;
    } while(retry && loop < 300);
    
    if(loop >= 150) {
        return -1;
    }

    sig_ptr = signature;
    
    memcpy(sig_ptr, c_tilde, 32);
    sig_ptr += 32;
    
    for(i=0; i<L; i++) {
        pack_z(sig_ptr, z[i]);
        sig_ptr += (3 * DILI_N);
    }
    
    if(pack_h(sig_ptr, h, K)) {
        return -1;
    }
    return 0;
}


// Verify
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk) {
    lp_ctx ctx;
    const int K = DILI_K;
    const int L = DILI_L;
    
    /* All variable declarations at the beginning */
    const uint8_t *c_tilde;
    const uint8_t *z_bytes;
    const uint8_t *h_bytes;
    poly z[L];
    uint8_t h[K][DILI_N];
    poly c, c_ntt;
    poly z_ntt[DILI_L];
    poly Az_ntt_vec[DILI_K];
    poly Az_ntt[K];
    poly t1_2d_ntt[K];
    poly t1_shifted;
    poly w_prime_ntt[K];
    poly ct1_ntt;
    poly w_prime[K];
    poly w1_prime[K];
    poly A[DILI_K * DILI_L];
    uint8_t tr[64];
    uint8_t pk_packed[32 + K * 416];
    uint8_t mu[64];
    uint8_t c_tilde_prime[32];
    uint8_t w1_prime_packed[64 + K * ((DILI_N / 4) * 3)];
    KECCAK1600_CTX shake_ctx;
    int rc_unpack;
    int offset;
    int i, j;
    int32_t t1_val;
    int64_t shifted;
    
    INIT_NTT_DILI(&ctx);
    
    c_tilde = signature;
    z_bytes = signature + 32;
    h_bytes = z_bytes + (L * 3 * DILI_N);
    
    for(i = 0; i < L; i++) {
        unpack_z(z[i], z_bytes + i * 3 * DILI_N);
    }
    
    for(i = 0; i < L; i++) {
        if(check_norm(z[i], DILI_GAMMA1 - DILI_BETA)) {
            return -1;
        }
    }
    
    dili_pack_pk(pk_packed, pk->rho, pk->t1, K);
    shake256(tr, 64, pk_packed, sizeof(pk_packed));
    
    shake256_init(&shake_ctx);
    sha3_update(&shake_ctx, tr, 64);
    sha3_update(&shake_ctx, message, message_len);
    shake_xof_final(&shake_ctx);
    shake_squeeze(&shake_ctx, mu, 64);

    polyvec_matrix_expand(A, DILI_K, DILI_L, pk->rho, &ctx);

    challenge(c, c_tilde, DILI_TAU);
    POLY_COPY(c_ntt, c, DILI_N);
    poly_ntt(c_ntt, &ctx);
    
    polyvec_copy(z_ntt, z, DILI_L);
    polyvec_ntt(z_ntt, DILI_L, &ctx);
    
    polyvec_matrix_pointwise_montgomery(Az_ntt_vec, A, DILI_K, z_ntt, DILI_L, &ctx);
    polyvec_copy(Az_ntt, Az_ntt_vec, DILI_K);
    
    for(i = 0; i < K; i++) {
        for(j = 0; j < DILI_N; j++) {
            t1_val = pk->t1[i][j];
            t1_val = t1_val % DILI_Q;
            if(t1_val < 0) t1_val += DILI_Q;
            
            shifted = ((int64_t)t1_val << DILI_D);
            t1_shifted[j] = (lssi_t)(shifted % DILI_Q);
        }
        poly_ntt(t1_shifted, &ctx);
        POLY_COPY(t1_2d_ntt[i], t1_shifted, DILI_N);
    }
    
    for(i = 0; i < K; i++) {
        poly_pointwise(ct1_ntt, c_ntt, t1_2d_ntt[i], &ctx);
        POLY_SUB(w_prime_ntt[i], Az_ntt[i], ct1_ntt, DILI_N);
    }
    
    polyvec_copy(w_prime, w_prime_ntt, K);
    polyvec_invntt_tomont(w_prime, K, &ctx);
    polyvec_caddq(w_prime, K);
    
    rc_unpack = unpack_h(h, h_bytes, K);
    if(rc_unpack != 0) {
        return -1;
    }

    for(i = 0; i < K; i++) {
        use_hint(w1_prime[i], h[i], w_prime[i], DILI_N, DILI_ALPHA);
    }

    memcpy(w1_prime_packed, mu, 64);
    offset = 64;
    for(i = 0; i < K; i++) {
        pack_w1(w1_prime_packed + offset, w1_prime[i]);
        offset += (DILI_N / 4) * 3;
    }

    shake256(c_tilde_prime, 32, w1_prime_packed, offset);
    
    if(compare_bytes(c_tilde, c_tilde_prime, 32) != 0) {
        return -1;
    }
    
    return 0;
}