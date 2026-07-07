#include "dili.h"
#include "poly.h"
#include "polyvec.h"
#include "pack.h"
#include "rand.h"
#include "dili_sample.h"
#include <string.h>

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
    ((a)&_c) | ((b)&(~_c)); \
}) // return cond?a:b;

#define REDUCE32(x) ({\
    int32_t _a = (x); \
    int32_t _t = (_a + (1<<22)) >> 23; \
    _t = _a - _t*DILI_Q; \
    _t += (DILI_Q & SGN_MASK(_t)); \
    _t; \
}) // int32 -> [0,q-1], 无除法

// Compare bytes: returns 0 if equal, non-zero if different
static int compare_bytes(const uint8_t *a, const uint8_t *b, int len) {
    int i;
    int diff = 0;
    for(i = 0; i < len; i++) {
        diff |= (a[i] ^ b[i]);
    }
    return diff;
}

static void power2round(poly *r1, poly *r0, poly *r, int k, int n, int d) {
    for(int j = 0; j < k; j++) {
        for(int i = 0; i < n; i++) {
            // r1 = round(r / 2^d)，与参考实现一致 (a + 2^(d-1) - 1) >> d
            r1[j][i] = (r[j][i] + (1 << (d - 1)) - 1) >> d;
            // r0 = r - r1 * 2^d
            r0[j][i] = r[j][i] - (r1[j][i] << d);
        }
    }
}

static void decompose(lssi_t *r1, lssi_t *r0, const lssi_t *r, int n, int alpha) {
    int i;
    int gamma2 = alpha / 2;
    lssi_t a, a0, a1;
    for(i=0;i<n;++i){
        a = REDUCE32(r[i]);
        a1 = (a + 127) >> 7;
        if(gamma2 == ((DILI_Q-1)/88)) {  // Dilithium2
            a1 = (a1 * 11275 + (1 << 23)) >> 24;
            a1 ^= ((43 - a1) >> 31) & a1;
        } else {                          // Dilithium3/5
            a1 = (a1 * 1025 + (1 << 21)) >> 22;
            a1 &= 15;
        }
        a0 = a - a1 * 2 * gamma2;
        a0 -= (((DILI_Q-1)/2 - a0) >> 31) & DILI_Q;
        r0[i] = a0;
        r1[i] = a1;
    }
}

static void high_bits(lssi_t *r1, const lssi_t *r, int n, int alpha) {
    int i;
    lssi_t a, r1_val;
    int gamma2 = alpha / 2;

    if(gamma2 == ((DILI_Q-1)/88)) {  // Dilithium2: gamma2 = 95232
        for(i=0;i<n;++i){
            a = REDUCE32(r[i]);
            r1_val = (a + 127) >> 7;
            r1_val = (r1_val * 11275 + (1 << 23)) >> 24;
            r1_val ^= ((43 - r1_val) >> 31) & r1_val;
            r1[i] = r1_val;
        }
    } else {  // Dilithium3/5: gamma2 = 261888
        for(i=0;i<n;++i){
            a = REDUCE32(r[i]);
            r1_val = (a + 127) >> 7;
            r1_val = (r1_val * 1025 + (1 << 21)) >> 22;
            r1_val &= 15;
            r1[i] = r1_val;
        }
    }
}

static int check_norm(const poly p, int bound){
    int fail, i;
    lssi_t val, abs_val, diff;
    fail = 0;
    for(i=0;i<DILI_N;++i){
        val = REDUCE32(p[i]);
        val -= (((DILI_Q >> 1) - val) >> 31) & DILI_Q;
        abs_val = ABS_VAL(val);
        diff = (bound-1)-abs_val;
        fail |= SGN_MASK(diff);
    }
    return (fail&1);
}

static void low_bits(lssi_t *r0, const lssi_t *r, int n, int alpha) {
    int i;
    int gamma2 = alpha / 2;
    lssi_t a, a0, a1;

    if(gamma2 == ((DILI_Q-1)/88)) {  // Dilithium2
        for(i=0;i<n;++i){
            a = REDUCE32(r[i]);
            a1 = (a + 127) >> 7;
            a1 = (a1 * 11275 + (1 << 23)) >> 24;
            a1 ^= ((43 - a1) >> 31) & a1;
            a0 = a - a1 * 2 * gamma2;
            a0 -= (((DILI_Q-1)/2 - a0) >> 31) & DILI_Q;
            r0[i] = a0;
        }
    } else {                          // Dilithium3/5
        for(i=0;i<n;++i){
            a = REDUCE32(r[i]);
            a1 = (a + 127) >> 7;
            a1 = (a1 * 1025 + (1 << 21)) >> 22;
            a1 &= 15;
            a0 = a - a1 * 2 * gamma2;
            a0 -= (((DILI_Q-1)/2 - a0) >> 31) & DILI_Q;
            r0[i] = a0;
        }
    }
}

static int make_hint_from_ct0(uint8_t *hint, const lssi_t *ct0, const lssi_t *b, int n, int alpha) {
    int weight = 0;
    int i;
    int gamma2 = alpha / 2;

    if(gamma2 == ((DILI_Q-1)/88)) {  // Dilithium2
        for(i = 0; i < n; ++i) {
            lssi_t a_val = REDUCE32(-ct0[i]);
            lssi_t b_val = REDUCE32(b[i]);
            lssi_t r1_b = (b_val + 127) >> 7;
            r1_b = (r1_b * 11275 + (1 << 23)) >> 24;
            r1_b ^= ((43 - r1_b) >> 31) & r1_b;

            lssi_t sum = REDUCE32(a_val + b_val);
            lssi_t r1_sum = (sum + 127) >> 7;
            r1_sum = (r1_sum * 11275 + (1 << 23)) >> 24;
            r1_sum ^= ((43 - r1_sum) >> 31) & r1_sum;

            if(r1_b != r1_sum) {
                hint[i] = 1;
                weight++;
            } else {
                hint[i] = 0;
            }
        }
    } else {  // Dilithium3/5
        for(i = 0; i < n; ++i) {
            lssi_t a_val = REDUCE32(-ct0[i]);
            lssi_t b_val = REDUCE32(b[i]);
            lssi_t r1_b = (b_val + 127) >> 7;
            r1_b = (r1_b * 1025 + (1 << 21)) >> 22;
            r1_b &= 15;

            lssi_t sum = REDUCE32(a_val + b_val);
            lssi_t r1_sum = (sum + 127) >> 7;
            r1_sum = (r1_sum * 1025 + (1 << 21)) >> 22;
            r1_sum &= 15;

            if(r1_b != r1_sum) {
                hint[i] = 1;
                weight++;
            } else {
                hint[i] = 0;
            }
        }
    }
    
    return weight;
}

static void use_hint(lssi_t *r1, const uint8_t *hint, const lssi_t *r, int n, int alpha) {
    int i;
    int gamma2 = alpha / 2;

    if(gamma2 == ((DILI_Q-1)/88)) {  // Dilithium2
        for(i = 0; i < n; ++i) {
            lssi_t a = REDUCE32(r[i]);
            lssi_t r1_val = (a + 127) >> 7;
            r1_val = (r1_val * 11275 + (1 << 23)) >> 24;
            r1_val ^= ((43 - r1_val) >> 31) & r1_val;

            if(hint[i] == 0) {
                r1[i] = r1_val;
                continue;
            }

            lssi_t r0 = a - r1_val * 2 * gamma2;
            r0 -= (((DILI_Q-1)/2 - r0) >> 31) & DILI_Q;

            if(r0 > 0) {
                r1_val = (r1_val == 43) ? 0 : (r1_val + 1);
            } else {
                r1_val = (r1_val == 0) ? 43 : (r1_val - 1);
            }

            r1[i] = r1_val;
        }
    } else {  // Dilithium3/5
        for(i = 0; i < n; ++i) {
            lssi_t a = REDUCE32(r[i]);
            lssi_t r1_val = (a + 127) >> 7;
            r1_val = (r1_val * 1025 + (1 << 21)) >> 22;
            r1_val &= 15;

            if(hint[i] != 0) {
                lssi_t r0 = a - r1_val * 2 * gamma2;
                r0 -= (((DILI_Q-1)/2 - r0) >> 31) & DILI_Q;

                if(r0 > 0) {
                    r1_val = (r1_val + 1) & 15;
                } else {
                    r1_val = (r1_val - 1) & 15;
                }
            }

            r1[i] = r1_val;
        }
    }
}

static void expand_s(poly *s1, poly *s2, const uint8_t *rho_prime, int k, int l, int eta, lp_ctx *ctx) {
    uint8_t seed[66];
    uint8_t buf[256];
    KECCAK1600_CTX kc;
    int i, got;

    for(i = 0; i < l + k; i++) {
        poly *dst = (i < l) ? &s1[i] : &s2[i - l];
        memcpy(seed, rho_prime, 64);          // rho' = 64 字节
        seed[64] = (uint8_t)(i & 0xFF);        // nonce (2 字节小端)
        seed[65] = (uint8_t)((i >> 8) & 0xFF);

        shake256_init(&kc);
        sha3_update(&kc, seed, 66);
        shake_xof_final(&kc);

        // 拒绝采样：持续 squeeze 直到填满全部系数
        got = 0;
        while(got < ctx->n) {
            shake_squeeze(&kc, buf, sizeof(buf));
            got += sample_eta((*dst) + got, ctx->n - got, buf, sizeof(buf), ctx, eta);
        }
    }
}

static void calc_t(poly *t, poly *A_mat, int k, int l, lp_ctx *ctx, poly *s1_ntt, poly *s2) {
    // t = A*s1 + s2 (计算在NTT域)
    polyvec_matrix_pointwise_montgomery(t, A_mat, k, s1_ntt, l, ctx);
    polyvec_reduce(t, k);            // 累加了 l 个乘积, invNTT 前必须归约防止 int32 溢出
    polyvec_invntt_tomont(t, k, ctx);
    polyvec_add(t, t, s2, k);
}

static void expand_mask(lssi_t *y, const uint8_t *seed, int seed_len, int len, int gamma1) {
    uint8_t buf[640];   // 最大 POLYZ_PACKEDBYTES (gamma1=2^19)
    KECCAK1600_CTX ctx;
    int nbytes = DILI_POLYZ_PACKEDBYTES(gamma1);
    (void)len;

    shake256_init(&ctx);
    sha3_update(&ctx, seed, seed_len);
    shake_xof_final(&ctx);
    shake_squeeze(&ctx, buf, nbytes);

    unpack_z(y, buf, gamma1);
}

// Challenge: hash to challenge polynomial c (SampleInBall), c_tilde 长度 = lambda/4
static void challenge(lssi_t *c, const uint8_t *hash, int tau, int c_tilde_bytes) {
    uint8_t buf[136];   // SHAKE256 rate
    KECCAK1600_CTX kc;
    uint64_t signs;
    int i, b, pos;

    shake256_init(&kc);
    sha3_update(&kc, hash, c_tilde_bytes);
    shake_xof_final(&kc);
    shake_squeeze(&kc, buf, sizeof(buf));

    signs = 0;
    for(i = 0; i < 8; ++i)
        signs |= (uint64_t)buf[i] << (8 * i);

    for(i = 0; i < DILI_N; ++i)
        c[i] = 0;

    pos = 8;
    for(i = DILI_N - tau; i < DILI_N; ++i) {
        do {
            if(pos >= (int)sizeof(buf)) {
                shake_squeeze(&kc, buf, sizeof(buf));  // 继续 XOF 流, 取新字节
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
    const int K_DILI = param.k;
    const int L_DILI = param.l;
    lp_ctx ctx;
    uint8_t seed[32];
    uint8_t seedbuf[34];
    uint8_t buf[128];
    uint8_t rho[32], rho_prime[64], K[32];
    uint8_t pk_packed[32 + K_DILI * DILI_POLYT1_PACKEDBYTES];
    poly matrix_A[K_DILI * L_DILI];
    poly s1_ntt[L_DILI];
    poly t[K_DILI];

    INIT_NTT_DILI(&ctx);
    rand_bytes(seed, 32);

    memcpy(seedbuf, seed, 32);
    seedbuf[32] = (uint8_t)K_DILI;         // 域分隔: IntegerToBytes(k,1)
    seedbuf[33] = (uint8_t)L_DILI;         //         IntegerToBytes(l,1)
    shake256(buf, 128, seedbuf, 34);
    memcpy(rho, buf, 32);
    memcpy(rho_prime, buf + 32, 64);   // rho' = 64 字节 (标准)
    memcpy(K, buf + 96, 32);

    polyvec_matrix_expand(matrix_A, K_DILI, L_DILI, rho, &ctx);
    
    expand_s(keypair->sk.s1, keypair->sk.s2, rho_prime, K_DILI, L_DILI, param.eta, &ctx);
    
    polyvec_copy(s1_ntt, keypair->sk.s1, L_DILI);
    polyvec_ntt(s1_ntt, L_DILI, &ctx);

    calc_t(t, matrix_A, K_DILI, L_DILI, &ctx, s1_ntt, keypair->sk.s2);
    
    polyvec_caddq(t, K_DILI);   // 归约到 [0,q), 使 t1 落在规范 10-bit 范围便于打包

    power2round(keypair->pk.t1, keypair->sk.t0, t, K_DILI, ctx.n, DILI_D);

    dili_pack_pk(pk_packed, rho, keypair->pk.t1, K_DILI);
    
    shake256(keypair->sk.tr, 64, pk_packed, sizeof(pk_packed));

    memcpy(keypair->pk.rho, rho, 32);
    memcpy(keypair->sk.rho, rho, 32);
    memcpy(keypair->sk.K, K, 32);
}

// Sign
int dili_sign(uint8_t* signature, const uint8_t* message, int message_len, const dili_sk* sk, dili_param param) {
    lp_ctx ctx;
    const int K = param.k;
    const int L = param.l;
    
    poly matrix_A[K * L];
    poly y[L];
    poly w[K], w1[K];
    poly cp_ntt;
    poly z[L];
    poly cs1[L], cs2[K], ct0[K];
    poly s1_ntt[L];
    poly s2_ntt[K], t0_ntt[K];
    poly y_ntt[L];
    poly w_cs2_ct0[K];
    poly r0;
    uint8_t h[K][DILI_N];
    uint8_t mu[64];
    uint8_t rnd[32];
    uint8_t c_tilde[64];
    uint8_t c_tilde_temp_buf[64 + K * DILI_N];
    uint8_t rho_prime_nonce[128];  // K(32) || rnd(32) || mu(64)
    uint8_t rhoprimepp[64];        // rho'' = H(K || rnd || mu)
    uint8_t seed_with_nonce[66];
    uint8_t *sig_ptr;
    KECCAK1600_CTX shake_ctx;
    uint16_t nonce;
    int retry, loop, offset;
    int z_check_fail, r0_check_fail, ct0_check_fail;
    int hint_ones;
    int i;
    
    INIT_NTT_DILI(&ctx);
    
    nonce = 0;
    retry = 1;
    loop = 0;

    polyvec_matrix_expand(matrix_A, K, L, sk->rho, &ctx);

    // mu = H(tr(64) || M'),  M' = 0x00 || |ctx| || ctx || M  (纯 ML-DSA, 空 ctx)
    {
        uint8_t domsep[2] = {0x00, 0x00};   // 域分隔字节 + ctx 长度(0)
        shake256_init(&shake_ctx);
        sha3_update(&shake_ctx, sk->tr, 64);
        sha3_update(&shake_ctx, domsep, 2);
        sha3_update(&shake_ctx, message, message_len);
        shake_xof_final(&shake_ctx);
        shake_squeeze(&shake_ctx, mu, 64);
    }

    rand_bytes(rnd, 32);   // hedged 签名: 新鲜 32 字节 (确定性变体为全 0)

    // rho'' = H(K(32) || rnd(32) || mu(64), 64), 循环外常量
    memcpy(rho_prime_nonce, sk->K, 32);
    memcpy(rho_prime_nonce + 32, rnd, 32);
    memcpy(rho_prime_nonce + 64, mu, 64);
    shake256(rhoprimepp, 64, rho_prime_nonce, 128);

    polyvec_copy(s1_ntt, sk->s1, L);
    polyvec_ntt(s1_ntt, L, &ctx);
    polyvec_copy(s2_ntt, sk->s2, K);
    polyvec_copy(t0_ntt, sk->t0, K);
    polyvec_ntt(s2_ntt, K, &ctx);
    polyvec_ntt(t0_ntt, K, &ctx);
    memcpy(seed_with_nonce, rhoprimepp, 64);
    memcpy(c_tilde_temp_buf, mu, 64);

    do{
        ++loop;

        for(i=0; i<L; i++) {
            seed_with_nonce[64] = nonce & 0xFF;
            seed_with_nonce[65] = (nonce >> 8) & 0xFF;
            nonce++;

            expand_mask(y[i], seed_with_nonce, 66, DILI_N, param.gamma1);
        }

        polyvec_copy(y_ntt, y, L);
        polyvec_ntt(y_ntt, L, &ctx);

        polyvec_matrix_pointwise_montgomery(w, matrix_A, K, y_ntt, L, &ctx);
        polyvec_reduce(w, K);
        polyvec_invntt_tomont(w, K, &ctx);

        for(i=0; i<K; i++) {
            high_bits(w1[i], w[i], DILI_N, 2 * param.gamma2);
        }

        offset = 64;
        for(i=0; i<K; i++) {
            pack_w1(c_tilde_temp_buf + offset, w1[i], param.gamma2);
            offset += param.w1_packed_bytes;
        }
        shake256(c_tilde, param.c_tilde_bytes, c_tilde_temp_buf, offset);

        challenge(cp_ntt, c_tilde, param.tau, param.c_tilde_bytes);
        poly_ntt(cp_ntt, &ctx);

        for(i=0; i<L; i++) {
            poly_pointwise(cs1[i], cp_ntt, s1_ntt[i], &ctx);
        }

        polyvec_invntt_tomont(cs1, L, &ctx);
        polyvec_add(z, y, cs1, L);
        
        z_check_fail = 0;
        for(i=0; i<L; i++) {
            if(check_norm(z[i], param.gamma1 - param.beta)) {
                z_check_fail = 1;
                break;
            }
        }
        if(z_check_fail) continue;

        for(i=0; i<K; i++) {
            poly_pointwise(cs2[i], cp_ntt, s2_ntt[i], &ctx);
            poly_pointwise(ct0[i], cp_ntt, t0_ntt[i], &ctx);
        }
        polyvec_invntt_tomont(cs2, K, &ctx);
        polyvec_invntt_tomont(ct0, K, &ctx);
        
        r0_check_fail = 0;
        for(i=0; i<K; i++) {
            POLY_SUB(w_cs2_ct0[i], w[i], cs2[i], DILI_N);
            low_bits(r0, w_cs2_ct0[i], DILI_N, 2 * param.gamma2);
            
            if(check_norm(r0, param.gamma2 - param.beta)) {
                r0_check_fail = 1;
                break;
            }
        }
        if(r0_check_fail) continue;

        ct0_check_fail = 0;
        for(i=0; i<K; i++) {
            if(check_norm(ct0[i], param.gamma2)) {
                ct0_check_fail = 1;
                break;
            }
        }
        if(ct0_check_fail) continue;

        for(i=0; i<K; i++) {
            POLY_ADD(w_cs2_ct0[i], w_cs2_ct0[i], ct0[i], DILI_N);
        }
        
        polyvec_caddq(w_cs2_ct0, K);
        
        hint_ones = 0;
        for(i=0; i<K; i++) {
            hint_ones += make_hint_from_ct0(h[i], ct0[i], w_cs2_ct0[i], DILI_N, 2 * param.gamma2);
        }

        if(hint_ones > param.omega) {
            continue;
        }

        retry = 0;
    } while(retry && loop < 814);

    if(retry) {
        return -1;
    }


    sig_ptr = signature;

    memcpy(sig_ptr, c_tilde, param.c_tilde_bytes);
    sig_ptr += param.c_tilde_bytes;
    
    for(i=0; i<L; i++) {
        pack_z(sig_ptr, z[i], param.gamma1);
        sig_ptr += DILI_POLYZ_PACKEDBYTES(param.gamma1);
    }
    
    if(pack_h(sig_ptr, h, K, param.omega)) {
        return -1;
    }
    return 0;
}


// Verify
int dili_verify(const uint8_t* signature, const uint8_t* message, int message_len, const dili_pk* pk, dili_param param) {
    lp_ctx ctx;
    const int K = param.k;
    const int L = param.l;
    
    const uint8_t *c_tilde;
    const uint8_t *z_bytes;
    const uint8_t *h_bytes;
    poly z[L];
    uint8_t h[K][DILI_N];
    poly c_ntt;
    poly z_ntt[L];
    poly Az_ntt[K];
    poly t1_shifted;
    poly ct1_ntt;
    poly w_prime[K];
    poly w1_prime[K];
    poly A[K * L];
    uint8_t tr[64];
    uint8_t pk_packed[32 + K * DILI_POLYT1_PACKEDBYTES];
    uint8_t mu[64];
    uint8_t c_tilde_prime[64];
    uint8_t w1_prime_packed[64 + K * ((DILI_N / 4) * 3)];
    KECCAK1600_CTX shake_ctx;
    int rc_unpack;
    int offset;
    int i, j;

    INIT_NTT_DILI(&ctx);

    c_tilde = signature;
    z_bytes = signature + param.c_tilde_bytes;
    h_bytes = z_bytes + (L * DILI_POLYZ_PACKEDBYTES(param.gamma1));

    for(i = 0; i < L; i++) {
        unpack_z(z[i], z_bytes + i * DILI_POLYZ_PACKEDBYTES(param.gamma1), param.gamma1);
    }

    for(i = 0; i < L; i++) {
        if(check_norm(z[i], param.gamma1 - param.beta)) {
            return -1;
        }
    }

    dili_pack_pk(pk_packed, pk->rho, pk->t1, K);
    shake256(tr, 64, pk_packed, sizeof(pk_packed));

    // mu = H(tr(64) || 0x00 || |ctx|(0) || M),  纯 ML-DSA 空 ctx
    {
        uint8_t domsep[2] = {0x00, 0x00};
        shake256_init(&shake_ctx);
        sha3_update(&shake_ctx, tr, 64);
        sha3_update(&shake_ctx, domsep, 2);
        sha3_update(&shake_ctx, message, message_len);
        shake_xof_final(&shake_ctx);
        shake_squeeze(&shake_ctx, mu, 64);
    }

    polyvec_matrix_expand(A, K, L, pk->rho, &ctx);

    challenge(c_ntt, c_tilde, param.tau, param.c_tilde_bytes);
    poly_ntt(c_ntt, &ctx);

    polyvec_copy(z_ntt, z, L);
    polyvec_ntt(z_ntt, L, &ctx);
    
    polyvec_matrix_pointwise_montgomery(Az_ntt, A, K, z_ntt, L, &ctx);
    
    for(i = 0; i < K; i++) {
        for(j = 0; j < DILI_N; j++) {
            t1_shifted[j] = pk->t1[i][j] << DILI_D;
        }
        poly_ntt(t1_shifted, &ctx);
        poly_pointwise(ct1_ntt, c_ntt, t1_shifted, &ctx);
        POLY_SUB(w_prime[i], Az_ntt[i], ct1_ntt, DILI_N);
    }
    
    polyvec_reduce(w_prime, K);     // invNTT 前归约, 防止累加溢出
    polyvec_invntt_tomont(w_prime, K, &ctx);
    polyvec_caddq(w_prime, K);
    
    rc_unpack = unpack_h(h, h_bytes, K, param.omega);
    if(rc_unpack != 0) {
        return -1;
    }

    for(i = 0; i < K; i++) {
        use_hint(w1_prime[i], h[i], w_prime[i], DILI_N, 2 * param.gamma2);
    }

    memcpy(w1_prime_packed, mu, 64);
    offset = 64;
    for(i = 0; i < K; i++) {
        pack_w1(w1_prime_packed + offset, w1_prime[i], param.gamma2);
        offset += param.w1_packed_bytes;
    }

    shake256(c_tilde_prime, param.c_tilde_bytes, w1_prime_packed, offset);

    if(compare_bytes(c_tilde, c_tilde_prime, param.c_tilde_bytes) != 0) {
        return -1;
    }
    
    return 0;
}
