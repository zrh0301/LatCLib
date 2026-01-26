#include "../inc/dili.h"
#include "../inc/poly.h"
#include "../inc/rand.h"
#include <string.h>

// ============================================
//  Helper Functions
// ============================================

// Power2Round: decompose r into (r1, r0) where r = r1*2^d + r0
static void power2round(lssi_t *r1, lssi_t *r0, const lssi_t *r, int n, int d) {
    // TODO: Implement power2round
    // r1[i] = (r[i] + (1 << (d-1))) >> d
    // r0[i] = r[i] - (r1[i] << d)
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
static void expand_a(lssi_t **A, const uint8_t *rho, int k, int l) {
    // TODO: Use SHAKE128 to expand rho into matrix A
    // A is k x l matrix of polynomials
}

// ExpandS: expand seed into secret vectors
static void expand_s(lssi_t *s, const uint8_t *seed, int len, int eta) {
    // TODO: Sample s from seed with coefficients in [-eta, eta]
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
void dili_keygen(dili_keypair_t* keypair) {
    // TODO: Implement key generation
    // 1. Generate random seed rho
    // 2. Expand A from rho
    // 3. Sample s1, s2 from uniform distribution
    // 4. Compute t = A*s1 + s2 (using NTT)
    // 5. Decompose t into t1 (high) and t0 (low)
    // 6. Store pk = (rho, t1), sk = (rho, s1, s2, t0)
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
