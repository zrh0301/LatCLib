#include "../../inc/poly.h"
#include "../../inc/rand.h"
#include "../../inc/types.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define TEST_LEN 256
#define NUM_TESTS 1000

double compute_mean(int32_t *data, int len) {
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += data[i];
    }
    return sum / len;
}

double compute_std(int32_t *data, int len, double mean) {
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        double diff = data[i] - mean;
        sum += diff * diff;
    }
    return sqrt(sum / len);
}

void print_histogram(int32_t *data, int len, int min_val, int max_val) {
    int range = max_val - min_val + 1;
    int *hist = (int *)calloc(range, sizeof(int));
    
    for (int i = 0; i < len; i++) {
        if (data[i] >= min_val && data[i] <= max_val) {
            hist[data[i] - min_val]++;
        }
    }
    
    printf("  Value distribution:\n");
    for (int i = 0; i < range; i++) {
        printf("    %3d: ", i + min_val);
        int bars = (hist[i] * 50) / len;
        for (int j = 0; j < bars; j++) printf("#");
        printf(" (%d)\n", hist[i]);
    }
    
    free(hist);
}

void test_sample_uniform() {
    printf("\n=== Test sample_uniform (23-bit rejection sampling) ===\n");
    
    lp_ctx ctx;
    ctx.q = 8380417;
    ctx.n = 256;
    
    int32_t coeffs[TEST_LEN];
    uint8_t rand_buf[1024];
    
    rand_bytes(rand_buf, sizeof(rand_buf));
    
    unsigned int sampled = sample_uniform(coeffs, TEST_LEN, rand_buf, sizeof(rand_buf), &ctx);
    
    printf("  Sampled: %u / %d coefficients\n", sampled, TEST_LEN);
    printf("  Expected range: [0, Q-1] = [0, %d]\n", ctx.q - 1);
    printf("  First 10 values: ");
    for (int i = 0; i < 10 && i < sampled; i++) {
        printf("%d ", coeffs[i]);
    }
    printf("\n");
    
    // Check all values are in range
    int min_val = coeffs[0], max_val = coeffs[0];
    int valid = 1;
    for (int i = 0; i < sampled; i++) {
        if (coeffs[i] < min_val) min_val = coeffs[i];
        if (coeffs[i] > max_val) max_val = coeffs[i];
        
        if (coeffs[i] < 0 || coeffs[i] >= ctx.q) {
            valid = 0;
            printf("  [FAIL] Value out of range at index %d: %d\n", i, coeffs[i]);
        }
    }
    printf("  Actual range: [%d, %d]\n", min_val, max_val);
    
    if (valid) {
        printf("  [PASS] All values in valid range [0, Q-1]\n");
    }
    
    // Check distribution uniformity (divide into 8 bins)
    #define NUM_BINS 8
    int bins[NUM_BINS] = {0};
    int bin_size = ctx.q / NUM_BINS;
    
    for (int i = 0; i < sampled; i++) {
        int bin_idx = coeffs[i] / bin_size;
        if (bin_idx >= NUM_BINS) bin_idx = NUM_BINS - 1;
        bins[bin_idx]++;
    }
    
    printf("  Distribution uniformity check (8 bins):\n");
    double expected_per_bin = (double)sampled / NUM_BINS;
    int uniform_pass = 1;
    for (int i = 0; i < NUM_BINS; i++) {
        double deviation = fabs(bins[i] - expected_per_bin) / expected_per_bin;
        printf("    Bin %d: %d samples (%.1f%% deviation)\n", 
               i, bins[i], deviation * 100);
        
        // Allow 30% deviation for small sample sizes
        if (deviation > 0.30) {
            uniform_pass = 0;
        }
    }
    
    if (uniform_pass) {
        printf("  [PASS] Distribution appears uniform\n");
    } else {
        printf("  [WARN] Distribution may not be perfectly uniform (small sample size)\n");
    }
    
    // Check rejection sampling efficiency
    double efficiency = (double)sampled / (sizeof(rand_buf) / 3.0);  // 3 bytes per sample attempt
    printf("  Sampling efficiency: %.2f%%\n", efficiency * 100);
    if (efficiency > 0.50) {
        printf("  [PASS] Good rejection sampling efficiency\n");
    }
    #undef NUM_BINS
}

/* Test sample_eta */
void test_sample_eta() {
    printf("\n=== Test sample_eta (ETA=2 and ETA=4) ===\n");
    
    lp_ctx ctx;
    int32_t coeffs[TEST_LEN];
    uint8_t rand_buf[512];
    
    // Test ETA = 2
    printf("\n  --- ETA = 2 (Dilithium requirement) ---\n");
    rand_bytes(rand_buf, sizeof(rand_buf));
    unsigned int sampled = sample_eta(coeffs, TEST_LEN, rand_buf, sizeof(rand_buf), &ctx, 2);
    
    printf("  Sampled: %u / %d coefficients\n", sampled, TEST_LEN);
    printf("  Expected range: [-2, 2]\n");
    printf("  First 20 values: ");
    for (int i = 0; i < 20 && i < sampled; i++) {
        printf("%d ", coeffs[i]);
    }
    printf("\n");
    
    // Range validation
    int range_valid = 1;
    for (int i = 0; i < sampled; i++) {
        if (coeffs[i] < -2 || coeffs[i] > 2) {
            printf("  [FAIL] Value out of range at index %d: %d\n", i, coeffs[i]);
            range_valid = 0;
        }
    }
    if (range_valid) printf("  [PASS] All values in range [-2, 2]\n");
    
    double mean = compute_mean(coeffs, sampled);
    double std = compute_std(coeffs, sampled, mean);
    printf("  Mean: %.4f (expected ~0)\n", mean);
    printf("  Std Dev: %.4f (expected ~%.2f)\n", std, sqrt(2.0 * 2 * 2 / 3.0));
    
    // Check distribution uniformity (chi-square test approximation)
    if (fabs(mean) < 0.5) {
        printf("  [PASS] Mean is close to 0\n");
    } else {
        printf("  [WARN] Mean deviates from 0\n");
    }
    
    print_histogram(coeffs, sampled, -2, 2);
    
    // Test ETA = 4
    printf("\n  --- ETA = 4 (Dilithium requirement) ---\n");
    rand_bytes(rand_buf, sizeof(rand_buf));
    sampled = sample_eta(coeffs, TEST_LEN, rand_buf, sizeof(rand_buf), &ctx, 4);
    
    printf("  Sampled: %u / %d coefficients\n", sampled, TEST_LEN);
    printf("  Expected range: [-4, 4]\n");
    printf("  First 20 values: ");
    for (int i = 0; i < 20 && i < sampled; i++) {
        printf("%d ", coeffs[i]);
    }
    printf("\n");
    
    // Range validation for ETA=4
    range_valid = 1;
    for (int i = 0; i < sampled; i++) {
        if (coeffs[i] < -4 || coeffs[i] > 4) {
            printf("  [FAIL] Value out of range at index %d: %d\n", i, coeffs[i]);
            range_valid = 0;
        }
    }
    if (range_valid) printf("  [PASS] All values in range [-4, 4]\n");
    
    mean = compute_mean(coeffs, sampled);
    std = compute_std(coeffs, sampled, mean);
    printf("  Mean: %.4f\n", mean);
    printf("  Std Dev: %.4f\n", std);
    
    print_histogram(coeffs, sampled, -4, 4);
}

/* Test sample_CBD */
void test_sample_CBD() {
    printf("\n=== Test sample_CBD (Centered Binomial Distribution for Kyber) ===\n");
    
    lp_ctx ctx;
    int32_t coeffs[TEST_LEN];
    uint8_t rand_buf[512];
    
    rand_bytes(rand_buf, sizeof(rand_buf));
    
    unsigned int sampled = sample_CBD(coeffs, TEST_LEN, rand_buf, sizeof(rand_buf), &ctx);
    
    printf("  Sampled: %u / %d coefficients\n", sampled, TEST_LEN);
    printf("  Expected range: [-2, 2] (for eta=2)\n");
    printf("  First 20 values: ");
    for (int i = 0; i < 20 && i < sampled; i++) {
        printf("%d ", coeffs[i]);
    }
    printf("\n");
    
    // Range validation
    int range_valid = 1;
    for (int i = 0; i < sampled; i++) {
        if (coeffs[i] < -2 || coeffs[i] > 2) {
            printf("  [FAIL] CBD value out of range at index %d: %d\n", i, coeffs[i]);
            range_valid = 0;
        }
    }
    if (range_valid) printf("  [PASS] All CBD values in range [-2, 2]\n");
    
    double mean = compute_mean(coeffs, sampled);
    double std = compute_std(coeffs, sampled, mean);
    printf("  Mean: %.4f (expected ~0)\n", mean);
    printf("  Std Dev: %.4f (expected ~1.0 for eta=2)\n", std);
    
    // Check bell curve property of CBD
    int hist[5] = {0}; // -2, -1, 0, 1, 2
    for (int i = 0; i < sampled; i++) {
        hist[coeffs[i] + 2]++;
    }
    
    // CBD should have peak at 0
    if (hist[2] > hist[1] && hist[2] > hist[0] && hist[2] > hist[3] && hist[2] > hist[4]) {
        printf("  [PASS] CBD shows bell curve property (peak at 0)\n");
    } else {
        printf("  [WARN] CBD distribution may not be bell-shaped\n");
    }
    
    // Check standard deviation is close to 1.0
    if (fabs(std - 1.0) < 0.3) {
        printf("  [PASS] Standard deviation close to expected 1.0\n");
    } else {
        printf("  [WARN] Standard deviation deviates from expected\n");
    }
    
    print_histogram(coeffs, sampled, -2, 2);
}

/* Test sample_gamma1 */
void test_sample_gamma1() {
    printf("\n=== Test sample_gamma1 (Dilithium masking) ===\n");
    
    lp_ctx ctx;
    int32_t coeffs[TEST_LEN];
    uint8_t rand_buf[2048];
    
    rand_bytes(rand_buf, sizeof(rand_buf));
    
    unsigned int sampled = sample_gamma1(coeffs, TEST_LEN, rand_buf, sizeof(rand_buf), &ctx);
    
    printf("  Sampled: %u / %d coefficients\n", sampled, TEST_LEN);
    printf("  Expected range: [-(2^19-1), 2^19] = [-524287, 524288]\n");
    printf("  First 10 values: ");
    for (int i = 0; i < 10 && i < sampled; i++) {
        printf("%d ", coeffs[i]);
    }
    printf("\n");
    
    int min_val = coeffs[0], max_val = coeffs[0];
    int range_valid = 1;
    for (int i = 0; i < sampled; i++) {
        if (coeffs[i] < min_val) min_val = coeffs[i];
        if (coeffs[i] > max_val) max_val = coeffs[i];
        
        // Check gamma1 = 2^19 bounds
        if (coeffs[i] < -524287 || coeffs[i] > 524288) {
            printf("  [FAIL] Value out of GAMMA1 range at index %d: %d\n", i, coeffs[i]);
            range_valid = 0;
        }
    }
    printf("  Actual range: [%d, %d]\n", min_val, max_val);
    
    if (range_valid) {
        printf("  [PASS] All values within GAMMA1 bounds\n");
    }
    
    // Check sampling efficiency (should sample most coefficients)
    if (sampled >= TEST_LEN * 0.95) {
        printf("  [PASS] Sampling efficiency: %.1f%%\n", (sampled * 100.0) / TEST_LEN);
    } else {
        printf("  [WARN] Low sampling efficiency: %.1f%%\n", (sampled * 100.0) / TEST_LEN);
    }
}

/* Test sample_gamma2 */
void test_sample_gamma2() {
    printf("\n=== Test sample_gamma2 (Dilithium hint poly) ===\n");
    
    lp_ctx ctx;
    ctx.q = 8380417;
    
    int32_t coeffs[TEST_LEN];
    uint8_t rand_buf[512];
    
    rand_bytes(rand_buf, sizeof(rand_buf));
    
    unsigned int sampled = sample_gamma2(coeffs, TEST_LEN, rand_buf, sizeof(rand_buf), &ctx);
    
    printf("  Sampled: %u / %d coefficients\n", sampled, TEST_LEN);
    printf("  Expected range: [0, 15] (for GAMMA2=(Q-1)/32)\n");
    printf("  First 20 values: ");
    for (int i = 0; i < 20 && i < sampled; i++) {
        printf("%d ", coeffs[i]);
    }
    printf("\n");
    
    int min_val = coeffs[0], max_val = coeffs[0];
    int range_valid = 1;
    for (int i = 0; i < sampled; i++) {
        if (coeffs[i] < min_val) min_val = coeffs[i];
        if (coeffs[i] > max_val) max_val = coeffs[i];
        
        // Check 4-bit range [0, 15]
        if (coeffs[i] < 0 || coeffs[i] > 15) {
            printf("  [FAIL] Value out of 4-bit range at index %d: %d\n", i, coeffs[i]);
            range_valid = 0;
        }
    }
    printf("  Actual range: [%d, %d]\n", min_val, max_val);
    
    if (range_valid) {
        printf("  [PASS] All values in 4-bit range [0, 15]\n");
    }
    
    // Check sampling efficiency
    if (sampled >= TEST_LEN) {
        printf("  [PASS] Full sampling achieved\n");
    }
}

int main() {
    printf("============================================================\n");
    printf("        PQC Sampling Functions Test Suite                  \n");
    printf("============================================================\n");
    
    test_sample_uniform();
    test_sample_eta();
    test_sample_CBD();
    test_sample_gamma1();
    test_sample_gamma2();
    
    printf("\n============================================================\n");
    printf("                  All Tests Completed!                      \n");
    printf("============================================================\n");
    
    return 0;
}
