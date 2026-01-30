#include <stdio.h>
#include <stdint.h>

#define Q 3329
#define ROOT 17
#define MONT_R 2285

int64_t mod_exp(int64_t base, int64_t exp, int64_t mod) {
    int64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}

int64_t to_mont(int64_t a, int64_t R, int64_t q) {
    int64_t result = (a * R) % q;
    if (result < 0) result += q;
    if (result > q / 2) result -= q;
    return result;
}

int bit_reverse(int x) {
    int result = 0;
    for (int i = 0; i < 8; i++) {
        result <<= 1;
        result |= (x & 1);
        x >>= 1;
    }
    return result;
}

int main() {
    printf("Kyber q = %d, root = %d, R = %d\n\n", Q, ROOT, MONT_R);
    int64_t root_powers_mont[256];
    int64_t root_power = 1;
    for (int i = 0; i < 256; i++) {
        root_powers_mont[i] = to_mont(root_power, MONT_R, Q);
        root_power = (root_power * ROOT) % Q;
        if (root_power < 0) root_power += Q;
    }

    printf("static const int16_t zetas_kyber_256[256] = {\n");

    for (int i = 0; i < 256; i++) {
        if (i % 8 == 0) {
            printf("    ");
        }

        int64_t value;
        if (i == 0) {
            value = 0;
        } else {
            int br_index = bit_reverse(i);
            value = root_powers_mont[br_index];
        }

        printf("%5lld", value);

        if (i < 255) {
            printf(",");
            if (i % 8 != 7) {
                printf(" ");
            }
        }

        if (i % 8 == 7 || i == 255) {
            printf("\n");
        }
    }

    printf("};\n");
    return 0;
}
