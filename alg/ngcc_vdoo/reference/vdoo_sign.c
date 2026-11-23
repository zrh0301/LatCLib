#include "ngcc/auxfunc.h"
#include "vdoo_sign.h"
#include "vdoo_keypair.h"
#include "gf_config.h"
#include <string.h>

static int solve_diagonal_eq(uint8_t *y, int diag_eq_idx,
                              const uint8_t *t, const sk_t *sk)
{
    int i = diag_eq_idx;
    int diag_var_idx = VDOO_V1 + i;

    uint8_t A_coeff = 0;
    uint8_t B_const = 0;

    // VD, linear terms (A)
    for (int j = 0; j < diag_var_idx; j++)
    {
        // Index of the monomial t_j * t_diag_var_idx
        int monom_idx = j * VDOO_N + diag_var_idx - j * (j + 1) / 2;
        uint8_t coeff = gfv_get_ele(sk->F, monom_idx * VDOO_M + i);
        if (coeff != 0)
            A_coeff ^= gfv_mul(coeff, gfv_get_ele(y, j));
    }

    // VV, constant terms (B)
    for (int j = 0; j < diag_var_idx; j++)
    {
        for (int k = j; k < diag_var_idx; k++)
        {
            // Index of the monomial t_j * t_k
            int monom_idx = j * VDOO_N + k - j * (j + 1) / 2;
            uint8_t coeff = gfv_get_ele(sk->F, monom_idx * VDOO_M + i);
            if (coeff != 0)
                B_const ^= gfv_mul(coeff, gfv_mul(gfv_get_ele(y, j), gfv_get_ele(y, k)));
        }
    }

    // A * t_diag + B = y_i
    uint8_t yi = gfv_get_ele(t, i);
    // if A = 0, B = y_i
    if (A_coeff == 0)
    {
        // if B != y_i, it fails, we retry
        if (B_const != yi)
            return -1;

        // if B = y_i, t_diag can have a random value
        uint8_t val;

        get_random_number(&pqc_ngcc_vdoo_drng, &val, 1 * 8ULL);
        if (VDOO_Q == 16)
            val &= 0x0f;

        gfv_set_ele(y, diag_var_idx, val);
    }
    else
        // If A =/= 0, t_diag = (yi - B) * A⁻¹
        // In GF(2^l), substraction <=> addition <=> XOR
        // So, t_diag = (yi ^ B) * A⁻¹
        gfv_set_ele(y, diag_var_idx, gfv_mul(yi ^ B_const, gfv_inv(A_coeff)));

    return 0;
}

static void substitute_fixed_vars(uint8_t *A_oil, uint8_t *b_oil,
                                   const uint8_t *y, const uint8_t *t,
                                   int eq_start, int oil_start,
                                   int n_oil, int n_fixed, const sk_t *sk)
{
    for (int j = 0; j < n_oil; j++)
    {
        int eq = eq_start + j;

        // linear coeffs, construction of the line j of the matrix A
        // XOR for each k of coeff(i, oil_k) * x_i*x_oil_k
        //                 = coeff(i, oil_k) * y[i]*x_oil_k (y[i] known)
        for (int k = 0; k < n_oil; k++)
        {
            int oil_k = oil_start + k;
            uint8_t coeff = 0;
            for (int i = 0; i < n_fixed; i++)
            {
                int monom_idx = i * VDOO_N + oil_k - i * (i + 1) / 2;
                uint8_t c = gfv_get_ele(sk->F, monom_idx * VDOO_M + eq);
                if (c != 0)
                    coeff ^= gfv_mul(c, gfv_get_ele(y, i));
            }
            A_oil[j * n_oil + k] = coeff;
        }

        // constant term b[j]
        uint8_t constant = 0;
        for (int i = 0; i < n_fixed; i++)
        {
            for (int k = i; k < n_fixed; k++)
            {
                // XOR for each k of coeff(i, k) * x_i*x_k
                //                 = coeff(i, k) * y[i]*y[k] (y[i] and y[ k] known)
                int monom_idx = i * VDOO_N + k - i * (i + 1) / 2;
                uint8_t c = gfv_get_ele(sk->F, monom_idx * VDOO_M + eq);
                if (c != 0)
                    constant ^= gfv_mul(c,
                                    gfv_mul(gfv_get_ele(y, i),
                                            gfv_get_ele(y, k)));
            }
        }
        // b[i] = t[i] - C = t[i] + C = t[i] ^ C, in GF(2^l)
        b_oil[j] = gfv_get_ele(t, eq) ^ constant;
    }
}

static int solve_oil_layer(uint8_t *y, size_t offset,
                           const uint8_t *matrix, const uint8_t *rhs, size_t n)
{
    uint8_t work[n * (n + 1)];
    uint8_t solution[n];
    if (gf_composite_matrix_solve(solution, matrix, rhs, n, work, VDOO_FIELD) != 0)
        return -1;
    for (size_t i = 0; i < n; ++i)
        gfv_set_ele(y, offset + i, solution[i]);
    return 0;
}

// VDOOCentPoly_Inversion
static int vdoo_inverse_F(uint8_t *y, const sk_t *sk, const uint8_t *t)
{
    uint8_t A1[VDOO_O1 * VDOO_O1];
    uint8_t b1[VDOO_O1];
    uint8_t A2[VDOO_O2 * VDOO_O2];
    uint8_t b2[VDOO_O2];

    int max_attemps = 100;
    for (int attempt = 0; attempt < max_attemps; attempt++)
    {
        memset(y, 0, VDOO_N_BYTE);
        for (int i = 0; i < VDOO_V1; i++)
        {
            uint8_t val;
            get_random_number(&pqc_ngcc_vdoo_drng, &val, 1 * 8ULL);
            if (VDOO_Q == 16) val
                &= 0x0f;
            gfv_set_ele(y, i, val);
        }

        // Resolution of the diagonal layer
        int ok = 1;
        for (int i = 0; i < VDOO_D && ok; i++)
            if (solve_diagonal_eq(y, i, t, sk) != 0)
                ok = 0;
        if (!ok)
            continue;


        memset(A1, 0, sizeof(A1));
        memset(b1, 0, sizeof(b1));

        // Substitution
        substitute_fixed_vars(A1, b1, y, t,
                              VDOO_D, VDOO_V2, VDOO_O1, VDOO_V2, sk);

        // Gaussian Elimination
        if (solve_oil_layer(y, VDOO_V2, A1, b1, VDOO_O1) != 0)
            continue;


        memset(A2, 0, sizeof(A2));
        memset(b2, 0, sizeof(b2));

        // Substitution
        substitute_fixed_vars(A2, b2, y, t,
                              VDOO_D + VDOO_O1, VDOO_V3, VDOO_O2, VDOO_V3, sk);

        // Gaussian Elimination
        if (solve_oil_layer(y, VDOO_V3, A2, b2, VDOO_O2) != 0)
            continue;

        return 0;
    }
    return -1;
}

int vdoo_sign(uint8_t *signature, const sk_t *sk, const uint8_t *_digest)
{
    uint8_t salt[SALT_BYTES];
    uint8_t digest_salt[HASH_LEN + SALT_BYTES];

    uint8_t d[VDOO_M_BYTE];
    uint8_t t[VDOO_M_BYTE];

    uint8_t y[VDOO_N_BYTE];
    uint8_t s[VDOO_N_BYTE];

    int max_attempts = 100;
    for (int attempt = 0; attempt < max_attempts; attempt++)
    {
        get_random_number(&pqc_ngcc_vdoo_drng, salt, SALT_BYTES * 8ULL);

        memcpy(digest_salt, _digest, HASH_LEN);
        memcpy(digest_salt + HASH_LEN, salt, SALT_BYTES);
        if (pseudoXOF(VDOO_M_BYTE * 8, digest_salt, sizeof(digest_salt) * 8, d) != 0)
            return -1;

        // t = invS * d
        gf_composite_col_matvec(t, sk->invS, d, VDOO_M, VDOO_M, VDOO_M_BYTE, VDOO_FIELD);

        // y = F⁻¹(t)
        if (vdoo_inverse_F(y, sk, t) != 0)
            continue;

        // s = invT * y
        gf_composite_col_matvec(s, sk->invT, y, VDOO_N, VDOO_N, VDOO_N_BYTE, VDOO_FIELD);

        // make signature (s, salt)
        memcpy(signature, s, VDOO_N_BYTE);
        memcpy(signature + VDOO_N_BYTE, salt, SALT_BYTES);

        return 0;
    }

    return -1;
}

