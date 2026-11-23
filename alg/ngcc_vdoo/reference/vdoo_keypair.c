#include "vdoo_keypair.h"

#include <time.h>

#include "gf/byte_matrix.h"

static void generate_M_invM(unsigned char *M, unsigned char *invM, int n)
{
    int nb_tri = SUM_K(n);

    unsigned char L[nb_tri];
    unsigned char U[nb_tri];

    // L (lower triangular matrix with diagonal elements != 0)
    //     [a 0 0]
    // L = [b c 0]
    //     [d e f]
    // stocked as an array: L = [a, b, c, d, e, f]
    int idx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            uint8_t val;
            do
            {
                get_random_number(&pqc_ngcc_vdoo_drng, &val, 1 * 8ULL);
                if (VDOO_Q == 16)
                    val &= 0x0f;
            } while (i == j && val == 0);

            L[idx++] = val;
        }
    }

    // U (upper triangular matriw with diagonal elements != 0)
    //     [a' b' c']
    // U = [0  d' e']
    //     [0  0  f']
    // stocked as an array: U = [a', b', c', d', e', f']
    idx = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            uint8_t val;
            do
            {
                get_random_number(&pqc_ngcc_vdoo_drng, &val, 1 * 8ULL);
                if (VDOO_Q == 16)
                    val &= 0x0f;
            } while (i == j && val == 0);

            U[idx++] = val;
        }
    }

    // M = L * U (invertible matrix thanks to L and U structure)
    // Stocked in column-major, to use with batch operations :
    //     [x00 x01 x02]
    // M = [x10 x11 x12]
    //     [x20 x21 x22]
    // M = [x00, x10, x20, x01, x11, x21, x02, x12, x22]
    gf_composite_triangular_product(M, L, U, n, 0, VDOO_FIELD);

    // invL and invU
    unsigned char invL[nb_tri];
    unsigned char invU[nb_tri];

    gf_composite_triangular_inv(invL, L, n, 0, VDOO_FIELD);
    gf_composite_triangular_inv(invU, U, n, 1, VDOO_FIELD);

    // invM = invU * invL
    // stocked in column-major
    gf_composite_triangular_product(invM, invL, invU, n, 1, VDOO_FIELD);
}

static int is_allowed_diag(int poly, int i, int j)
{
    int diag_var_idx = VDOO_V1 + poly;

    // Vinegar × Vinegar
    if (i < diag_var_idx && j < diag_var_idx)
        return 1;

    // Vinegar × Diagonal
    if (i < diag_var_idx && j == diag_var_idx)
        return 1;

    return 0;
}

static int is_allowed_oil1(int poly, int i, int j)
{
    int oil_var_idx = VDOO_V2 + poly;
    int vinegar_max = VDOO_V2;

    // Vinegar × Vinegar
    if (i < vinegar_max && j < vinegar_max)
        return 1;

    // Vinegar × Oil
    if (i < vinegar_max && j == oil_var_idx)
        return 1;

    return 0;
}

static int is_allowed_oil2(int poly, int i, int j)
{
    int oil_var_idx = VDOO_V3 + poly;
    int vinegar_max = VDOO_V3;

    // Vinegar × Vinegar
    if (i < vinegar_max && j < vinegar_max)
        return 1;

    // Vinegar × Oil
    if (i < vinegar_max && j == oil_var_idx)
        return 1;

    return 0;
}

static int is_term_allowed(int poly, int i, int j)
{
    if (poly < VDOO_D)
        return is_allowed_diag(poly, i, j);
    if (poly < VDOO_D + VDOO_O1)
        return is_allowed_oil1(poly - VDOO_D, i, j);
    if (poly < VDOO_D + VDOO_O1 + VDOO_O2)
        return is_allowed_oil2(poly - VDOO_D - VDOO_O1, i, j);

    return -1;
}

static void set_coeff_to_zero(unsigned char *F, int poly, int monom_idx, int m)
{
    int coeff_idx = monom_idx * m + poly;
    gfv_set_ele(F, coeff_idx, 0);
}

static void set_zero_forbidden_terms(unsigned char *F)
{
    for (int poly = 0; poly < VDOO_M; poly++)
    {
        for (int i = 0; i < VDOO_N; i++)
        {
            for (int j = i; j < VDOO_N; j++)
            {
                int monom_idx = i * VDOO_N - i*(i-1)/2 + j - i;
                if (!is_term_allowed(poly, i, j))
                    set_coeff_to_zero(F, poly, monom_idx, VDOO_M);
            }
        }
    }
}

static void generate_central_map(sk_t *sk)
{
    get_random_number(&pqc_ngcc_vdoo_drng, sk->F, PK_BYTES * 8ULL);
    set_zero_forbidden_terms(sk->F);
}


static void generate_private_key(sk_t *sk, const unsigned char *sk_seed, unsigned char *S, unsigned char *T)
{
    memcpy(sk->sk_seed, sk_seed, LEN_SKSEED);

    // S and invS
    unsigned char invS[VDOO_M_BYTE * VDOO_M];
    generate_M_invM(S, invS, VDOO_M);

    // T and invT
    unsigned char invT[VDOO_N_BYTE * VDOO_N];
    generate_M_invM(T, invT, VDOO_N);

    // Only invS and invT are stocked, as T and S are only use to generate the central map
    memcpy(sk->invS, invS, VDOO_M_BYTE * VDOO_M);
    memcpy(sk->invT, invT, VDOO_N_BYTE * VDOO_N);

    generate_central_map(sk);

    memset(invS, 0, VDOO_M_BYTE * VDOO_M);
    memset(invT, 0, VDOO_N_BYTE * VDOO_N);
}

// each polynomial can be written such as:
// f_k(T*w) = (T*w)^T A_k (T*w)
//        = (w^T * T^T) A_k (T * w)
//        = w^T (T^T * A_k * T) w
// Thus, we need to compute A_k' = (T^T * A_k * T) for each polynomial in
// order to find A' the matrix with the coefficients of the monomials of F(T*w).
// We use batch matrix functions to compute A_k' for every polynomial at the same time.
void apply_T_to_F(unsigned char *FT, const unsigned char *F,
                          const unsigned char *T)
{
    static unsigned char tmp[VDOO_N * VDOO_N * VDOO_M_BYTE];
    static unsigned char A_prime[VDOO_N * VDOO_N * VDOO_M_BYTE];

    memset(tmp, 0, sizeof(tmp));
    memset(A_prime, 0, sizeof(A_prime));

    // tmp = F * T
    gf_composite_batch_tri_madd(tmp, F, T, VDOO_N, VDOO_N_BYTE, VDOO_N,
                                VDOO_M_BYTE, VDOO_FIELD);
    // A' = T^T * tmp
    gf_composite_batch_transpose_madd(A_prime, T, VDOO_N, VDOO_N_BYTE, VDOO_N,
                                      tmp, VDOO_N, VDOO_M_BYTE, VDOO_FIELD);

    // "Trianglize" the obtained matrix A'
    static unsigned char A_prime_tri[VDOO_M_BYTE * SUM_K(VDOO_N)];
    memset(A_prime_tri, 0, sizeof(A_prime_tri));
    gf_byte_batch_fold_upper_add(A_prime_tri, A_prime, VDOO_N, VDOO_M_BYTE);
    memcpy(FT, A_prime_tri, VDOO_M_BYTE * SUM_K(VDOO_N));

    memset(tmp, 0, VDOO_N * VDOO_N * VDOO_M_BYTE);
    memset(A_prime, 0, VDOO_N * VDOO_N * VDOO_M_BYTE);
}

static void apply_S_to_FT(uint8_t *out, const uint8_t *ft, const uint8_t *s)
{
    for (size_t k = 0; k < SUM_K(VDOO_N); ++k)
        gf_composite_col_matvec(out + k * VDOO_M_BYTE, s, ft + k * VDOO_M_BYTE,
                                VDOO_M, VDOO_M, VDOO_M_BYTE, VDOO_FIELD);
}

static void generate_public_key(pk_t *pk, const sk_t *sk, unsigned char *S, unsigned char *T)
{
    static unsigned char FT[PK_BYTES];

    // FT = F ° T
    // FT w = F(T w)
    apply_T_to_F(FT, sk->F, T);

    // P = S ° FT (= S ° F ° T)
    // P w = S(FT w)
    apply_S_to_FT(pk->pk, FT, S);

    memset(FT, 0, PK_BYTES);
}

int generate_keypair(pk_t *pk, sk_t *sk, const unsigned char *sk_seed)
{
    if (!pk || !sk || !sk_seed)
        return -1;

    unsigned char S[VDOO_M_BYTE * VDOO_M];
    unsigned char T[VDOO_N_BYTE * VDOO_N];

    generate_private_key(sk, sk_seed, S, T);
    generate_public_key(pk, sk, S, T);

    memset(S, 0, VDOO_M_BYTE * VDOO_M);
    memset(T, 0, VDOO_N_BYTE * VDOO_N);

    return 0;
}