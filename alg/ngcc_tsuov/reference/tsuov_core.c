#include <string.h>
#include <stdlib.h>
#include "tsuov_core.h"
#include "util/bytes.h"
#include "util/secure_clear.h"

/* =====================================================================
   F_{q^L} linear algebra
   ---------------------------------------------------------------------
   The field, vector and planar-tower matrix operators live in the library
   (inc/gf/gf31.h and inc/gf/gf31_matrix.h). Only the TSUOV-specific integer
   inner products used by the E operator remain here.
   ===================================================================== */

static uint64_t vector_V_dot_vector_V(const vector_V A, const vector_V B)
{
    uint64_t C = 0;
    for (int j = 0; j < TSUOV_V; j++)
        C += (uint64_t)A[j] * B[j];
    return C;
}

static uint64_t vector_v_dot_vector_v(const VECTOR_V A, const VECTOR_V B)
{
    uint64_t ci = 0;
    for (int l = 0; l < TSUOV_L; l++)
        ci += vector_V_dot_vector_V(A[l], B[l]);
    return ci;
}

static void MATRIX_TRANSPOSE_m2xLO(const Fq A[TSUOV_m1][TSUOV_k][TSUOV_o],
                                 Fq C[TSUOV_k * TSUOV_o][TSUOV_m2])
{
    for (int i = 0; i < TSUOV_m1; i++)
        for (int k = 0; k < TSUOV_k; k++)
            for (int lo = 0; lo < TSUOV_o; lo++)
                C[k * TSUOV_o + lo][i] = A[i][k][lo];
    for (int k = 0; k < TSUOV_k; k++)
        for (int lo = 0; lo < TSUOV_o; lo++)
            memset(C[k * TSUOV_o + lo] + TSUOV_m1, 0, (TSUOV_m2 - TSUOV_m1) * sizeof(Fq));
}

static void whipk_accumulate_o_terms(vector_V acc, const MATRIX_OxV A, int l1,
                                     const MATRIX_kxO a, int i, int l2)
{
    vector_V tmp;
#if (TSUOV_O <= 8)
    for (int j = 0; j < TSUOV_O; j++) {
        gf31_vector_mul_scalar(A[j][l1], a[i][l2][j], BLOCK_ALIGNED_BYTES(TSUOV_V), tmp);
        gf31_vector_add_no_mod(acc, tmp, BLOCK_ALIGNED_BYTES(TSUOV_V), acc);
    }
    gf31_vector_reduce(acc, BLOCK_ALIGNED_BYTES(TSUOV_V), acc);
#else
    for (int j = 0; j < TSUOV_O; j++) {
        gf31_vector_mul_scalar(A[j][l1], a[i][l2][j], BLOCK_ALIGNED_BYTES(TSUOV_V), tmp);
        gf31_vector_add_mod(acc, tmp, BLOCK_ALIGNED_BYTES(TSUOV_V), acc);
    }
#endif
}

static void MATRIX_VxO_dot_VECTOR_O_whipk(const MATRIX_OxV A, const MATRIX_kxO a, MATRIX_kxV res)
{
    vector_V T[2 * TSUOV_L - 1];
    for (int i = 0; i < TSUOV_k; i++) {
        for (int l = 0; l < 2 * TSUOV_L - 1; l++)
            memset(T[l], 0, sizeof(T[l]));
        for (int l1 = 0; l1 < TSUOV_L; l1++)
            for (int l2 = 0; l2 < TSUOV_L; l2++)
                whipk_accumulate_o_terms(T[l1 + l2], A, l1, a, i, l2);
#if (TSUOV_L == 2)
        gf31_vector_mul_scalar(T[2], 3, BLOCK_ALIGNED_BYTES(TSUOV_V), T[2]);
        gf31_vector_add_mod(T[0], T[2], BLOCK_ALIGNED_BYTES(TSUOV_V), T[0]);
#else
#  error "MATRIX_VxO_dot_VECTOR_O_whipk: unsupported TSUOV_L"
#endif
        for (int l = 0; l < TSUOV_L; l++)
            memcpy(res[i][l], T[l], sizeof(T[l]));
    }
}

static void VECTOR_V_sub_VECTOR_V_whipk_save_in_sign(const MATRIX_kxV a, const MATRIX_kxV b,
                                                      MATRIX_kxN sign)
{
    for (int i = 0; i < TSUOV_k; i++)
        for (int l = 0; l < TSUOV_L; l++)
            gf31_vector_sub_mod(a[i][l], b[i][l], BLOCK_ALIGNED_BYTES(TSUOV_V), sign[i][l]);
}

/* =====================================================================
   scalar multiply_E (tables from root multiply_E.c)
   ===================================================================== */

#if (TSUOV_security_strength_category == 1)
static const uint8_t c0_T0_data[32] = {
    0,4,8,12,16,20,24,28,1,5,9,13,17,21,25,29,
    0,4,8,12,16,20,24,28,1,5,9,13,17,21,25,29};
static const uint8_t c0_T1_data[32] = {
    2,6,10,14,18,22,26,30,3,7,11,15,19,23,27,0,
    2,6,10,14,18,22,26,30,3,7,11,15,19,23,27,0};
static const uint8_t c1_T0_data[32] = {
    0,3,6,9,12,15,18,21,24,27,30,2,5,8,11,14,
    0,3,6,9,12,15,18,21,24,27,30,2,5,8,11,14};
static const uint8_t c1_T1_data[32] = {
    17,20,23,26,29,1,4,7,10,13,16,19,22,25,28,0,
    17,20,23,26,29,1,4,7,10,13,16,19,22,25,28,0};
static const uint8_t c2_T0_data[32] = {
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
static const uint8_t c2_T1_data[32] = {
    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,0,
    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,0};
#elif (TSUOV_security_strength_category == 3)
static const uint8_t c0_T0_data[32] = {
  0,6,12,18,24,30,5,11,17,23,29,4,10,16,22,28,
  0,6,12,18,24,30,5,11,17,23,29,4,10,16,22,28};
static const uint8_t c0_T1_data[32] = {
  3,9,15,21,27,2,8,14,20,26,1,7,13,19,25,0,
  3,9,15,21,27,2,8,14,20,26,1,7,13,19,25,0};
static const uint8_t c1_T0_data[32] = {
  0,4,8,12,16,20,24,28,1,5,9,13,17,21,25,29,
  0,4,8,12,16,20,24,28,1,5,9,13,17,21,25,29};
static const uint8_t c1_T1_data[32] = {
  2,6,10,14,18,22,26,30,3,7,11,15,19,23,27,0,
  2,6,10,14,18,22,26,30,3,7,11,15,19,23,27,0};
static const uint8_t c2_T0_data[32] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
static const uint8_t c2_T1_data[32] = {
  16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,0,
  16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,0};
#elif (TSUOV_security_strength_category == 5)
static const uint8_t c0_T0_data[32] = {
  0,15,30,14,29,13,28,12,27,11,26,10,25,9,24,8,
  0,15,30,14,29,13,28,12,27,11,26,10,25,9,24,8};
static const uint8_t c0_T1_data[32] = {
  23,7,22,6,21,5,20,4,19,3,18,2,17,1,16,0,
  23,7,22,6,21,5,20,4,19,3,18,2,17,1,16,0};
static const uint8_t c1_T0_data[32] = {
  0,12,24,5,17,29,10,22,3,15,27,8,20,1,13,25,
  0,12,24,5,17,29,10,22,3,15,27,8,20,1,13,25};
static const uint8_t c1_T1_data[32] = {
  6,18,30,11,23,4,16,28,9,21,2,14,26,7,19,0,
  6,18,30,11,23,4,16,28,9,21,2,14,26,7,19,0};
static const uint8_t c2_T0_data[32] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
static const uint8_t c2_T1_data[32] = {
  16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,0,
  16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,0};
#endif

static void init_mul_table_mod31(void) { }

static void multiply_E_vec_scalar(const int m, uint8_t *u, uint8_t *u_new, int ell)
{
    memset(u_new, 0, (size_t)ell);
    memcpy(u_new + ell, u, (size_t)(TSUOV_m2 - ell));

    for (int j = 0; j < ell + 2; j++) {
        int idx = m - ell + j;
        uint8_t u0 = (j < ell) ? u[idx] : 0;
        uint8_t u1 = (j >= 1 && j - 1 < ell) ? u[idx - 1] : 0;
        uint8_t u2 = (j >= 2 && j - 2 < ell) ? u[idx - 2] : 0;

        u_new[j] += (uint8_t)((u0 < 16) ? c0_T0_data[u0] : c0_T1_data[u0]);
        u_new[j] += (uint8_t)((u1 < 16) ? c1_T0_data[u1] : c1_T1_data[u1]);
        u_new[j] += (uint8_t)((u2 < 16) ? c2_T0_data[u2] : c2_T1_data[u2]);
    }

    for (int j = 0; j < TSUOV_m2; j++)
        u_new[j] = (uint8_t)(u_new[j] % 31);
}

static void multiply_E_add_m1_m2(Fq *u, Fq *u_acc, int ell)
{
    Fq u_m2[TSUOV_m2];
    memcpy(u_m2, u, TSUOV_m1);
    memset(u_m2 + TSUOV_m1, 0, (size_t)(TSUOV_m2 - TSUOV_m1) * sizeof(Fq));
    Fq u_new_temp[TSUOV_m2];
    multiply_E_vec_scalar(TSUOV_m2, u_m2, u_new_temp, ell);
    gf31_vector_add(u_acc, u_new_temp, TSUOV_m2);
}

static void multiply_E_mat_m2(const uint8_t *M, uint8_t *M_new, int ctr, int m_col)
{
    int ell = ctr;
    int shift = TSUOV_m2 - ell;

    for (int col = 0; col < m_col; col++) {
        const uint8_t *u = M + col * TSUOV_m2;
        uint8_t *u_new = M_new + col * TSUOV_m2;

        memset(u_new, 0, (size_t)ell);
        memcpy(u_new + ell, u, (size_t)(TSUOV_m2 - ell));

        for (int j = 0; j < ell + 2; j++) {
            uint8_t u0 = (j < ell) ? u[shift + j] : 0;
            uint8_t u1 = (j >= 1 && j - 1 < ell) ? u[shift + j - 1] : 0;
            uint8_t u2 = (j >= 2 && j - 2 < ell) ? u[shift + j - 2] : 0;

            u_new[j] += (uint8_t)((u0 < 16) ? c0_T0_data[u0] : c0_T1_data[u0]);
            u_new[j] += (uint8_t)((u1 < 16) ? c1_T0_data[u1] : c1_T1_data[u1]);
            u_new[j] += (uint8_t)((u2 < 16) ? c2_T0_data[u2] : c2_T1_data[u2]);
        }

        for (int j = 0; j < TSUOV_m2; j++)
            u_new[j] = (uint8_t)(u_new[j] % 31);
    }
}

static void multiply_E_add_mat_m2(Fq *M, Fq *M_new, int ctr, int m_col)
{
    Fq temp[TSUOV_o][TSUOV_m2];
    multiply_E_mat_m2((const uint8_t *)M, (uint8_t *)temp, ctr, m_col);
    gf31_vector_add((uint8_t *)M_new, (uint8_t *)temp, TSUOV_o * TSUOV_m2);
}

/* =====================================================================
   PRG / sampling (from root tsuov.c)
   ===================================================================== */

typedef struct {
    const uint8_t *src;
    const uint8_t *end;
    uint32_t       bit_buf;
    uint8_t        bit_cnt;
} rejsamp_pack_reader;

static inline void rejsamp_pack_reader_init(rejsamp_pack_reader *rd,
                                            const uint8_t *src, unsigned int tau)
{
    rd->src = src;
    rd->end = src + tau;
    rd->bit_buf = 0;
    rd->bit_cnt = 0;
}

static inline int rejsamp_pack_reader_refill(rejsamp_pack_reader *rd)
{
    if (rd->src >= rd->end)
        return -1;
    rd->bit_buf |= ((uint32_t)*rd->src++) << rd->bit_cnt;
    rd->bit_cnt += 8;
    return 0;
}

static inline int rejsamp_pack_reader_next(rejsamp_pack_reader *rd, uint8_t *out)
{
    for (;;) {
        while (rd->bit_cnt < 5) {
            if (rejsamp_pack_reader_refill(rd) != 0)
                return -1;
        }
        uint8_t v = (uint8_t)(rd->bit_buf & 0x1Fu);
        rd->bit_buf >>= 5;
        rd->bit_cnt = (uint8_t)(rd->bit_cnt - 5);
        if (v != TSUOV_q) {
            *out = v;
            return 0;
        }
    }
}

static inline void RejSamp(const unsigned int length, const unsigned int tau, uint8_t *dst)
{
    uint8_t xof[tau];
    memcpy(xof, dst, tau);
    rejsamp_pack_reader rd;
    rejsamp_pack_reader_init(&rd, xof, tau);
    for (unsigned int i = 0; i < length; i++) {
        if (rejsamp_pack_reader_next(&rd, &dst[i]) != 0)
            dst[i] = 0;
    }
    memset(xof, 0, tau);
}

static void secure_zero(void *ptr, size_t len)
{
    pqc_secure_clear(ptr, len);
}

MGF_CTX_s *PRG_init_sm3(const uint8_t seed[TSUOV_SEED_LEN])
{
    MGF_CTX_s *ctx = (MGF_CTX_s *)malloc(sizeof(MGF_CTX));
    if (ctx == NULL)
        return NULL;
    MGF_init(seed, TSUOV_SEED_LEN, ctx);
    return ctx;
}

void PRG_yield_sm3(MGF_CTX_s *ctx, int length, uint8_t *dst)
{
    MGF_yield(ctx, dst, (size_t)length);
}

void PRG_final_sm3(MGF_CTX_s *ctx)
{
    if (ctx == NULL)
        return;
    MGF_final(ctx);
    free(ctx);
}

MGF_CTX_s *PRG_copy_sm3(MGF_CTX_s *src)
{
    MGF_CTX_s *dst = (MGF_CTX_s *)malloc(sizeof(MGF_CTX));
    if (dst == NULL)
        return NULL;
    MGF_CTX_copy(src, dst);
    return dst;
}

int RejSampPRG_sm3(MGF_CTX_s *ctx, const uint64_t index, const unsigned int length,
                   const unsigned int tau, uint8_t *dst)
{
    uint8_t msg[TSUOV_SEED_LEN + 2];
    memcpy(msg, ctx->in, TSUOV_SEED_LEN);
    pqc_store_be(msg + TSUOV_SEED_LEN, 2, (uint64_t)index);
    /* sign-31-8: propagate a failed placeholder-XOF call instead of sampling
     * from an unwritten (possibly stale) buffer. */
    if (pseudoXOF((unsigned long long)tau * 8, msg,
                  (unsigned long long)(TSUOV_SEED_LEN + 2) * 8, dst) != 0) {
        secure_zero(msg, sizeof(msg));
        return -1;
    }
    RejSamp(length, tau, dst);
    secure_zero(msg, sizeof(msg));
    return 0;
}

/* sign-31-3 and sign-31-5: detect the size overflow and report a failed
 * allocation instead of hashing a stale stack buffer. */
static int Expand_mu(const TSUOV_SEED seed_pk, const uint8_t message[],
                     const size_t message_length, uint8_t mu[TSUOV_MU_LEN])
{
    size_t buf_len;
    uint8_t *buf;

    if (message_length > ((SIZE_MAX >> 3) - TSUOV_SEED_LEN))
        return -1;
    buf_len = TSUOV_SEED_LEN + message_length;
    buf = (uint8_t *)malloc(buf_len);
    if (buf == NULL)
        return -1;
    memcpy(buf, seed_pk, TSUOV_SEED_LEN);
    if (message_length != 0)
        memcpy(buf + TSUOV_SEED_LEN, message, message_length);
    int rc = pseudohash(512, buf, (unsigned long long)buf_len * 8, mu);
    secure_zero(buf, buf_len);
    free(buf);
    return rc == 0 ? 0 : -1;
}

static int Hash(const uint8_t mu[TSUOV_MU_LEN], const TSUOV_SALT salt, uint8_t dst[TSUOV_m2])
{
    uint8_t tmp[TSUOV_m2 > TSUOV_t_RS_LEN ? TSUOV_m2 : TSUOV_t_RS_LEN];
    uint8_t msg[TSUOV_MU_LEN + TSUOV_SALT_LEN];
    memcpy(msg, mu, TSUOV_MU_LEN);
    memcpy(msg + TSUOV_MU_LEN, salt, TSUOV_SALT_LEN);
    int rc = pseudoXOF((unsigned long long)TSUOV_t_RS_LEN * 8, msg,
                       (unsigned long long)(TSUOV_MU_LEN + TSUOV_SALT_LEN) * 8, tmp);
    if (rc == 0) {
        RejSamp(TSUOV_m2, TSUOV_t_RS_LEN, tmp);
        memcpy(dst, tmp, TSUOV_m2);
    }
    secure_zero(tmp, sizeof(tmp));
    secure_zero(msg, sizeof(msg));
    return rc == 0 ? 0 : -1;
}

static void ExpandMatrixOxV(uint8_t *src, MATRIX_OxV A)
{
    uint8_t *s = src;
    for (int i = 0; i < TSUOV_O; i++) {
        for (int j = 0; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++)
                A[i][k][j] = *s++;
        VECTOR_V_CLEAR_TAIL(A[i]);
    }
}

static void ExpandMatrixkxV(uint8_t *src, MATRIX_kxV A)
{
    uint8_t *s = src;
    for (int i = 0; i < TSUOV_k; i++) {
        for (int j = 0; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++)
                A[i][k][j] = *s++;
        VECTOR_V_CLEAR_TAIL(A[i]);
    }
}

static void ExpandSymmetricMatrixVxV(uint8_t *src, MATRIX_VxV A)
{
    uint8_t *s = src;
    for (int i = 0; i < TSUOV_V; i++) {
        for (int j = i; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++) {
                A[i][k][j] = *s++;
                A[j][k][i] = A[i][k][j];
            }
        VECTOR_V_CLEAR_TAIL(A[i]);
    }
}

static int Expand_sk(const TSUOV_SEED seed_sk, MATRIX_OxV SdT)
{
    const int n2 = TSUOV_Pi2_LEN;
    uint8_t r2[TSUOV_Pi2_LEN];
    TSUOV_PRG_CTX *ctx = PRG_init(seed_sk);
    if (ctx == NULL)
        return -1;
    int rc = RejSampPRG(ctx, 0, n2, TSUOV_Pi2_RS_LEN, r2);
    PRG_final(ctx);
    if (rc != 0) {
        secure_zero(r2, TSUOV_Pi2_LEN);
        return -1;
    }
    uint8_t *ptr_r2 = r2;
    for (int i = 0; i < TSUOV_O; i++) {
        for (int j = 0; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++)
                SdT[i][k][j] = *ptr_r2++;
        VECTOR_V_CLEAR_TAIL(SdT[i]);
    }
    secure_zero(r2, TSUOV_Pi2_LEN);
    return 0;
}

static int Expand_S(const TSUOV_SEED seed_sk, MATRIX_OxN SdT)
{
    const int n2 = TSUOV_Pi2_LEN;
    uint8_t r2[TSUOV_Pi2_LEN];
    TSUOV_PRG_CTX *ctx = PRG_init(seed_sk);
    if (ctx == NULL)
        return -1;
    int rc = RejSampPRG(ctx, 0, n2, TSUOV_Pi2_RS_LEN, r2);
    PRG_final(ctx);
    if (rc != 0) {
        secure_zero(r2, TSUOV_Pi2_LEN);
        return -1;
    }

    uint8_t *ptr_r2 = r2;
    for (int i = 0; i < TSUOV_O; i++)
        for (int j = 0; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++)
                SdT[i][k][j] = *ptr_r2++;

    int N_PADDED = (TSUOV_N + 31) & ~31;
    for (int i = 0; i < TSUOV_O; i++)
        for (int k = 0; k < TSUOV_L; k++)
            memset(&SdT[i][k][TSUOV_V], 0, (N_PADDED - TSUOV_V) * sizeof(uint8_t));
    for (int i = 0; i < TSUOV_O; i++)
        SdT[i][0][TSUOV_V + i] = 1;
    secure_zero(r2, TSUOV_Pi2_LEN);
    return 0;
}

static int Expand_pk(TSUOV_PRG_CTX *ctx0, const uint64_t index,
                     MATRIX_VxV Pi1, MATRIX_OxV Pi2T)
{
    const int n1 = TSUOV_Pi1_LEN;
    const int n2 = TSUOV_Pi2_LEN;
    uint8_t r1[TSUOV_Pi1_LEN];
    uint8_t r2[TSUOV_Pi2_LEN];

    TSUOV_PRG_CTX *ctx1 = PRG_copy(ctx0);
    if (ctx1 == NULL)
        return -1;
    int rc = RejSampPRG(ctx1, 2 * index, n1, TSUOV_Pi1_RS_LEN, r1);
    if (rc == 0)
        ExpandSymmetricMatrixVxV(r1, Pi1);
    PRG_final(ctx1);
    if (rc != 0)
        return -1;

    TSUOV_PRG_CTX *ctx2 = PRG_copy(ctx0);
    if (ctx2 == NULL)
        return -1;
    rc = RejSampPRG(ctx2, 2 * index + 1, n2, TSUOV_Pi2_RS_LEN, r2);
    if (rc == 0)
        ExpandMatrixOxV(r2, Pi2T);
    PRG_final(ctx2);
    return rc;
}

static int Expand_PK(const Fq P3i[TSUOV_O][TSUOV_L][TSUOV_O], TSUOV_PRG_CTX *ctx0,
                     const uint64_t index, MATRIX_NxN Pi)
{
    const int n1 = TSUOV_Pi1_LEN;
    const int n2 = TSUOV_Pi2_LEN;
    uint8_t r1[TSUOV_Pi1_LEN];
    uint8_t r2[TSUOV_Pi2_LEN];

    TSUOV_PRG_CTX *ctx1 = PRG_copy(ctx0);
    if (ctx1 == NULL)
        return -1;
    int rc = RejSampPRG(ctx1, 2 * index, n1, TSUOV_Pi1_RS_LEN, r1);
    PRG_final(ctx1);
    if (rc != 0)
        return -1;

    TSUOV_PRG_CTX *ctx2 = PRG_copy(ctx0);
    if (ctx2 == NULL)
        return -1;
    rc = RejSampPRG(ctx2, 2 * index + 1, n2, TSUOV_Pi2_RS_LEN, r2);
    PRG_final(ctx2);
    if (rc != 0)
        return -1;

    uint8_t *ptr_r1 = r1;
    uint8_t *ptr_r2 = r2;

    for (int i = 0; i < TSUOV_V; i++)
        for (int j = i; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++) {
                uint8_t val = TSUOV_q - *ptr_r1++;
                Pi[i][k][j] = val;
                Pi[j][k][i] = val;
            }

    for (int i = TSUOV_V; i < TSUOV_N; i++)
        for (int j = 0; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++) {
                uint8_t val = *ptr_r2++;
                Pi[i][k][j] = val;
                Pi[j][k][i] = val;
            }

    for (int i = TSUOV_V; i < TSUOV_N; i++)
        for (int k = 0; k < TSUOV_L; k++)
            memcpy(&Pi[i][k][TSUOV_V], &P3i[i - TSUOV_V][k][0], TSUOV_O * sizeof(uint8_t));

    int N_PADDED = (TSUOV_N + 31) & ~31;
    if (N_PADDED == TSUOV_N)
        return 0;
    for (int i = 0; i < TSUOV_N; i++)
        for (int k = 0; k < TSUOV_L; k++)
            memset(&Pi[i][k][TSUOV_N], 0, (N_PADDED - TSUOV_N) * sizeof(uint8_t));
    return 0;
}

static int Expand_P_MATRIX(TSUOV_PRG_CTX *ctx0, const uint64_t index, MATRIX_NxN Pi)
{
    const int n1 = TSUOV_Pi1_LEN;
    const int n2 = TSUOV_Pi2_LEN;
    uint8_t r1[TSUOV_Pi1_LEN];
    uint8_t r2[TSUOV_Pi2_LEN];

    TSUOV_PRG_CTX *ctx1 = PRG_copy(ctx0);
    if (ctx1 == NULL)
        return -1;
    int rc = RejSampPRG(ctx1, 2 * index, n1, TSUOV_Pi1_RS_LEN, r1);
    PRG_final(ctx1);
    if (rc != 0)
        return -1;

    TSUOV_PRG_CTX *ctx2 = PRG_copy(ctx0);
    if (ctx2 == NULL)
        return -1;
    rc = RejSampPRG(ctx2, 2 * index + 1, n2, TSUOV_Pi2_RS_LEN, r2);
    PRG_final(ctx2);
    if (rc != 0)
        return -1;

    uint8_t *ptr_r1 = r1;
    uint8_t *ptr_r2 = r2;

    for (int i = 0; i < TSUOV_V; i++)
        for (int j = i; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++) {
                uint8_t val = *ptr_r1++;
                Pi[i][k][j] = val;
                Pi[j][k][i] = val;
            }

    for (int i = TSUOV_V; i < TSUOV_N; i++)
        for (int j = 0; j < TSUOV_V; j++)
            for (int k = 0; k < TSUOV_L; k++) {
                uint8_t val = *ptr_r2++;
                Pi[i][k][j] = val;
                Pi[j][k][i] = val;
            }

    for (int i = TSUOV_V; i < TSUOV_N; i++)
        for (int k = 0; k < TSUOV_L; k++)
            memset(&Pi[i][k][TSUOV_V], 0, (TSUOV_N - TSUOV_V) * sizeof(uint8_t));

    int N_PADDED = (TSUOV_N + 31) & ~31;
    if (N_PADDED == TSUOV_N)
        return 0;
    for (int i = 0; i < TSUOV_N; i++)
        for (int k = 0; k < TSUOV_L; k++)
            memset(&Pi[i][k][TSUOV_N], 0, (N_PADDED - TSUOV_N) * sizeof(uint8_t));
    return 0;
}

static int Expand_V(const TSUOV_SEED seed_V, MATRIX_kxV V)
{
    const int n3 = TSUOV_kv_LEN;
    uint8_t r3[TSUOV_kv_LEN];
    TSUOV_PRG_CTX *ctx = PRG_init(seed_V);
    if (ctx == NULL)
        return -1;
    int rc = RejSampPRG(ctx, 0, n3, TSUOV_kv_RS_LEN, r3);
    if (rc == 0)
        ExpandMatrixkxV(r3, V);
    PRG_final(ctx);
    return rc;
}

static int Expand_sol(const TSUOV_SEED seed_sol, uint8_t dst[TSUOV_m2])
{
    const int n4 = TSUOV_t_LEN;
    uint8_t r4[TSUOV_m2 > TSUOV_t_RS_LEN ? TSUOV_m2 : TSUOV_t_RS_LEN];
    TSUOV_PRG_CTX *ctx = PRG_init(seed_sol);
    if (ctx == NULL)
        return -1;
    int rc = RejSampPRG(ctx, 0, n4, TSUOV_t_RS_LEN, r4);
    PRG_final(ctx);
    if (rc != 0)
        return -1;
    memcpy(dst, r4, (size_t)n4);
    return 0;
}

/* =====================================================================
   linear algebra (echelon form)
   ===================================================================== */

TYPEDEF_STRUCT(ROW,
  Fq *col;
  int original_row_id;
);

TYPEDEF_STRUCT(ECHELON_FORM,
  ROW_s row[TSUOV_m2];
  ROW_s *eqn[TSUOV_m2];
  int rank;
  int index[TSUOV_m2];
);

static void echelon_form_init(Fq mat[TSUOV_m2][TSUOV_m2], ECHELON_FORM echelon_form)
{
    for (int i = 0; i < TSUOV_m2; i++) {
        echelon_form->row[i].col = mat[i];
        echelon_form->row[i].original_row_id = i;
        echelon_form->eqn[i] = echelon_form->row + i;
    }
    echelon_form->rank = 0;
    memset(echelon_form->index, -1, sizeof(echelon_form->index));
}

static void row_swap(ROW_s *eqn[TSUOV_m2], int i, int j)
{
    ROW_s *tmp = eqn[i];
    eqn[i] = eqn[j];
    eqn[j] = tmp;
}

#define EQN(I, J) (eqn[I]->col[J])

static void LU_decompose(Fq A[TSUOV_m2][TSUOV_m2], ECHELON_FORM echelon_form)
{
    echelon_form_init(A, echelon_form);
    ROW_s **eqn = echelon_form->eqn;

    int c = -1;
    for (int i = 0; i < TSUOV_m2; i++) {
        c++;
        if (c >= TSUOV_m2)
            return;
        int j = i;
        while (EQN(j, c) == 0) {
            j++;
            if (j >= TSUOV_m2) {
                c++;
                if (c >= TSUOV_m2)
                    return;
                j = i;
            }
        }

        row_swap(eqn, i, j);
        echelon_form->index[echelon_form->rank++] = c;

        Fq pivot = EQN(i, c);
        Fq inv = Fq_inv(pivot);

        EQN(i, i) = pivot;
        for (int k = c + 1; k < TSUOV_m2; k++)
            EQN(i, k) = Fq_mul(inv, EQN(i, k));

        for (j = i + 1; j < TSUOV_m2; j++) {
            Fq mul = EQN(j, c);
            EQN(j, i) = mul;
            for (int k = c + 1; k < TSUOV_m2; k++)
                EQN(j, k) = Fq_sub(EQN(j, k), Fq_mul(mul, EQN(i, k)));
        }
    }
}

static void Fq_mxm_identity(Fq A[TSUOV_m2][TSUOV_m2])
{
    memset(A, 0, sizeof(Fq) * TSUOV_m2 * TSUOV_m2);
    for (int i = 0; i < TSUOV_m2; i++)
        A[i][i] = 1;
}

static void L_inverse(ECHELON_FORM echelon_form, Fq R[TSUOV_m2][TSUOV_m2])
{
    ROW_s **eqn = echelon_form->eqn;
    int rank = echelon_form->rank;

    Fq_mxm_identity(R);

    for (int i = 0; i < rank; i++) {
        Fq inv = Fq_inv(EQN(i, i));
        for (int k = 0; k <= i; k++)
            R[i][k] = Fq_mul(R[i][k], inv);
        for (int j = i + 1; j < TSUOV_m2; j++)
            for (int k = 0; k <= i; k++)
                R[j][k] = Fq_sub(R[j][k], Fq_mul(EQN(j, i), R[i][k]));
    }
}

static int consistent(ECHELON_FORM echelon_form, Fq b[TSUOV_m2], int *cacheR,
                      Fq R[TSUOV_m2][TSUOV_m2])
{
    ROW_s **eqn = echelon_form->eqn;
    int rank = echelon_form->rank;
    if (rank == TSUOV_m2)
        return 1;

    if (*cacheR == 0) {
        L_inverse(echelon_form, R);
        *cacheR = 1;
    }

    for (int i = rank; i < TSUOV_m2; i++) {
        uint64_t t = 0;
        for (int j = 0; j < rank; j++) {
            int k = eqn[j]->original_row_id;
            t += ((uint64_t)R[i][j]) * ((uint64_t)b[k]);
        }
        int k = eqn[i]->original_row_id;
        t += ((uint64_t)R[i][i]) * ((uint64_t)b[k]);
        t %= TSUOV_q;
        if (t)
            return 0;
    }
    return 1;
}

static int sample_a_solution(const TSUOV_SEED seed_sol, ECHELON_FORM echelon_form,
                             Fq b[TSUOV_m2], Fq x[TSUOV_m2], Fq b2[TSUOV_m2])
{
    int rank = echelon_form->rank;
    int *index = echelon_form->index;
    ROW_s **eqn = echelon_form->eqn;

    uint8_t random_buff[TSUOV_m2];
    uint8_t *random_element = random_buff;
    if (rank < TSUOV_m2 && Expand_sol(seed_sol, random_buff) != 0)
        return -1;

    for (int i = 0; i < rank; i++) {
        uint64_t t = 0;
        for (int j = 0; j < i; j++)
            t += ((uint64_t)EQN(i, j)) * ((uint64_t)b2[j]);
        Fq tmp = (Fq)(t % TSUOV_q);
        int k = eqn[i]->original_row_id;
        tmp = Fq_sub(b[k], tmp);
        b2[i] = Fq_mul(tmp, Fq_inv(EQN(i, i)));
    }

    int i = TSUOV_m2 - 1;
    for (int j = rank - 1; j >= 0; j--) {
        int k = index[j];
        for (; i > k; i--)
            x[i] = *random_element++;
        uint64_t t = 0;
        for (k++; k < TSUOV_m2; k++)
            t += ((uint64_t)EQN(j, k)) * ((uint64_t)(x[k]));
        x[i] = Fq_sub(b2[j], (Fq)(t % TSUOV_q));
        i--;
    }
    for (; i >= 0; i--)
        x[i] = *random_element++;
    return 0;
}

static void pack_0(Fq oil_u[TSUOV_m2], MATRIX_kxO oil)
{
    Fq *s = oil_u;
    for (int k = 0; k < TSUOV_k; k++)
        for (int i = 0; i < TSUOV_O; i++) {
            uint8_t b0 = *s++;
            uint8_t b1 = *s++;
            oil[k][1][i] = Fq_mul(b1, 4);
            oil[k][0][i] = Fq_sub(b0, oil[k][1][i]);
        }
}

static void SIG_GEN(const MATRIX_kxO oil, const MATRIX_OxV SdT, const MATRIX_kxV vineger,
                    TSUOV_SIGNATURE sig)
{
    MATRIX_kxV t;
    MATRIX_VxO_dot_VECTOR_O_whipk(SdT, oil, t);
    VECTOR_V_sub_VECTOR_V_whipk_save_in_sign(vineger, t, sig->s);
    for (int k = 0; k < TSUOV_k; k++)
        for (int l = 0; l < TSUOV_L; l++)
            memcpy(&sig->s[k][l][TSUOV_V], &oil[k][l][0], TSUOV_O * sizeof(uint8_t));
}

/* =====================================================================
   KeyGen / Sign / Verify
   ===================================================================== */

int TSUOV_KeyGen(const TSUOV_SEED seed_sk, const TSUOV_SEED seed_pk, TSUOV_P3 P3)
{
    MATRIX_OxN SdT;
    MATRIX_NxN Pi;
    MATRIX_OxN TMP;
    int rc = -1;

    if (Expand_S(seed_sk, SdT) != 0)
        goto out_s;
    TSUOV_PRG_CTX *ctx = PRG_init(seed_pk);
    if (ctx == NULL)
        goto out_s;
    for (int i = 0; i < TSUOV_m1; i++) {
        if (Expand_P_MATRIX(ctx, i, Pi) != 0)
            goto out_ctx;
        gf31_t2_matrix_mul((uint8_t *)Pi, (uint8_t *)SdT, (uint8_t *)TMP, TSUOV_N, TSUOV_N,
                          TSUOV_O, (TSUOV_N + 31) & ~31);
        gf31_t2_matrix_mul((uint8_t *)SdT, (uint8_t *)TMP, (uint8_t *)P3[i], TSUOV_O, TSUOV_N,
                          TSUOV_O, TSUOV_O);
    }
    rc = 0;

out_ctx:
    PRG_final(ctx);
out_s:
    secure_zero(TMP, sizeof(TMP));
    secure_zero(SdT, sizeof(MATRIX_OxV));
    return rc;
}

int TSUOV_Sign(const TSUOV_SEED seed_sk, const TSUOV_SEED seed_pk,
               const TSUOV_SEED seed_v, const TSUOV_SEED seed_r, const TSUOV_SEED seed_sol,
               const uint8_t message[], const size_t message_length, TSUOV_SIGNATURE sig)
{
    ECHELON_FORM echelon_form;
    MATRIX_kxO oil;
    Fq oil_u[TSUOV_m2];
    Fq b[TSUOV_m2];
    Fq b2[TSUOV_m2];
    vector_m1 U[(TSUOV_k * (TSUOV_k + 1)) / 2];
    Fq c[TSUOV_m2];
    Fq M[TSUOV_m1][TSUOV_k][TSUOV_o];
    Fq eqn[TSUOV_m2][TSUOV_k * TSUOV_o];
    Fq R[TSUOV_m2][TSUOV_m2];
    int cacheR = 0;
    Fq msg[TSUOV_m2];

    MATRIX_OxV SdT;
    MATRIX_OxV Fi2T;
    MATRIX_VxV Pi1;
    MATRIX_kxV y;
    VECTOR_V *vT = y;
    MATRIX_kxV vT_Pi1;
    VECTOR_V *vT_Fi1 = vT_Pi1;
    MATRIX_kxV phi_vT;
    Fq vT_Fi2[TSUOV_k][TSUOV_L][TSUOV_O];
    Fq vTk_Pi1_Sd[TSUOV_L][TSUOV_O];

    if (Expand_sk(seed_sk, SdT) != 0)
        return -1;
    if (Expand_V(seed_v, vT) != 0) {
        secure_zero(SdT, sizeof(MATRIX_OxV));
        return -1;
    }
    for (int i = 0; i < TSUOV_k; i++)
        phi_t2_q31_planar((uint8_t *)vT[i], (uint8_t *)phi_vT[i], BLOCK_ALIGNED_BYTES(TSUOV_V));

    TSUOV_PRG_CTX *ctx_pk = PRG_init(seed_pk);
    if (ctx_pk == NULL) {
        secure_zero(SdT, sizeof(MATRIX_OxV));
        return -1;
    }

    for (int i = 0; i < TSUOV_m1; i++) {
        if (Expand_pk(ctx_pk, i, Pi1, Fi2T) != 0) {
            PRG_final(ctx_pk);
            secure_zero(SdT, sizeof(MATRIX_OxV));
            return -1;
        }

        for (int k = 0; k < TSUOV_k; k++) {
            int TSUOV_Vpadded = (TSUOV_V + 31) & ~31;
            gf31_t2_matrix_vec_mul((uint8_t *)vT[k], (uint8_t *)Pi1, (uint8_t *)vT_Pi1[k], TSUOV_V,
                              TSUOV_V, TSUOV_Vpadded);
            gf31_t2_matrix_vec_mul((uint8_t *)(vT_Pi1[k]), (uint8_t *)SdT, (uint8_t *)vTk_Pi1_Sd,
                              TSUOV_V, TSUOV_O, TSUOV_O);
            gf31_t2_matrix_vec_mul((uint8_t *)vT[k], (uint8_t *)Fi2T, (uint8_t *)vT_Fi2[k], TSUOV_V,
                              TSUOV_O, TSUOV_O);
            gf31_vector_add((uint8_t *)vT_Fi2[k], (uint8_t *)vTk_Pi1_Sd, TSUOV_O * TSUOV_L);

            for (int j = 0; j < TSUOV_O; j++) {
                Fq u0 = vT_Fi2[k][0][j];
                Fq u1 = vT_Fi2[k][1][j];
                M[i][k][TSUOV_L * j] = Fq_add(u0, u1);
                M[i][k][TSUOV_L * j + 1] = Fq_mul(8, u1);
            }

            size_t rowk_start_index = (size_t)(k * TSUOV_k - (k * (k + 1)) / 2);
            phi_t2_q31_planar((uint8_t *)vT_Fi1[k], (uint8_t *)vT_Fi1[k],
                              BLOCK_ALIGNED_BYTES(TSUOV_V));
            for (int j = k; j < TSUOV_k; j++) {
                uint64_t cv = vector_v_dot_vector_v(vT_Fi1[k], phi_vT[j]);
                U[rowk_start_index + j][i] = (Fq)(cv % TSUOV_q);
            }
        }
    }

    PRG_final(ctx_pk);

    Fq MT[TSUOV_k * TSUOV_o][TSUOV_m2];
    MATRIX_TRANSPOSE_m2xLO(M, MT);

    Fq eqnT[TSUOV_k * TSUOV_o][TSUOV_m2];
    memset(eqnT, 0, sizeof(eqnT));
    memset(c, 0, sizeof(c));

    init_mul_table_mod31();

    int ctr = 0;
    for (int i = 0; i < TSUOV_k; i++) {
        int rowi_start_index = i * TSUOV_k - (i * (i + 1)) / 2;
        for (int j = TSUOV_k - 1; j >= i; j--) {
            multiply_E_add_m1_m2(U[rowi_start_index + j], c, ctr);
            ctr++;
        }
    }

    int ctr1 = 0, ctr2 = TSUOV_k * (TSUOV_k + 1) / 2 - 3;
    Fq sum[TSUOV_o][TSUOV_m2];
    memset(sum, 0, sizeof(sum));
    for (int j = TSUOV_k - 1; j >= 1; j--) {
        multiply_E_add_mat_m2((Fq *)MT[j * TSUOV_o], (Fq *)sum, ctr1, TSUOV_o);
        multiply_E_add_mat_m2((Fq *)sum, (Fq *)eqnT[(j - 1) * TSUOV_o], ctr2, TSUOV_o);
        ctr1++;
        ctr2 -= TSUOV_k - j + 2;
    }

    ctr1 = 0;
    ctr2 = TSUOV_k - 1;
    memset(sum, 0, sizeof(sum));
    for (int j = 0; j < TSUOV_k; j++) {
        Fq tmp[TSUOV_o][TSUOV_m2];
        multiply_E_mat_m2((Fq *)MT[j * TSUOV_o], (Fq *)tmp, ctr1, TSUOV_o);
        gf31_vector_add((uint8_t *)sum, (uint8_t *)tmp, TSUOV_o * TSUOV_m2);
        gf31_vector_add((uint8_t *)tmp, (uint8_t *)sum, TSUOV_o * TSUOV_m2);
        multiply_E_add_mat_m2((Fq *)tmp, (Fq *)eqnT[j * TSUOV_o], ctr2, TSUOV_o);
        ctr1 += TSUOV_k - j;
        ctr2--;
    }

    for (int i = 0; i < TSUOV_m2; i++)
        for (int j = 0; j < TSUOV_o * TSUOV_k; j++)
            eqn[i][j] = (Fq)(eqnT[j][i]);

    LU_decompose(eqn, echelon_form);

    uint8_t mu[TSUOV_MU_LEN];
    if (Expand_mu(seed_pk, message, message_length, mu) != 0) {
        secure_zero(SdT, sizeof(MATRIX_OxV));
        return -1;
    }

    TSUOV_PRG2_CTX *ctx_r = PRG2_init(seed_r);
    if (ctx_r == NULL) {
        secure_zero(SdT, sizeof(MATRIX_OxV));
        return -1;
    }
    int rc;
    do {
        PRG2_yield(ctx_r, TSUOV_SALT_LEN, sig->r);
        rc = Hash(mu, sig->r, msg);
        if (rc != 0)
            break;
        for (int i = 0; i < TSUOV_m2; i++)
            b[i] = Fq_add(msg[i], c[i]);
    } while (!consistent(echelon_form, b, &cacheR, R));
    PRG2_final(ctx_r);
    if (rc != 0) {
        secure_zero(SdT, sizeof(MATRIX_OxV));
        return -1;
    }

    if (sample_a_solution(seed_sol, echelon_form, b, oil_u, b2) != 0) {
        secure_zero(SdT, sizeof(MATRIX_OxV));
        return -1;
    }
    pack_0(oil_u, oil);
    SIG_GEN(oil, SdT, y, sig);
    secure_zero(SdT, sizeof(MATRIX_OxV));
    return 0;
}

int TSUOV_Verify(const TSUOV_SEED seed_pk, const TSUOV_P3 P3, const uint8_t message[],
                 const size_t message_length, const TSUOV_SIGNATURE sig)
{
    uint8_t mu[TSUOV_MU_LEN];
    if (Expand_mu(seed_pk, message, message_length, mu) != 0)
        return -1;
    Fq msg[TSUOV_m2];
    if (Hash(mu, sig->r, msg) != 0)
        return -1;

    MATRIX_NxN Pi;
    vector_m1 U[(TSUOV_k * (TSUOV_k + 1)) / 2];
    MATRIX_kxN phi_s;

    for (int i = 0; i < TSUOV_k; i++)
        phi_t2_q31_planar((uint8_t *)sig->s[i], (uint8_t *)phi_s[i], BLOCK_ALIGNED_BYTES(TSUOV_N));

    int verify = 1;
    TSUOV_PRG_CTX *ctx_pk = PRG_init(seed_pk);
    if (ctx_pk == NULL)
        return -1;

    for (int i = 0; i < TSUOV_m1; i++) {
        if (Expand_PK(P3[i], ctx_pk, i, Pi) != 0) {
            PRG_final(ctx_pk);
            return -1;
        }
        for (int k1 = 0; k1 < TSUOV_k; k1++) {
            size_t rowk_start_index = (size_t)(k1 * TSUOV_k - (k1 * (k1 + 1)) / 2);
            VECTOR_N sk1_Pi;
            VECTOR_N phi_sk1_Pi;
            gf31_t2_matrix_vec_mul((uint8_t *)sig->s[k1], (uint8_t *)Pi, (uint8_t *)sk1_Pi, TSUOV_N,
                              TSUOV_N, BLOCK_ALIGNED_BYTES(TSUOV_N));
            phi_t2_q31_planar((uint8_t *)sk1_Pi, (uint8_t *)phi_sk1_Pi,
                              BLOCK_ALIGNED_BYTES(TSUOV_N));
            for (int k2 = k1; k2 < TSUOV_k; k2++)
                U[rowk_start_index + k2][i] =
                    gf31_vector_dot((uint8_t *)phi_sk1_Pi, (uint8_t *)phi_s[k2],
                               2 * BLOCK_ALIGNED_BYTES(TSUOV_N));
        }
    }

    Fq acc[TSUOV_m2];
    memset(acc, 0, sizeof(acc));
    init_mul_table_mod31();

    int ctr = 0;
    for (int i = 0; i < TSUOV_k; i++) {
        size_t rowi_start_index = (size_t)(i * TSUOV_k - (i * (i + 1)) / 2);
        for (int j = TSUOV_k - 1; j >= i; j--) {
            multiply_E_add_m1_m2(U[rowi_start_index + j], acc, ctr);
            ctr++;
        }
    }

    for (int i = 0; i < TSUOV_m2; i++) {
        verify &= (acc[i] == msg[i]);
        if (!verify)
            break;
    }
    PRG_final(ctx_pk);
    return verify;
}

void store_TSUOV_P3(const TSUOV_P3 P3, uint8_t *pool, size_t *pool_bits)
{
    for (int i = 0; i < TSUOV_m1; i++)
        for (int j = 0; j < TSUOV_O; j++)
            for (int k = 0; k < TSUOV_O; k++)
                if (k >= j)
                    for (int n = 0; n < TSUOV_L; n++)
                        store_Fq(P3[i][j][n][k], pool, pool_bits);
}

int restore_TSUOV_P3(const uint8_t *pool, size_t *pool_bits, TSUOV_P3 P3)
{
    for (int i = 0; i < TSUOV_m1; i++)
        for (int j = 0; j < TSUOV_O; j++)
            for (int k = 0; k < TSUOV_O; k++)
                if (k < j)
                    for (int n = 0; n < TSUOV_L; n++)
                        P3[i][j][n][k] = P3[i][k][n][j];
                else
                    for (int n = 0; n < TSUOV_L; n++)
                        if (restore_Fq(pool, pool_bits, &P3[i][j][n][k]) != 0)
                            return -1;

    /* The public key ends with 5-bit coefficients followed by zero padding;
     * non-zero padding bits would give a second valid encoding (sign-31-6). */
    size_t used_bits = *pool_bits;
    int rem = (int)(used_bits & 7);
    if (rem != 0) {
        size_t last = (used_bits - 1) >> 3;
        if ((pool[last] >> rem) != 0)
            return -1;
    }
    return 0;
}
