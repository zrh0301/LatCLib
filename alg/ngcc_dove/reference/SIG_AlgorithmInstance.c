/*
The software is provided by the Institute of Commercial Cryptography Standards
(ICCS), and is used for algorithm submissions in the Next-generation Commercial
Cryptographic Algorithms Program (NGCC).

ICCS doesn't represent or warrant that the operation of the software will be
uninterrupted or error-free in all cases. ICCS will take no responsibility for
the use of the software or the results thereof, if the software is used for any
other purposes.
*/

#include "SIG_AlgorithmInstance.h"
#include "ngcc/drng.h"
#include "ngcc/auxfunc.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
// DRNG_ctx for generating pseudorandom numbers within the SIG scheme
extern DRNG_ctx drng_algorithm;




int pseudoXOFwithSeedAndCounter(unsigned long long output_len_bytes, const unsigned char *msg, unsigned long long msg_len_bytes, uint32_t counter,unsigned char *output){
    uint8_t *buf = malloc(msg_len_bytes + 2);
    memcpy(buf, msg, msg_len_bytes);
    buf[msg_len_bytes] = counter & 0xff;
    buf[msg_len_bytes + 1] = (counter >> 8) & 0xff;
    pseudoXOF(output_len_bytes*8, buf, (msg_len_bytes+2)*8, output);
    free(buf);
    return 0;
}

// The following should be used to get pseudorandom numbers
// get_random_number(&drng_algorithm, random_number, random_number_len_bits);

unsigned long long sig_get_pk_len_bytes()
{
	return DOVE_PK_SEED_BYTES + DOVE_M * UPPER_SIZE(DOVE_O);

}

unsigned long long sig_get_sk_len_bytes()
{
	 return DOVE_SK_SEED_BYTES;
}

unsigned long long sig_get_sn_len_bytes()
{
	 return DOVE_N + DOVE_SALT_BYTES;
}

//int pseudoXOF(unsigned long long output_len_bits, const unsigned char *msg, unsigned long long msg_len_bits, unsigned char *output)

// 从 seed_pk 生成第 k 个上三角矩阵 A_k
static void generate_Ak(gf256_matrix_t *Ak, const uint8_t *seed_pk, uint32_t k) {
    size_t out_len = UPPER_SIZE(DOVE_V);
    uint8_t *out = malloc(out_len);
    pseudoXOFwithSeedAndCounter(out_len,seed_pk, DOVE_PK_SEED_BYTES, k, out);
    matrix_fill_upper(Ak, out);
    free(out);
}

static void generate_Pi2(gf256_matrix_t *Pi2, const uint8_t *seed_pk, uint32_t i) {
    size_t out_len = DOVE_V * DOVE_O;
    uint32_t idx = DOVE_L + i;
    uint8_t *out = malloc(out_len);
	pseudoXOFwithSeedAndCounter(out_len,seed_pk, DOVE_PK_SEED_BYTES,idx,out);
    matrix_fill(Pi2, out);
    free(out);
}

static void derive_seedpk_T(uint8_t *seed_pk, gf256_matrix_t *T, const uint8_t *seed_sk) {
    size_t out_len = DOVE_PK_SEED_BYTES + DOVE_V * DOVE_O;
    uint8_t *out = malloc(out_len);
    pseudoXOF(out_len*8, seed_sk, DOVE_SK_SEED_BYTES*8, out);

    memcpy(seed_pk, out, DOVE_PK_SEED_BYTES);
    matrix_fill(T, out + DOVE_PK_SEED_BYTES);
    free(out);
}

int sig_keygen(
	unsigned char *pk, unsigned long long *pk_len_bytes,
	unsigned char *sk, unsigned long long *sk_len_bytes)
{
    // uint8_t seed_sk[DOVE_SK_SEED_BYTES]={0x2c, 0xdc, 0xe3, 0xd6, 0x97, 0xaf, 0x37, 0x16, 0xa9, 0xb3, 0xcd, 0xf0, 0x68, 0xb4, 0x3e, 0x51, 0x38, 0x46, 0xe1, 0x7c, 0xc9, 0xfd, 0x42, 0x79,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
    uint8_t seed_sk[DOVE_SK_SEED_BYTES];
    get_random_number(&drng_algorithm, seed_sk, (unsigned long long)DOVE_SK_SEED_BYTES*8 );
    uint8_t seed_pk[DOVE_PK_SEED_BYTES];
    gf256_matrix_t *T = matrix_create(DOVE_V, DOVE_O);
    if (!T) return -1;
    derive_seedpk_T(seed_pk, T, seed_sk);

    gf256_matrix_t *Ak[DOVE_L];
    gf256_matrix_t *Bk[DOVE_L]; // O×V
    gf256_matrix_t *Ck[DOVE_L]; // V×O
    gf256_matrix_t *Tt = matrix_create(DOVE_O, DOVE_V);
    matrix_transpose(Tt, T);

    for (int k = 0; k < DOVE_L; k++) {
        Ak[k] = matrix_create(DOVE_V, DOVE_V);
        generate_Ak(Ak[k], seed_pk, k);
        Bk[k] = matrix_create(DOVE_O, DOVE_V);
        matrix_mul(Bk[k], Tt, Ak[k]);
        Ck[k] = matrix_create(DOVE_V, DOVE_O);
        matrix_mul(Ck[k], Ak[k], T);
    }

    // 组装公钥
    uint8_t *pk_ptr = pk;
    memcpy(pk_ptr, seed_pk, DOVE_PK_SEED_BYTES);
    pk_ptr += DOVE_PK_SEED_BYTES;

    gf256_matrix_t *Pi2 = matrix_create(DOVE_V, DOVE_O);
    gf256_matrix_t *tmp1 = matrix_create(DOVE_O, DOVE_O);
    gf256_matrix_t *tmp2 = matrix_create(DOVE_O, DOVE_O);
    gf256_matrix_t *Pi3 = matrix_create(DOVE_O, DOVE_O);

    for (int i = 0; i < DOVE_M; i++) {
        generate_Pi2(Pi2, seed_pk, i);
        
        int k1 = i / DOVE_L;
        int k2 = i % DOVE_L;

        matrix_mul(tmp1, Bk[k1], Ck[k2]);
        matrix_mul(tmp2, Tt, Pi2);
        matrix_add(Pi3, tmp1, tmp2);
        matrix_upper(Pi3);
        matrix_dump_upper(Pi3, pk_ptr);
        pk_ptr += UPPER_SIZE(DOVE_O);
    }

    // 组装私钥
    uint8_t *sk_ptr = sk;
    memcpy(sk_ptr, seed_sk, DOVE_SK_SEED_BYTES);
    // sk_ptr += DOVE_SK_SEED_BYTES;
    // matrix_dump(T, sk_ptr);

    // 释放资源
    matrix_free(T);
    matrix_free(Tt);
    matrix_free(Pi2);
    matrix_free(tmp1);
    matrix_free(tmp2);
    matrix_free(Pi3);
    for (int k = 0; k < DOVE_L; k++) {
        matrix_free(Ak[k]);
        matrix_free(Bk[k]);
        matrix_free(Ck[k]);
    }

    *pk_len_bytes = sig_get_pk_len_bytes();
    *sk_len_bytes = sig_get_sk_len_bytes();

    return 0;
}

int sig_sign(
	unsigned char *sk, unsigned long long sk_len_bytes,
	unsigned char *m, unsigned long long m_len_bytes,
	unsigned char *sn, unsigned long long *sn_len_bytes)
{
    /* The reference read DOVE_V*DOVE_O bytes past the end of sk here and then
     * overwrote T immediately in derive_seedpk_T(). Removed: it changed no
     * output, only the memory reads. */
    if (sk == NULL || sn == NULL || sn_len_bytes == NULL ||
        sk_len_bytes < DOVE_SK_SEED_BYTES ||
        *sn_len_bytes < sig_get_sn_len_bytes() ||
        (m == NULL && m_len_bytes != 0))
        return -1;

    const uint8_t *seed_sk = sk;
    gf256_matrix_t *T = matrix_create(DOVE_V, DOVE_O);

    uint8_t seed_pk[DOVE_PK_SEED_BYTES];
    derive_seedpk_T(seed_pk, T, seed_sk);

    // 计算消息哈希 t
    uint8_t t[DOVE_M];
    size_t input_len = m_len_bytes + DOVE_PK_SEED_BYTES;
    uint8_t *input = malloc(input_len);
    memcpy(input, m, m_len_bytes);
    memcpy(input + m_len_bytes, seed_pk, DOVE_PK_SEED_BYTES);
	pseudoXOF(DOVE_M*8,input, input_len*8, t);
    free(input);

    // 计算 salt
    uint8_t salt[DOVE_SALT_BYTES];
    input_len = DOVE_SK_SEED_BYTES + DOVE_M;
    input = malloc(input_len);
    memcpy(input, seed_sk, DOVE_SK_SEED_BYTES);
    memcpy(input + DOVE_SK_SEED_BYTES, t, DOVE_M);
	pseudoXOF(DOVE_SALT_BYTES*8,input, input_len*8, salt);
    free(input);

#ifdef DOVE_BIND_SALT
    /* sign-09-1: recompute the target with the salt bound in. Verification
     * mirrors this, so the salt is authenticated. */
    input_len = m_len_bytes + DOVE_PK_SEED_BYTES + DOVE_SALT_BYTES;
    input = malloc(input_len);
    if (input == NULL)
        return -1;
    memcpy(input, m, m_len_bytes);
    memcpy(input + m_len_bytes, seed_pk, DOVE_PK_SEED_BYTES);
    memcpy(input + m_len_bytes + DOVE_PK_SEED_BYTES, salt, DOVE_SALT_BYTES);
    pseudoXOF(DOVE_M*8, input, input_len*8, t);
    free(input);
#endif

    gf256_matrix_t *Tt = matrix_create(DOVE_O, DOVE_V);

    matrix_transpose(Tt, T);

    int ret = -1;
    uint8_t v[DOVE_V];
    gf256_matrix_t *Ak[DOVE_L];
    uint8_t *xk[DOVE_L], *yk[DOVE_L];
    gf256_matrix_t *Bk[DOVE_L], *Ck[DOVE_L];

    for (int k = 0; k < DOVE_L; k++) {
        Ak[k] = matrix_create(DOVE_V, DOVE_V);
        xk[k] = malloc(DOVE_V);
        yk[k] = malloc(DOVE_V);
        Bk[k] = matrix_create(DOVE_O, DOVE_V);
        Ck[k] = matrix_create(DOVE_V, DOVE_O);
    }

    gf256_matrix_t *L_mat = matrix_create(DOVE_M, DOVE_O);
    uint8_t w[DOVE_M], t_minus_w[DOVE_M], z[DOVE_M], s_v[DOVE_V];
    gf256_matrix_t *Pi2 = matrix_create(DOVE_V, DOVE_O);
    uint8_t term1[DOVE_O], term2[DOVE_O], term3[DOVE_O];

    for (int ctr = 0; ctr < 256; ctr++) {
        // 生成醋变量 v
        input_len = DOVE_M + DOVE_SK_SEED_BYTES;
        input = malloc(input_len);
        memcpy(input, t, DOVE_M);
        memcpy(input + DOVE_M, seed_sk, DOVE_SK_SEED_BYTES);
        pseudoXOFwithSeedAndCounter(DOVE_V,input, input_len,ctr,v);


        free(input);

        for (int k = 0; k < DOVE_L; k++) {

            generate_Ak(Ak[k], seed_pk, k);
            row_vec_mul_matrix(xk[k], v, Ak[k]);
            matrix_mul_col_vec(yk[k], Ak[k], v);
            matrix_mul(Bk[k], Tt, Ak[k]);
            matrix_mul(Ck[k], Ak[k], T);
        }

        // 构造线性系统 L * z = t - w
        memset(L_mat->data, 0, DOVE_M * DOVE_O);
        for (int i = 0; i < DOVE_M; i++) {
            generate_Pi2(Pi2, seed_pk, i);
            int k1 = i / DOVE_L;
            int k2 = i % DOVE_L;

            w[i] = gf256_dot(xk[k1], yk[k2], DOVE_V);

            row_vec_mul_matrix(term1, v, Pi2);
            row_vec_mul_matrix(term2, xk[k1], Ck[k2]);

            uint8_t tmp[DOVE_O];

            matrix_mul_col_vec(tmp, Bk[k1], yk[k2]);
            memcpy(term3, tmp, DOVE_O);
            


            uint8_t *row = &L_mat->data[i * DOVE_O];
            for (int j = 0; j < DOVE_O; j++) {
                row[j] = term1[j] ^ term2[j] ^ term3[j];
            }


        }

        // 计算 t - w（必须在 matrix_solve 之前完成）
        for (int i = 0; i < DOVE_M; i++) {
            t_minus_w[i] = t[i] ^ w[i];
        }

        // 求解线性系统 L * z = t - w
        if (matrix_solve(z, L_mat, t_minus_w) == 0) {

            // 构造签名 s = (v + Tz) || z
            uint8_t Tz[DOVE_V];

            matrix_mul_col_vec(Tz, T, z);
            for (int i = 0; i < DOVE_V; i++) {
                s_v[i] = v[i] ^ Tz[i];
            }

            uint8_t *sig_ptr = sn;
            memcpy(sig_ptr, s_v, DOVE_V);
            sig_ptr += DOVE_V;
            memcpy(sig_ptr, z, DOVE_M);
            sig_ptr += DOVE_M;
            memcpy(sig_ptr, salt, DOVE_SALT_BYTES);

            *sn_len_bytes = sig_get_sn_len_bytes();
            ret = 0;
            break;
        }
    }

    // 释放资源
    matrix_free(T);
    matrix_free(Tt);
    matrix_free(L_mat);
    matrix_free(Pi2);
    for (int k = 0; k < DOVE_L; k++) {

        matrix_free(Ak[k]);
        free(xk[k]);
        free(yk[k]);
        matrix_free(Bk[k]);
        matrix_free(Ck[k]);
    }
    return ret;
}

int sig_verify(
	unsigned char *pk, unsigned long long pk_len_bytes,
	unsigned char *sn, unsigned long long sn_len_bytes,
	unsigned char *m, unsigned long long m_len_bytes)
{
    if (pk == NULL || sn == NULL ||
        pk_len_bytes < sig_get_pk_len_bytes() ||
        sn_len_bytes < sig_get_sn_len_bytes() ||
        (m == NULL && m_len_bytes != 0))
        return -1;

    const uint8_t *seed_pk = pk;
    const uint8_t *pi3_data = pk + DOVE_PK_SEED_BYTES;


    const uint8_t *s = sn;
    const uint8_t *v = s;
    const uint8_t *z = s + DOVE_V;

    // 计算消息哈希 t
    uint8_t t[DOVE_M];
    size_t input_len = m_len_bytes + DOVE_PK_SEED_BYTES;
    uint8_t *input = malloc(input_len);
    memcpy(input, m, m_len_bytes);
    memcpy(input + m_len_bytes, seed_pk, DOVE_PK_SEED_BYTES);
	pseudoXOF(DOVE_M*8,input, input_len*8, t);
    free(input);

#ifdef DOVE_BIND_SALT
    /* sign-09-1: absorb the salt into the target so a different salt yields a
     * different target and verification rejects it. This is a hardening
     * variant, not the submitted construction; it is off by default so the
     * output stays byte-identical to the official KAT. */
    {
        const uint8_t *salt_in = sn + DOVE_N;   /* layout: v || z || salt */
        input_len = m_len_bytes + DOVE_PK_SEED_BYTES + DOVE_SALT_BYTES;
        input = malloc(input_len);
        if (input == NULL)
            return -1;
        memcpy(input, m, m_len_bytes);
        memcpy(input + m_len_bytes, seed_pk, DOVE_PK_SEED_BYTES);
        memcpy(input + m_len_bytes + DOVE_PK_SEED_BYTES, salt_in,
               DOVE_SALT_BYTES);
        pseudoXOF(DOVE_M*8, input, input_len*8, t);
        free(input);
    }
#endif

    gf256_matrix_t *Ak[DOVE_L];
    uint8_t *xk[DOVE_L], *yk[DOVE_L];
    for (int k = 0; k < DOVE_L; k++) {
        Ak[k] = matrix_create(DOVE_V, DOVE_V);
        generate_Ak(Ak[k], seed_pk, k);
        xk[k] = malloc(DOVE_V);
        yk[k] = malloc(DOVE_V);
        row_vec_mul_matrix(xk[k], v, Ak[k]);
        matrix_mul_col_vec(yk[k], Ak[k], v);
    }

    gf256_matrix_t *Pi2 = matrix_create(DOVE_V, DOVE_O);
    gf256_matrix_t *Pi3 = matrix_create(DOVE_O, DOVE_O);
    uint8_t w[DOVE_M];

    for (int i = 0; i < DOVE_M; i++) {
        generate_Pi2(Pi2, seed_pk, i);
        matrix_fill_upper(Pi3, pi3_data + i * UPPER_SIZE(DOVE_O));



        int k1 = i / DOVE_L;
        int k2 = i % DOVE_L;

        uint8_t term1 = gf256_dot(xk[k1], yk[k2], DOVE_V);

        uint8_t pi2_z[DOVE_V];
        matrix_mul_col_vec(pi2_z, Pi2, z);
        uint8_t term2 = gf256_dot(v, pi2_z, DOVE_V);

        uint8_t term3 = quadratic_form(z, Pi3);

        w[i] = term1 ^ term2 ^ term3;
    }

    // 验证二次型结果等于消息哈希 t
    int valid = 0;
    for (int i = 0; i < DOVE_M; i++) {

        if (w[i] != t[i]) {
            valid = -1;
            break;
        }
    }

    matrix_free(Pi2);
    matrix_free(Pi3);
    for (int k = 0; k < DOVE_L; k++) {
        matrix_free(Ak[k]);
        free(xk[k]);
        free(yk[k]);
    }
    return valid;
}


//   void sign_fprintlen(FILE *file_output, char *identifier, unsigned long long len)
// {
// 	if (OUTPUT_BLANK_TEST_VECTORS)
// 		fprintf(file_output, "%s\n", identifier);
// 	else
// 		fprintf(file_output, "%s%llu\n", identifier, len);
// }

//   void sign_fprintstr(FILE *file_output, char *identifier, unsigned char *msg, unsigned long long len)
// {
// 	fprintf(file_output, "%s", identifier);
// 	for (unsigned long long i = 0; i < len; i++)
// 		fprintf(file_output, "%02X", msg[i]);
// 	fprintf(file_output, "\n");
// }