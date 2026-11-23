/*
The software is provided by the Institute of Commercial Cryptography Standards
(ICCS), and is used for algorithm submissions in the Next-generation Commercial
Cryptographic Algorithms Program (NGCC).

ICCS doesn't represent or warrant that the operation of the software will be
uninterrupted or error-free in all cases. ICCS will take no responsibility for
the use of the software or the results thereof, if the software is used for any
other purposes.
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ngcc/drng.h"
#include "sym/sm3.h"
#include "util/bytes.h"

#define OUTLEN (32)
#define DRNG_SUCCESS 0
#define DRNG_MEMORY_ALLOCATION_FAILED -2

#define MIN_INT(a, b) ((a) < (b) ? (a) : (b))
#define MAX_INT(a, b) ((a) > (b) ? (a) : (b))
#define DIVISION_ROUND_UP(dividend, divisor, result) \
    do                                               \
    {                                                \
        (result) = (dividend) / (divisor);           \
        if ((dividend) % (divisor))                  \
            (result)++;                              \
    } while (0)

#define HIGH_N_BIT_MASK(N) ((unsigned char)((~0U) << (8 - (N))))

static void inc_Big_Number(unsigned char *BN, unsigned long long BN_len_bytes)
{
    for (; BN_len_bytes != 0; BN_len_bytes--)
    {
        BN[BN_len_bytes - 1] += 1;
        if (BN[BN_len_bytes - 1])
            break;
    }
}

static void plus_Big_Number(unsigned char *BN1, const unsigned char *BN2, const unsigned char *BN3, const unsigned char *BN4, unsigned long long BN_len_bytes)
{
    unsigned int sum;
    unsigned char carry = 0;

    for (; BN_len_bytes != 0; BN_len_bytes--)
    {
        sum = BN1[BN_len_bytes - 1] + BN2[BN_len_bytes - 1] + BN3[BN_len_bytes - 1] + BN4[BN_len_bytes - 1] + carry;
        carry = sum / (0xFFU + 1);
        BN1[BN_len_bytes - 1] = (unsigned char)(sum & 0xFFU);
    }
}

static int SM3_df(unsigned char *input_string, unsigned long long input_string_len_bytes)
{
    unsigned char temp[(1 + PQC_NGCC_DRNG_SEEDLEN / OUTLEN) * OUTLEN];
    unsigned long long len;
    unsigned char counter;
    unsigned char *data_to_SM3 = NULL;
    unsigned long long data_to_SM3_len_bytes;
    unsigned char sm3_dgst[OUTLEN];
    unsigned char number_of_bits_to_return_big_endian[4];

    memset(temp, 0, sizeof(temp));
    DIVISION_ROUND_UP(PQC_NGCC_DRNG_SEEDLEN, OUTLEN, len);
    counter = 0x01;
    for (unsigned long long i = 0; i < len; i++)
    {
        data_to_SM3_len_bytes = sizeof(counter) + sizeof(number_of_bits_to_return_big_endian) + input_string_len_bytes;
        data_to_SM3 = (unsigned char *)malloc(data_to_SM3_len_bytes);
        if (NULL == data_to_SM3)
        {
            fprintf(stderr, "ERROR: Memory allocation failed at %s, line %d. \n", __FILE__, __LINE__);
            return DRNG_MEMORY_ALLOCATION_FAILED;
        }
        data_to_SM3[0] = counter;
        pqc_store_be(number_of_bits_to_return_big_endian, sizeof(number_of_bits_to_return_big_endian),
                     8 * PQC_NGCC_DRNG_SEEDLEN);
        memcpy(data_to_SM3 + sizeof(counter), number_of_bits_to_return_big_endian, sizeof(number_of_bits_to_return_big_endian));
        memcpy(data_to_SM3 + sizeof(counter) + sizeof(number_of_bits_to_return_big_endian), input_string, input_string_len_bytes);
        sm3_hash(sm3_dgst, data_to_SM3, data_to_SM3_len_bytes);
        memcpy(temp + i * OUTLEN, sm3_dgst, OUTLEN);
        free(data_to_SM3);
        data_to_SM3 = NULL;
        counter++;
    }
    memcpy(input_string, temp, PQC_NGCC_DRNG_SEEDLEN);

    return DRNG_SUCCESS;
}

static int SM3_DRNG_Instantiate(DRNG_ctx *drng, const unsigned char *nonce, unsigned long long nonce_len_bytes)
{
    unsigned char *seed_material = NULL;
    unsigned char seed[PQC_NGCC_DRNG_SEEDLEN];
    unsigned char padded_V[1 + sizeof(drng->V)];

    memset(drng, 0, sizeof(*drng));
    seed_material = (unsigned char *)malloc(MAX_INT(nonce_len_bytes, PQC_NGCC_DRNG_SEEDLEN));
    if (NULL == seed_material)
    {
        fprintf(stderr, "ERROR: Memory allocation failed at %s, line %d. \n", __FILE__, __LINE__);
        return DRNG_MEMORY_ALLOCATION_FAILED;
    }
    memset(seed_material, 0, MAX_INT(nonce_len_bytes, PQC_NGCC_DRNG_SEEDLEN));
    memcpy(seed_material, nonce, nonce_len_bytes);
    SM3_df(seed_material, nonce_len_bytes);
    memcpy(seed, seed_material, sizeof(seed));
    memcpy(drng->V, seed, sizeof(drng->V));
    padded_V[0] = 0x00;
    memcpy(padded_V + 1, drng->V, sizeof(drng->V));
    SM3_df(padded_V, sizeof(padded_V));
    memcpy(drng->C, padded_V, sizeof(drng->C));
    inc_Big_Number(drng->reseed_counter, PQC_NGCC_DRNG_SEEDLEN);

    free(seed_material);
    return DRNG_SUCCESS;
}

static int SM3_DRNG_Generate(DRNG_ctx *drng, unsigned long long requested_no_of_bits, unsigned char *return_bits)
{
    unsigned long long m;
    unsigned char data[PQC_NGCC_DRNG_SEEDLEN];
    unsigned char w[OUTLEN];
    unsigned char padded_V[1 + sizeof(drng->V)];
    unsigned char H[PQC_NGCC_DRNG_SEEDLEN];
    unsigned long long requested_no_of_Byte;
    unsigned long long remainder_Byte;

    DIVISION_ROUND_UP(requested_no_of_bits, OUTLEN * 8, m);
    memcpy(data, drng->V, PQC_NGCC_DRNG_SEEDLEN);
    DIVISION_ROUND_UP(requested_no_of_bits, 8, requested_no_of_Byte);
    remainder_Byte = requested_no_of_Byte;
    for (unsigned long long i = 0; i < m; i++)
    {
        sm3_hash(w, data, sizeof(data));
        if (remainder_Byte >= sizeof(w))
        {
            memcpy(return_bits + OUTLEN * i, w, sizeof(w));
            remainder_Byte -= sizeof(w);
        }
        else
        {
            for (unsigned long long j = 0; j < sizeof(w); j++)
            {
                return_bits[OUTLEN * i + j] = w[j];
                remainder_Byte--;
                if (!remainder_Byte)
                {
                    break;
                }
            }
        }

        inc_Big_Number(data, PQC_NGCC_DRNG_SEEDLEN);
    }

    if (requested_no_of_Byte >= 1)
    {
        return_bits[requested_no_of_Byte - 1] &= HIGH_N_BIT_MASK(8 - (8 * requested_no_of_Byte - requested_no_of_bits));
    }
    memset(H, 0, sizeof(H));
    padded_V[0] = 0x03;
    memcpy(padded_V + 1, drng->V, sizeof(drng->V));
    sm3_hash(H + (PQC_NGCC_DRNG_SEEDLEN - OUTLEN), padded_V, sizeof(padded_V));
    plus_Big_Number(drng->V, H, drng->C, drng->reseed_counter, PQC_NGCC_DRNG_SEEDLEN);
    inc_Big_Number(drng->reseed_counter, PQC_NGCC_DRNG_SEEDLEN);

    return DRNG_SUCCESS;
}

int init_random_number(DRNG_ctx *drng, const unsigned char *seed, unsigned long long seed_len_bytes)
{
    return SM3_DRNG_Instantiate(drng, seed, seed_len_bytes);
}

int get_random_number(DRNG_ctx *drng, unsigned char *random_number, unsigned long long random_number_len_bits)
{
    return SM3_DRNG_Generate(drng, random_number_len_bits, random_number);
}
