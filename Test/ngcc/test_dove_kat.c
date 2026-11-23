/* Reproduces the ICCS KAT_SIG schedule for DOVE and compares it byte for byte
 * with the vector file shipped in the submission package.
 *
 * The schedule is the one implemented by the package's KAT_SIG.c: both the
 * seed and message generators are the ICCS DRNG seeded with fixed nonces, ten
 * records are produced and the message length grows by 8 bytes per record. */

#include "SIG_AlgorithmInstance.h"
#include "ngcc/drng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SEED_LEN_BYTES 64

/* Defined by alg/ngcc_dove/adapter.c for SIG_AlgorithmInstance.c. */
extern DRNG_ctx drng_algorithm;

static void print_len(FILE *f, const char *id, unsigned long long len)
{
    fprintf(f, "%s%llu\n", id, len);
}

static void print_hex(FILE *f, const char *id, const unsigned char *msg,
                      unsigned long long len)
{
    fprintf(f, "%s", id);
    for (unsigned long long i = 0; i < len; i++)
        fprintf(f, "%02X", msg[i]);
    fprintf(f, "\n");
}

static unsigned char *read_file(const char *path, size_t *size)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL)
        return NULL;
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return NULL;
    }
    long n = ftell(f);
    if (n < 0 || fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        return NULL;
    }
    unsigned char *buf = malloc((size_t)n + 1);
    if (buf == NULL) {
        fclose(f);
        return NULL;
    }
    if (fread(buf, 1, (size_t)n, f) != (size_t)n) {
        free(buf);
        fclose(f);
        return NULL;
    }
    fclose(f);
    *size = (size_t)n;
    return buf;
}

static int generate(FILE *out)
{
    unsigned char nonce1[SEED_LEN_BYTES];
    unsigned char nonce2[SEED_LEN_BYTES];
    DRNG_ctx drng_seed;
    DRNG_ctx drng_msg;

    for (int i = 0; i < SEED_LEN_BYTES / 4; i++)
        memcpy(nonce1 + 4 * i, "seed", 4);
    for (int i = 0; i < SEED_LEN_BYTES / 3; i++)
        memcpy(nonce2 + 3 * i, "msg", 3);
    memcpy(nonce2 + SEED_LEN_BYTES - 1, "m", 1);
    if (init_random_number(&drng_seed, nonce1, SEED_LEN_BYTES) != 0 ||
        init_random_number(&drng_msg, nonce2, SEED_LEN_BYTES) != 0)
        return -1;

    unsigned long long pk_len = sig_get_pk_len_bytes();
    unsigned long long sk_len = sig_get_sk_len_bytes();
    unsigned long long sn_len = sig_get_sn_len_bytes();
    unsigned char *pk = calloc((size_t)pk_len, 1);
    unsigned char *sk = calloc((size_t)sk_len, 1);
    unsigned char *sn = calloc((size_t)sn_len, 1);
    unsigned char seed[SEED_LEN_BYTES];
    unsigned char m[128];
    int m_len = 56;
    int rc = 0;

    if (pk == NULL || sk == NULL || sn == NULL) {
        rc = -1;
        goto out;
    }

    for (int i = 0; i < 10; i++) {
        fprintf(out, "Count = %d\n", i);
        if (get_random_number(&drng_seed, seed, SEED_LEN_BYTES * 8) != 0 ||
            get_random_number(&drng_msg, m, (unsigned long long)m_len * 8) != 0) {
            rc = -1;
            goto out;
        }
        fprintf(out, "Seed_Len = %d\n", SEED_LEN_BYTES);
        print_hex(out, "Seed = ", seed, SEED_LEN_BYTES);

        if (init_random_number(&drng_algorithm, seed, SEED_LEN_BYTES) != 0) {
            rc = -1;
            goto out;
        }
        if (sig_keygen(pk, &pk_len, sk, &sk_len) != 0) {
            rc = -1;
            goto out;
        }
        print_len(out, "PK_Len = ", pk_len);
        print_hex(out, "PK = ", pk, pk_len);
        print_len(out, "SK_Len = ", sk_len);
        print_hex(out, "SK = ", sk, sk_len);
        fprintf(out, "M_Len = %d\n", m_len);
        print_hex(out, "M = ", m, (unsigned long long)m_len);

        if (sig_sign(sk, sk_len, m, (unsigned long long)m_len, sn, &sn_len) != 0) {
            rc = -1;
            goto out;
        }
        print_len(out, "Sn_Len = ", sn_len);
        print_hex(out, "Sn = ", sn, sn_len);

        if (sig_verify(pk, pk_len, sn, sn_len, m,
                       (unsigned long long)m_len) != 0) {
            rc = -1;
            goto out;
        }
        fprintf(out, "\n");
        m_len += 8;
    }

out:
    free(pk);
    free(sk);
    free(sn);
    return rc;
}

int main(void)
{
    FILE *produced = tmpfile();
    if (produced == NULL) {
        fprintf(stderr, "cannot create temporary file\n");
        return 1;
    }
    if (generate(produced) != 0) {
        fprintf(stderr, "KAT generation failed\n");
        return 1;
    }

    if (fflush(produced) != 0 || fseek(produced, 0, SEEK_END) != 0) {
        fprintf(stderr, "cannot rewind temporary file\n");
        return 1;
    }
    long produced_len = ftell(produced);
    unsigned char *got = malloc((size_t)produced_len);
    if (got == NULL || fseek(produced, 0, SEEK_SET) != 0 ||
        fread(got, 1, (size_t)produced_len, produced) != (size_t)produced_len) {
        fprintf(stderr, "cannot read back generated KAT\n");
        return 1;
    }

    size_t want_len = 0;
    unsigned char *want = read_file(PQC_DOVE_KAT_FILE, &want_len);
    if (want == NULL) {
        fprintf(stderr, "cannot open %s\n", PQC_DOVE_KAT_FILE);
        return 1;
    }

    if ((size_t)produced_len != want_len ||
        memcmp(got, want, want_len) != 0) {
        size_t limit = (size_t)produced_len < want_len ? (size_t)produced_len : want_len;
        size_t i = 0;
        while (i < limit && got[i] == want[i])
            i++;
        unsigned long line = 1;
        for (size_t k = 0; k < i; k++)
            if (want[k] == '\n')
                line++;
        fprintf(stderr,
                "KAT mismatch at line %lu (offset %zu): generated %zu bytes, "
                "expected %zu bytes\n",
                line, i, (size_t)produced_len, want_len);
        if (i < (size_t)produced_len)
            fprintf(stderr, "generated: %02X\n", got[i]);
        if (i < want_len)
            fprintf(stderr, "expected : %02X\n", want[i]);
        return 1;
    }

    printf("DOVE KAT matches %s (%zu bytes)\n", PQC_DOVE_KAT_FILE, want_len);
    return 0;
}
