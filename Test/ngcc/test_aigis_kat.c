/* Reproduces the ICCS KAT_KEM schedule for Aigis-Enc+ and compares it byte for
 * byte with the vector file shipped in the submission package.
 *
 * The reference draws every random byte from a file-scope DRNG context that
 * the KAT driver seeds before each count. The adapter provides that context as
 * pqc_ngcc_aigis_drng; this test seeds it exactly as the submission driver
 * does and calls the reference entry points directly. */

#include "api.h"
#include "ngcc/drng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SEED_LEN_BYTES 64

/* Defined by the adapter. */
extern DRNG_ctx pqc_ngcc_aigis_drng;

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
    unsigned char nonce[SEED_LEN_BYTES];
    DRNG_ctx drng_seed;
    unsigned char seed[SEED_LEN_BYTES];
    unsigned char *pk;
    unsigned char *sk;
    unsigned char *ct;
    unsigned char *ss;
    unsigned char *ss1;

    for (int i = 0; i < SEED_LEN_BYTES / 4; i++)
        memcpy(nonce + 4 * i, "seed", 4);
    if (init_random_number(&drng_seed, nonce, SEED_LEN_BYTES) != 0)
        return -1;

    pk = calloc(KEM_PUBLICKEYBYTES, 1);
    sk = calloc(KEM_SECRETKEYBYTES, 1);
    ct = calloc(KEM_CIPHERTEXTBYTES, 1);
    ss = calloc(KEM_BYTES, 1);
    ss1 = calloc(KEM_BYTES, 1);
    if (pk == NULL || sk == NULL || ct == NULL || ss == NULL || ss1 == NULL) {
        free(pk);
        free(sk);
        free(ct);
        free(ss);
        free(ss1);
        return -1;
    }

    for (int i = 0; i < 10; i++) {
        fprintf(out, "Count = %d\n", i);
        if (get_random_number(&drng_seed, seed, SEED_LEN_BYTES * 8) != 0) {
            fprintf(stderr, "seed generation failed\n");
            return -1;
        }
        fprintf(out, "Seed_Len = %d\n", SEED_LEN_BYTES);
        print_hex(out, "Seed = ", seed, SEED_LEN_BYTES);

        if (init_random_number(&pqc_ngcc_aigis_drng, seed, SEED_LEN_BYTES) != 0 ||
            mkem_keygen(pk, sk) != 0) {
            fprintf(stderr, "Aigis keygen failed\n");
            return -1;
        }
        print_len(out, "PK_Len = ", KEM_PUBLICKEYBYTES);
        print_hex(out, "PK = ", pk, KEM_PUBLICKEYBYTES);
        print_len(out, "SK_Len = ", KEM_SECRETKEYBYTES);
        print_hex(out, "SK = ", sk, KEM_SECRETKEYBYTES);

        DRNG_ctx encaps_rng = pqc_ngcc_aigis_drng;
        uint8_t coins[KEM_BYTES], ct_direct[KEM_CIPHERTEXTBYTES];
        if (get_random_number(&encaps_rng, coins, sizeof(coins) * 8) != 0 ||
            mkem_enc(pk, ss, ct) != 0 ||
            mdkem_enc(pk, coins, ss1, ct_direct) != 0 ||
            memcmp(ct, ct_direct, sizeof(ct_direct)) != 0 ||
            memcmp(ss, ss1, KEM_BYTES) != 0) {
            fprintf(stderr, "Aigis encapsulation failed\n");
            return -1;
        }
        print_len(out, "CT_Len = ", KEM_CIPHERTEXTBYTES);
        print_hex(out, "CT = ", ct, KEM_CIPHERTEXTBYTES);
        print_len(out, "SS_Len = ", KEM_BYTES);
        print_hex(out, "SS = ", ss, KEM_BYTES);

        if (mkem_dec(sk, ct, ss1) != 0 || memcmp(ss, ss1, KEM_BYTES) != 0) {
            fprintf(stderr, "Aigis decapsulation mismatch\n");
            return -1;
        }
        fprintf(out, "\n");
    }

    free(pk);
    free(sk);
    free(ct);
    free(ss);
    free(ss1);
    return 0;
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
    unsigned char *want = read_file(PQC_AIGIS_KAT_FILE, &want_len);
    if (want == NULL) {
        fprintf(stderr, "cannot open %s\n", PQC_AIGIS_KAT_FILE);
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

    printf("Aigis-Enc+ KAT matches %s (%zu bytes)\n",
           PQC_AIGIS_KAT_FILE, want_len);
    return 0;
}
