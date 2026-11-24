/* 按 ICCS 的 KAT_KEM 调度复现 Scloud+-128-AES 测试向量，并与提交包中
 * 的向量文件逐字节比对。
 *
 * 参考实现的所有随机字节都来自文件级 DRNG；本库由 adapter 提供
 * pqc_ngcc_scloud_drng，并且 randombytes 直接读取该上下文。这里像官方驱动
 * 一样先播种，再依次调用密钥生成、封装、解封装。 */

#include "kem.h"
#include "ngcc/drng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SEED_LEN_BYTES 64

/* 由 adapter 定义。 */
extern DRNG_ctx pqc_ngcc_scloud_drng;

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

    pk = calloc(SCLOUDPLUS_PUBLIC_KEY_BYTES, 1);
    sk = calloc(SCLOUDPLUS_SECRET_KEY_BYTES, 1);
    ct = calloc(SCLOUDPLUS_CIPHERTEXT_BYTES, 1);
    ss = calloc(SCLOUDPLUS_SHARED_SECRET_BYTES, 1);
    ss1 = calloc(SCLOUDPLUS_SHARED_SECRET_BYTES, 1);
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

        if (init_random_number(&pqc_ngcc_scloud_drng, seed, SEED_LEN_BYTES) != 0 ||
            scloud_kemkeygen(pk, sk) != 0) {
            fprintf(stderr, "Scloud+ keygen failed\n");
            return -1;
        }
        print_len(out, "PK_Len = ", SCLOUDPLUS_PUBLIC_KEY_BYTES);
        print_hex(out, "PK = ", pk, SCLOUDPLUS_PUBLIC_KEY_BYTES);
        print_len(out, "SK_Len = ", SCLOUDPLUS_SECRET_KEY_BYTES);
        print_hex(out, "SK = ", sk, SCLOUDPLUS_SECRET_KEY_BYTES);

        if (scloud_kemencaps(pk, ct, ss) != 0) {
            fprintf(stderr, "Scloud+ encapsulation failed\n");
            return -1;
        }
        print_len(out, "CT_Len = ", SCLOUDPLUS_CIPHERTEXT_BYTES);
        print_hex(out, "CT = ", ct, SCLOUDPLUS_CIPHERTEXT_BYTES);
        print_len(out, "SS_Len = ", SCLOUDPLUS_SHARED_SECRET_BYTES);
        print_hex(out, "SS = ", ss, SCLOUDPLUS_SHARED_SECRET_BYTES);

        if (scloud_kemdecaps(sk, ct, ss1) != 0 ||
            memcmp(ss, ss1, SCLOUDPLUS_SHARED_SECRET_BYTES) != 0) {
            fprintf(stderr, "Scloud+ decapsulation mismatch\n");
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
        fclose(produced);
        return 1;
    }

    if (fflush(produced) != 0 || fseek(produced, 0, SEEK_END) != 0) {
        fprintf(stderr, "cannot rewind temporary file\n");
        fclose(produced);
        return 1;
    }
    long produced_len = ftell(produced);
    unsigned char *got = malloc((size_t)produced_len);
    if (got == NULL || fseek(produced, 0, SEEK_SET) != 0 ||
        fread(got, 1, (size_t)produced_len, produced) != (size_t)produced_len) {
        fprintf(stderr, "cannot read back generated KAT\n");
        free(got);
        fclose(produced);
        return 1;
    }

    size_t want_len = 0;
    unsigned char *want = read_file(PQC_SCLOUD_KAT_FILE, &want_len);
    if (want == NULL) {
        fprintf(stderr, "cannot open %s\n", PQC_SCLOUD_KAT_FILE);
        free(got);
        fclose(produced);
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
        free(want);
        free(got);
        fclose(produced);
        return 1;
    }

    printf("Scloud+ KAT matches %s (%zu bytes)\n",
           PQC_SCLOUD_KAT_FILE, want_len);
    free(want);
    free(got);
    fclose(produced);
    return 0;
}
