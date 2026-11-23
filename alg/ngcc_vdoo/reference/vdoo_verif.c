#include "ngcc/auxfunc.h"
#include "vdoo_verif.h"
#include "util/ct.h"

void vdoo_optimized_evaluation(unsigned char *y, const unsigned char *pk, const unsigned char *x)
{
    uint8_t scratch[VDOO_Q * VDOO_M_BYTE];
    gf_composite_batch_quad(y, pk, x, VDOO_N, VDOO_M_BYTE, scratch, VDOO_FIELD);
}


int vdoo_verify(const uint8_t *digest, const uint8_t *signature, const pk_t *pk)
{
    unsigned char computed[VDOO_M_BYTE];
    unsigned char correct[VDOO_M_BYTE];
    unsigned char digest_salt[HASH_LEN + SALT_BYTES];

    const uint8_t *salt = signature + VDOO_N_BYTE;
    memcpy(digest_salt, digest, HASH_LEN);
    memcpy(digest_salt + HASH_LEN, salt, SALT_BYTES);

    if (pseudoXOF(VDOO_M_BYTE * 8, digest_salt, sizeof(digest_salt) * 8, correct) != 0)
        return -1;

    vdoo_optimized_evaluation(computed, pk->pk, signature);

    return pqc_ct_equal(computed, correct, VDOO_M_BYTE) ? 0 : -1;
}
