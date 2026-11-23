#ifndef POLY_H
#define POLY_H

#include "types.h"
#include "ntt.h"

#define MAX_POLY_COEFS 256
typedef lssi_t poly256_i32[MAX_POLY_COEFS];
typedef poly256_i32 poly; /* compatibility spelling, explicitly fixed-size */


#define POLY_ZERO(a, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (a)[_i] = 0; \
} while(0)

#define POLY_COPY(dst, src, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (dst)[_i] = (src)[_i]; \
} while(0)

#define POLY_ADD(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i] + (b)[_i]; \
} while(0)

#define POLY_SUB(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i] - (b)[_i]; \
} while(0)

#define POLY_LSH_LOW(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i]<<(b); \
} while(0)

#define POLY_RSH_LOW(c, a, b, coefs) do { \
    int _i; \
    for(_i=0; _i<(coefs); ++_i) (c)[_i] = (a)[_i]>>(b); \
} while(0)

#define POLY_EQ(a, b, coefs) ({ \
    int _i; \
    int _result = 1; \
    for(_i=0; _i<(coefs); ++_i) { \
        if((a)[_i] != (b)[_i]) { \
            _result = 0; \
            break; \
        } \
    } \
    _result; \
})

#define POLY_IS_ZERO(a, coefs) ({ \
    int _i; \
    int _result = 1; \
    for(_i=0; _i<(coefs); ++_i) { \
        if((a)[_i] != 0) { \
            _result = 0; \
            break; \
        } \
    } \
    _result; \
})

/* Transitional name for the old combined context. New code can use the
 * independent fp32_ctx, ring_ctx and ntt_plan headers directly. */
typedef ntt_plan lp_ctx;

/* Compatibility presets implemented outside the algebra target. */
void ntt_init_dili(lp_ctx *Q);
void ntt_init_kyber(lp_ctx *Q);

#define INIT_NTT_DILI(Q) ntt_init_dili(Q)
#define INIT_NTT_KYBER(Q) ntt_init_kyber(Q)

int poly_mod_equal(const lssi_t *a, const lssi_t *b, const lssi_t q, int qbits, int coefs);
int poly_mod_is_zero(const lssi_t *a, const lssi_t q, int qbits, int coefs);
void poly_neg(lssi_t *c, const lssi_t *a, const lssi_t q, int qbits, int coefs);

void ntt_vector(poly* vec, int count, lp_ctx* Q);
void ntt_matrix(poly* mat, int rows, int cols, lp_ctx* Q);
void poly_mod_barrett_lssi(lssi_t* c, const lssi_t* a, const lssi_t q, int qbits, int coefs);
void poly_mod_barrett_lsdi(lssi_t* c, const lsdi_t* a, const lssi_t q, int qbits, int coefs);
/* n=256 backend, initialized plan required. In-place operations. Forward
 * accepts |a[i]| < q and returns |a[i]| < (layers+1)*q. Inverse accepts
 * arbitrary int32 residues, centers before its butterflies, and returns a
 * lazy Montgomery result in [-q,q]. See ntt.h for the scaling contract. */
void poly_ntt(lssi_t* a, lp_ctx* Q);
void poly_invntt(lssi_t* a, lp_ctx* Q);

/* Multiply compatible transform representations; input magnitudes must be
 * less than (layers+1)*q (the forward transform output bound).
 * Output has one less Montgomery factor and is canonical in [0,q).
 * Supports exact out==a or out==b aliasing; no partial overlap. */
void poly_mul_ntt(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q);

/* Low-level scalar/quadratic kernels require matching block_degree and
 * return lazy residues. The same input bounds and aliasing apply. */
void poly_base_pointwise(lssi_t* r, const lssi_t* a, const lssi_t* b, lp_ctx* Q);
void poly_pointwise(lssi_t* c, const lssi_t* a, const lssi_t* b, lp_ctx* Q);
void poly_mul_small(lssi_t* c, const lssi_t* a, const lssi_t* b, const lssi_t U, lp_ctx* Q);
void poly_mul_small_pub(lssi_t* u, const lssi_t* c, const lssi_t* a, lp_ctx* Q);


#endif // POLY_H
