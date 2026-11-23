/// @file blas_matrix_ref.h
/// @brief linear algebra functions for matrix op.
///
#ifndef _BLAS_MATRIX_REF_H_
#define _BLAS_MATRIX_REF_H_

#include <stdint.h>



#ifdef  __cplusplus
extern  "C" {
#endif


///////////////// Section: multiplications  ////////////////////////////////

// matrix-vector

void gf16mat_prod_ref(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

void gf256mat_prod_ref(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

/////////////////////////////////////////////////////

unsigned gf256mat_gaussian_elim_ref(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len);

void gf256mat_back_substitute_ref( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len);

unsigned gf16mat_gaussian_elim_ref(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len);

void gf16mat_back_substitute_ref( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len);



#ifdef  __cplusplus
}
#endif

#endif  // _BLAS_MATRIX_REF_H_

