#include "gf/gf256_matrix.h"
#include <string.h>
static int valid(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t r,size_t k,size_t c){return o&&a&&b&&r&&k&&c&&r<=SIZE_MAX/c&&r<=SIZE_MAX/k&&k<=SIZE_MAX/c;}
int gf256_matrix_mul_add(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t r,size_t k,size_t c){if(!valid(o,a,b,r,k,c))return -1;for(size_t i=0;i<r;i++)for(size_t j=0;j<k;j++)gf256_vector_madd(o+i*c,b+j*c,a[i*k+j],c);return 0;}
int gf256_matrix_mul(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t r,size_t k,size_t c){if(!valid(o,a,b,r,k,c))return -1;memset(o,0,r*c);return gf256_matrix_mul_add(o,a,b,r,k,c);}

void gf256_matrix_mul_vec(uint8_t *out, const uint8_t *matrix, const uint8_t *v,
                          size_t rows, size_t cols)
{
    for (size_t i = 0; i < rows; ++i)
        out[i] = gf256_vector_dot(matrix + i * cols, v, cols);
}

uint8_t gf256_matrix_quadratic(const uint8_t *v, const uint8_t *matrix, size_t n)
{
    uint8_t result = 0;
    for (size_t i = 0; i < n; ++i)
        result ^= gf256_mul(v[i], gf256_vector_dot(matrix + i * n, v, n));
    return result;
}
