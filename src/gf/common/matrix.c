#include "gf256_matrix.h"
#include <string.h>
static int valid(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t r,size_t k,size_t c){return o&&a&&b&&r&&k&&c&&r<=SIZE_MAX/c&&r<=SIZE_MAX/k&&k<=SIZE_MAX/c;}
int gf256_matrix_mul_add(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t r,size_t k,size_t c){if(!valid(o,a,b,r,k,c))return -1;for(size_t i=0;i<r;i++)for(size_t j=0;j<k;j++)gf256_vector_madd(o+i*c,b+j*c,a[i*k+j],c);return 0;}
int gf256_matrix_mul(uint8_t *o,const uint8_t *a,const uint8_t *b,size_t r,size_t k,size_t c){if(!valid(o,a,b,r,k,c))return -1;memset(o,0,r*c);return gf256_matrix_mul_add(o,a,b,r,k,c);}
