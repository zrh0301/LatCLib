#include "gf/gf256_matrix.h"

#include <stdlib.h>
#include <string.h>

/* Row-major augmented matrix. A zero pivot row is conditionally added to
 * until its pivot is nonzero. All row visits depend only on dimensions. */
static int row_reduce(uint8_t *matrix, size_t n, size_t stride, int inverse)
{
    unsigned nonsingular = 1;
    for (size_t col = 0; col < n; ++col) {
        uint8_t *pivot = matrix + col * stride + col;
        const size_t width = stride - col;
        for (size_t row = col + 1; row < n; ++row) {
            const uint8_t *candidate = matrix + row * stride + col;
            const uint8_t mask = (uint8_t)(0u - (unsigned)(pivot[0] == 0));
            for (size_t j = 0; j < width; ++j)
                pivot[j] ^= candidate[j] & mask;
        }
        nonsingular &= pivot[0] != 0;
        const uint8_t scale = gf256_inv(pivot[0]);
        for (size_t j = 0; j < width; ++j)
            pivot[j] = gf256_mul(pivot[j], scale);
        for (size_t row = inverse ? 0 : col + 1; row < n; ++row) {
            if (row != col) {
                uint8_t *target = matrix + row * stride + col;
                gf256_vector_madd(target, pivot, target[0], width);
            }
        }
    }
    return nonsingular ? 0 : -1;
}

int gf256_matrix_solve(uint8_t *out, const uint8_t *matrix, const uint8_t *rhs,
                       size_t n)
{
    if (!out || !matrix || !rhs || n == 0 || n == SIZE_MAX || n > SIZE_MAX / (n + 1))
        return -1;
    const size_t stride = n + 1;
    uint8_t *work = malloc(n * stride);
    if (!work)
        return -1;

    for (size_t row = 0; row < n; ++row) {
        memcpy(work + row * stride, matrix + row * n, n);
        work[row * stride + n] = rhs[row];
    }

    const int result = row_reduce(work, n, stride, 0);
    if (result == 0) {
        for (size_t row = n; row-- > 0;) {
            const uint8_t *coeffs = work + row * stride;
            out[row] = coeffs[n] ^ gf256_vector_dot(coeffs + row + 1,
                                                    out + row + 1, n - row - 1);
        }
    }
    free(work);
    return result;
}

int gf256_matrix_inv(uint8_t *out, const uint8_t *matrix, size_t n)
{
    if (!out || !matrix || n == 0 || n > SIZE_MAX / 2 || n > SIZE_MAX / (2 * n))
        return -1;
    const size_t stride = 2 * n;
    uint8_t *work = calloc(n, stride);
    if (!work)
        return -1;

    for (size_t row = 0; row < n; ++row) {
        memcpy(work + row * stride, matrix + row * n, n);
        work[row * stride + n + row] = 1;
    }
    const int result = row_reduce(work, n, stride, 1);
    if (result == 0) {
        for (size_t row = 0; row < n; ++row)
            memcpy(out + row * n, work + row * stride + n, n);
    }
    free(work);
    return result;
}
