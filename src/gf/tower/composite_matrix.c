#include "gf/composite_matrix.h"
#include "gf/byte_matrix.h"
#include "gf/gf16_vector.h"

#include <string.h>

static uint8_t element(const uint8_t *v, size_t i, gf_composite_field field)
{
    return field == GF_COMPOSITE_16 ? gf16_packed_get(v, i) : v[i];
}

static void set_element(uint8_t *v, size_t i, uint8_t value, gf_composite_field field)
{
    if (field == GF_COMPOSITE_16)
        gf16_packed_set(v, i, value);
    else
        v[i] = value;
}

static size_t vector_bytes(size_t n, gf_composite_field field)
{
    return field == GF_COMPOSITE_16 ? n / 2 + n % 2 : n;
}

static void (*madd_for(gf_composite_field field))(uint8_t *, const uint8_t *, uint8_t, size_t)
{
    return field == GF_COMPOSITE_16 ? gf16_composite_packed_madd : gf256_composite_madd;
}

static uint8_t (*mul_for(gf_composite_field field))(uint8_t, uint8_t)
{
    return field == GF_COMPOSITE_16 ? gf16_composite_mul : gf256_composite_mul;
}

static uint8_t (*inv_for(gf_composite_field field))(uint8_t)
{
    return field == GF_COMPOSITE_16 ? gf16_composite_inv : gf256_composite_inv;
}

static size_t lower_index(size_t row, size_t col)
{
    return row * (row + 1) / 2 + col;
}

static size_t upper_index(size_t row, size_t col, size_t n)
{
    return row * n - row * (row + 1) / 2 + col;
}

void gf_composite_triangular_inv(uint8_t *out, const uint8_t *in, size_t n,
                                 int upper, gf_composite_field field)
{
    uint8_t (*mul)(uint8_t, uint8_t) = mul_for(field);
    uint8_t (*inv)(uint8_t) = inv_for(field);

    for (size_t i = 0; i < n; ++i) {
        const size_t row = upper ? n - 1 - i : i;
        const size_t d = upper ? upper_index(row, row, n) : lower_index(row, row);
        out[d] = inv(in[d]);
        if (upper) {
            for (size_t col = row + 1; col < n; ++col) {
                uint8_t sum = 0;
                for (size_t k = row + 1; k <= col; ++k)
                    sum ^= mul(in[upper_index(row, k, n)], out[upper_index(k, col, n)]);
                out[upper_index(row, col, n)] = mul(sum, out[d]);
            }
        } else {
            for (size_t col = 0; col < row; ++col) {
                uint8_t sum = 0;
                for (size_t k = col; k < row; ++k)
                    sum ^= mul(in[lower_index(row, k)], out[lower_index(k, col)]);
                out[lower_index(row, col)] = mul(sum, out[d]);
            }
        }
    }
}

void gf_composite_triangular_product(uint8_t *out, const uint8_t *lower,
                                      const uint8_t *upper, size_t n, int upper_first,
                                      gf_composite_field field)
{
    const size_t stride = vector_bytes(n, field);
    uint8_t (*mul)(uint8_t, uint8_t) = mul_for(field);
    memset(out, 0, n * stride);

    for (size_t col = 0; col < n; ++col) {
        for (size_t row = 0; row < n; ++row) {
            uint8_t sum = 0;
            if (upper_first) {
                const size_t start = row > col ? row : col;
                for (size_t k = start; k < n; ++k)
                    sum ^= mul(upper[upper_index(row, k, n)], lower[lower_index(k, col)]);
            } else {
                const size_t end = row < col ? row : col;
                for (size_t k = 0; k <= end; ++k)
                    sum ^= mul(lower[lower_index(row, k)], upper[upper_index(k, col, n)]);
            }
            set_element(out + col * stride, row, sum, field);
        }
    }
}

int gf_composite_matrix_solve(uint8_t *x, const uint8_t *a, const uint8_t *b,
                              size_t n, uint8_t *work, gf_composite_field field)
{
    const size_t stride = n + 1;
    uint8_t (*mul)(uint8_t, uint8_t) = mul_for(field);
    uint8_t (*inv)(uint8_t) = inv_for(field);
    for (size_t row = 0; row < n; ++row) {
        memcpy(work + row * stride, a + row * n, n);
        work[row * stride + n] = b[row];
    }
    for (size_t col = 0; col < n; ++col) {
        uint8_t *pivot = work + col * stride + col;
        if (*pivot == 0) {
            for (size_t row = col + 1; row < n; ++row)
                gf_byte_add(pivot, pivot, work + row * stride + col, stride - col);
            if (*pivot == 0)
                return -1;
        }
        const uint8_t scale = inv(*pivot);
        for (size_t j = 0; j < stride - col; ++j)
            pivot[j] = mul(pivot[j], scale);
        for (size_t row = 0; row < n; ++row) {
            if (row != col) {
                uint8_t *target = work + row * stride + col;
                const uint8_t factor = *target;
                for (size_t j = 0; j < stride - col; ++j)
                    target[j] ^= mul(pivot[j], factor);
            }
        }
    }
    for (size_t row = 0; row < n; ++row)
        x[row] = work[row * stride + n];
    return 0;
}

void gf_composite_col_matvec(uint8_t *out, const uint8_t *a, const uint8_t *x,
                             size_t rows, size_t cols, size_t stride, gf_composite_field field)
{
    void (*madd)(uint8_t *, const uint8_t *, uint8_t, size_t) = madd_for(field);
    const size_t len = vector_bytes(rows, field);
    memset(out, 0, len);
    for (size_t col = 0; col < cols; ++col)
        madd(out, a + col * stride, element(x, col, field), len);
    if (field == GF_COMPOSITE_16 && rows % 2 != 0)
        out[len - 1] &= 15u;
}

void gf_composite_batch_tri_madd(uint8_t *out, const uint8_t *upper,
                                 const uint8_t *scalars, size_t n, size_t scalar_stride,
                                 size_t cols, size_t batch_bytes, gf_composite_field field)
{
    void (*madd)(uint8_t *, const uint8_t *, uint8_t, size_t) = madd_for(field);
    for (size_t row = 0; row < n; ++row) {
        for (size_t col = 0; col < cols; ++col) {
            for (size_t k = row; k < n; ++k)
                madd(out, upper + (k - row) * batch_bytes,
                     element(scalars + col * scalar_stride, k, field), batch_bytes);
            out += batch_bytes;
        }
        upper += (n - row) * batch_bytes;
    }
}

void gf_composite_batch_transpose_madd(uint8_t *out, const uint8_t *scalars,
                                       size_t inner, size_t scalar_stride, size_t rows,
                                       const uint8_t *batch, size_t cols, size_t batch_bytes,
                                       gf_composite_field field)
{
    void (*madd)(uint8_t *, const uint8_t *, uint8_t, size_t) = madd_for(field);
    const size_t row_bytes = cols * batch_bytes;
    for (size_t row = 0; row < rows; ++row) {
        for (size_t k = 0; k < inner; ++k)
            madd(out, batch + k * row_bytes, element(scalars + row * scalar_stride, k, field),
                 row_bytes);
        out += row_bytes;
    }
}

void gf_composite_batch_quad(uint8_t *out, const uint8_t *upper, const uint8_t *x,
                             size_t n, size_t batch_bytes, uint8_t *scratch,
                             gf_composite_field field)
{
    void (*madd)(uint8_t *, const uint8_t *, uint8_t, size_t) = madd_for(field);
    uint8_t (*mul)(uint8_t, uint8_t) = mul_for(field);
    const unsigned order = field == GF_COMPOSITE_16 ? 16 : 256;
    memset(scratch, 0, order * batch_bytes);
    for (size_t row = 0; row < n; ++row) {
        const uint8_t xi = element(x, row, field);
        for (size_t col = row; col < n; ++col) {
            const size_t term = mul(xi, element(x, col, field));
            uint8_t *sum = scratch + term * batch_bytes;
            gf_byte_add(sum, sum, upper, batch_bytes);
            upper += batch_bytes;
        }
    }
    memcpy(out, scratch + batch_bytes, batch_bytes);
    for (unsigned i = 2; i < order; ++i)
        madd(out, scratch + i * batch_bytes, (uint8_t)i, batch_bytes);
}
