#include "gf/gf31_matrix.h"

#include <string.h>

void gf31_vector_add(gf31_t *c, const gf31_t *d, int len)
{
    for (int i = 0; i < len; ++i)
        c[i] = gf31_add(c[i], d[i]);
}

gf31_t gf31_vector_dot(const gf31_t *a, const gf31_t *b, size_t len)
{
    uint32_t total = 0;
    for (size_t i = 0; i < len; ++i)
        total += (uint32_t)a[i] * b[i];
    return (gf31_t)(total % 31u);
}

void gf31_vector_mul_scalar(const gf31_t *v, gf31_t scalar, size_t len, gf31_t *res)
{
    for (size_t i = 0; i < len; ++i)
        res[i] = gf31_mul(v[i], scalar);
}

void gf31_vector_add_no_mod(const gf31_t *a, const gf31_t *b, size_t len, gf31_t *res)
{
    for (size_t i = 0; i < len; ++i)
        res[i] = (gf31_t)(a[i] + b[i]);
}

void gf31_vector_reduce(const gf31_t *v, size_t len, gf31_t *res)
{
    for (size_t i = 0; i < len; ++i)
        res[i] = gf31_reduce((int)v[i]);
}

void gf31_vector_add_mod(const gf31_t *a, const gf31_t *b, size_t len, gf31_t *res)
{
    for (size_t i = 0; i < len; ++i)
        res[i] = gf31_reduce((int)a[i] + (int)b[i]);
}

void gf31_vector_sub_mod(const gf31_t *a, const gf31_t *b, size_t len, gf31_t *res)
{
    for (size_t i = 0; i < len; ++i)
        res[i] = gf31_sub(a[i], b[i]);
}

void gf31_t2_matrix_mul(const uint8_t *a, const uint8_t *b, uint8_t *c,
                        int a_rows, int a_cols, int b_cols, int c_rows)
{
    const int padded = (a_cols + 31) & ~31;

    for (int i = 0; i < a_rows; ++i) {
        const uint8_t *a0 = &a[i * 2 * padded];
        const uint8_t *a1 = &a[i * 2 * padded + padded];
        for (int j = 0; j < b_cols; ++j) {
            uint32_t acc00 = 0, acc11 = 0, acc01 = 0, acc10 = 0;
            for (int k = 0; k < padded; ++k) {
                acc00 += (uint32_t)a0[k] * b[j * 2 * padded + k];
                acc11 += (uint32_t)a1[k] * b[j * 2 * padded + padded + k];
                acc01 += (uint32_t)a0[k] * b[j * 2 * padded + padded + k];
                acc10 += (uint32_t)a1[k] * b[j * 2 * padded + k];
            }
            c[j * 2 * c_rows + i] = (uint8_t)((acc00 + 3u * acc11) % 31u);
            c[j * 2 * c_rows + c_rows + i] = (uint8_t)((acc01 + acc10) % 31u);
        }
    }
}

void gf31_t2_matrix_vec_mul(const uint8_t *v, const uint8_t *m, uint8_t *c,
                            int v_cols, int m_cols, int c_cols)
{
    const int padded = (v_cols + 31) & ~31;
    const uint8_t *v0 = v;
    const uint8_t *v1 = &v[padded];

    for (int j = 0; j < m_cols; ++j) {
        uint32_t acc00 = 0, acc11 = 0, acc01 = 0, acc10 = 0;
        const uint8_t *m0 = &m[j * 2 * padded];
        const uint8_t *m1 = &m[j * 2 * padded + padded];
        for (int k = 0; k < padded; ++k) {
            acc00 += (uint32_t)v0[k] * m0[k];
            acc11 += (uint32_t)v1[k] * m1[k];
            acc01 += (uint32_t)v0[k] * m1[k];
            acc10 += (uint32_t)v1[k] * m0[k];
        }
        c[j] = (uint8_t)((acc00 + 3u * acc11) % 31u);
        c[c_cols + j] = (uint8_t)((acc01 + acc10) % 31u);
    }
}

void gf31_t2_phi_planar(uint8_t *data, uint8_t *dst, size_t n_pairs)
{
    for (size_t i = 0; i < n_pairs; ++i) {
        const uint8_t a0 = data[i];
        const uint8_t a1 = data[i + n_pairs];
        dst[i] = gf31_add(a0, a1);
        dst[i + n_pairs] = gf31_mul(8, a1);
    }
}
