#include "gf/byte_matrix.h"

#include <string.h>

void gf_byte_add(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t len)
{
    for (size_t i = 0; i < len; ++i)
        out[i] = a[i] ^ b[i];
}

void gf_byte_transpose(uint8_t *out, const uint8_t *in, size_t rows, size_t cols)
{
    for (size_t row = 0; row < rows; ++row)
        for (size_t col = 0; col < cols; ++col)
            out[col * rows + row] = in[row * cols + col];
}

void gf_byte_fold_upper(uint8_t *matrix, size_t n)
{
    for (size_t row = 0; row < n; ++row)
        for (size_t col = row + 1; col < n; ++col) {
            matrix[row * n + col] ^= matrix[col * n + row];
            matrix[col * n + row] = 0;
        }
}

void gf_byte_pack_upper(uint8_t *out, const uint8_t *in, size_t n)
{
    for (size_t row = 0; row < n; ++row) {
        memcpy(out, in + row * n + row, n - row);
        out += n - row;
    }
}

void gf_byte_unpack_upper(uint8_t *out, const uint8_t *in, size_t n)
{
    for (size_t row = 0; row < n; ++row) {
        memset(out, 0, row);
        memcpy(out + row, in, n - row);
        in += n - row;
        out += n;
    }
}

void gf_byte_batch_fold_upper_add(uint8_t *out, const uint8_t *in, size_t n,
                                  size_t batch_bytes)
{
    for (size_t row = 0; row < n; ++row) {
        gf_byte_add(out, out, in + (row * n + row) * batch_bytes,
                    (n - row) * batch_bytes);
        out += batch_bytes;
        for (size_t col = row + 1; col < n; ++col) {
            gf_byte_add(out, out, in + (col * n + row) * batch_bytes, batch_bytes);
            out += batch_bytes;
        }
    }
}
