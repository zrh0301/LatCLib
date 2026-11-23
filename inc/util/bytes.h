#ifndef PQC_BYTES_H
#define PQC_BYTES_H

#include <stddef.h>
#include <stdint.h>

/* Big-endian byte strings. Store zero-extends or truncates to len bytes.
 * Load accepts up to eight bytes. Zero length permits a null pointer. */
static inline void pqc_store_be(uint8_t *out, size_t len, uint64_t value)
{
    while (len > 0) {
        out[--len] = (uint8_t)value;
        value >>= 8;
    }
}

static inline uint64_t pqc_load_be(const uint8_t *in, size_t len)
{
    uint64_t value = 0;
    for (size_t i = 0; i < len; ++i)
        value = (value << 8) | in[i];
    return value;
}

#endif
