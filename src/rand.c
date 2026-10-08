#include "rand.h"
#include <stdlib.h>

#if defined(_WIN32)
#include <windows.h>
#include <bcrypt.h>
#elif defined(__linux__)
#include <errno.h>
#include <sys/random.h>
#elif defined(__unix__) || defined(__APPLE__)
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#endif

static void clear_output(uint8_t *out, size_t len) {
    volatile uint8_t *p = out;
    while(len-- != 0) {
        *p++ = 0;
    }
}

int pqc_system_random(void *rng_ctx, uint8_t *out, size_t out_len) {
    size_t offset = 0;
    (void)rng_ctx;

    if(out_len == 0) {
        return 0;
    }
    if(out == NULL) {
        return -1;
    }

#if defined(_WIN32)
    while(offset < out_len) {
        ULONG chunk = (ULONG)((out_len - offset > 256) ? 256 : out_len - offset);
        if(BCryptGenRandom(NULL, out + offset, chunk,
                           BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
            clear_output(out, out_len);
            return -1;
        }
        offset += chunk;
    }
#elif defined(__linux__)
    while(offset < out_len) {
        size_t chunk = (out_len - offset > 256) ? 256 : out_len - offset;
        ssize_t received = getrandom(out + offset, chunk, 0);
        if(received < 0 && errno == EINTR) {
            continue;
        }
        if(received <= 0) {
            clear_output(out, out_len);
            return -1;
        }
        offset += (size_t)received;
    }
#elif defined(__unix__) || defined(__APPLE__)
    int fd;
    int flags = O_RDONLY;
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
    do {
        fd = open("/dev/urandom", flags);
    } while(fd < 0 && errno == EINTR);
    if(fd < 0) {
        clear_output(out, out_len);
        return -1;
    }
    while(offset < out_len) {
        size_t chunk = (out_len - offset > 256) ? 256 : out_len - offset;
        ssize_t received = read(fd, out + offset, chunk);
        if(received < 0 && errno == EINTR) {
            continue;
        }
        if(received <= 0) {
            close(fd);
            clear_output(out, out_len);
            return -1;
        }
        offset += (size_t)received;
    }
    close(fd);
#else
    (void)offset;
    clear_output(out, out_len);
    return -1;
#endif
    return 0;
}

int pqc_randombytes(uint8_t *out, size_t out_len, pqc_rng_func rng, void *rng_ctx) {
    if(out_len == 0) {
        return 0;
    }
    if(out == NULL) {
        return -1;
    }
    if(rng == NULL) {
        rng = pqc_system_random;
    }
    if(rng(rng_ctx, out, out_len) != 0) {
        clear_output(out, out_len);
        return -1;
    }
    return 0;
}

void rand_bytes(uint8_t *buf, int size) {
    if(size < 0 || pqc_randombytes(buf, (size_t)size, NULL, NULL) != 0) {
        abort();
    }
}

void rand_bits(uint8_t *buf, int bits) {
    size_t bytes;
    if(bits < 0) {
        abort();
    }
    bytes = (size_t)bits / 8 + (bits % 8 != 0);
    if(pqc_randombytes(buf, bytes, NULL, NULL) != 0) {
        abort();
    }
    if(bits % 8 != 0) {
        buf[bytes - 1] &= (uint8_t)((1u << (bits % 8)) - 1u);
    }
}
