#include "secure_clear.h"

void pqc_secure_clear(void *ptr, size_t len) {
    volatile unsigned char *p = ptr;

    if(ptr == NULL || len == 0) {
        return;
    }
    while(len-- != 0) {
        *p++ = 0;
    }
}
