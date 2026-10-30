#ifndef PQC_SECURE_CLEAR_H
#define PQC_SECURE_CLEAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 将 ptr 指向的 len 字节置零; 经 volatile 指针写入, 即使开全优化也不会被当作死代码消除。 */
void pqc_secure_clear(void *ptr, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* PQC_SECURE_CLEAR_H */
