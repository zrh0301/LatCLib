#ifndef SM3_H
#define SM3_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#define SM3_HASH_SIZE 32  /*!< bytes of SM3 algorithm hash value */
#define SM3_HASH_BITS 256 /*!< bits of SM3 algorithm hash value */
#define SM3_BLOCK_SIZE 64 /*!< bytes of SM3 algorithm block (512 bits) */

#if defined(USE_NAMESPACE)
#define sm3_init   NAMESPACE(sm3_init)
#define sm3_update NAMESPACE(sm3_update)
#define sm3_final  NAMESPACE(sm3_final)
#define sm3_hash   NAMESPACE(sm3_hash)
#endif

#ifdef __cplusplus
extern "C"
{
#endif
    /**
     * @brief          SM3 context structure
     *
     * Streaming API with the same shape as the shared Keccak/SHAKE layer
     * (keccak.h): init -> update* -> final(md, ctx).
     */
    typedef struct _sm3_context_st
    {
        uint32_t bytes[2];  /*!< number of bytes processed */
        uint32_t state[8];  /*!< intermediate digest state  */
        uint8_t buffer[64]; /*!< data block being processed */
    } sm3_ctx;

    /**
     * @brief          SM3 context setup (reset to the IV)
     *
     * @param  ctx     context to be initialized
     * @return         1 on success, 0 on error
     */
    int sm3_init(sm3_ctx *ctx);

    /**
     * @brief          SM3 process buffer
     *
     * @param  ctx     SM3 context
     * @param  input   buffer holding the input data (may be NULL when ilen == 0)
     * @param  ilen    length of the input data in bytes
     * @return         1 on success, 0 on error
     */
    int sm3_update(sm3_ctx *ctx, const void *input, size_t ilen);

    /**
     * @brief          SM3 final digest: applies the Merkle-Damgard padding and
     *                 writes the 32-byte digest (@see sha3_final)
     *
     * @param  md      SM3 hash result (SM3_HASH_SIZE bytes)
     * @param  ctx     SM3 context
     * @return         1 on success, 0 on error
     */
    int sm3_final(uint8_t md[SM3_HASH_SIZE], sm3_ctx *ctx);

    /**
     * @brief          Output = SM3( input buffer )
     *
     * @param  md      SM3 hash result
     * @param  input   buffer holding the input data
     * @param  ilen    length of the input data in bytes
     */
    void sm3_hash(uint8_t md[SM3_HASH_SIZE], const uint8_t *input, size_t ilen);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* !SM3_H */
