#include <stdio.h>
#include <string.h>
#include "keccak.h"

void print_hex(const char *label, const unsigned char *data, size_t len) {
    printf("%s: ", label);
    for (size_t i = 0; i < len; i++) {
        printf("%02x", data[i]);
    }
    printf("\n");
}

int main() {
    // 测试 1: SHAKE128 (空输入)
    // 预期输出 (前16字节): 7f9c2ba4e88f827d616045507605853e
    unsigned char out1[128];
    shake128(out1, sizeof(out1), NULL, 0);
    print_hex("SHAKE128(Empty)", out1, sizeof(out1));

    // 测试 2: SHAKE256 (输入 "abc")
    // 预期输出 (前16字节): 4e03657aea857861940160937a346e29
    unsigned char out2[256];
    const char *msg = "abc";
    shake256(out2, sizeof(out2), (const unsigned char*)msg, strlen(msg));
    print_hex("SHAKE256('abc')", out2, sizeof(out2));

    // 测试 3: 流式接口测试
    KECCAK1600_CTX ctx;
    unsigned char out3[16];
    shake128_init(&ctx);
    sha3_update(&ctx, (const unsigned char*)"a", 1);
    sha3_update(&ctx, (const unsigned char*)"bc", 2);
    shake_xof_final(&ctx);
    shake_squeeze(&ctx, out3, sizeof(out3));
    
    // 应该与测试2的结果(如果测试2也是SHAKE128)一致，或者验证"abc"在SHAKE128下的值
    // SHAKE128("abc") -> 5881092dd818bf5cf8a3ddb793fbcba7
    print_hex("SHAKE128_Stream", out3, sizeof(out3));

    return 0;
}