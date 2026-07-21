#include "test.h"
#include <stdio.h>

int main(void) {
    int ok = 1;

    printf("\n");
    printf("========================================\n");
    printf("  Test Suite: Pack Functions\n");
    printf("========================================\n");

    for(size_t i = 0; i < test_pack_suite_count; i++) {
        const test_pack_suite *suite = test_pack_suites[i];
        printf("\n[Suite] %s\n", suite->name);
        ok &= suite->run();
    }

    printf("\n========================================\n");
    printf("  All Pack Tests Complete\n");
    printf("========================================\n\n");

    return ok ? 0 : 1;
}
