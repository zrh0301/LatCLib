#include "core/pqc_conf.h"

#include <stdio.h>

const char *pqc_conf_version(void) {
    return PQC_CONF_VERSION_STRING;
}

void pqc_conf_print(void) {
    printf("pqc %s\n", PQC_CONF_VERSION_STRING);
    printf("  architecture: %s\n", PQC_CONF_ARCH);
    printf("  backend:      %s\n", PQC_CONF_BACKEND);
    printf("  gf backend:   %s\n", PQC_CONF_GF_BACKEND);
    printf("  aes backend:  %s\n", PQC_CONF_AES_BACKEND);
    printf("  libraries:    %s%s\n",
           PQC_CONF_BUILD_SHARED ? "shared " : "",
           PQC_CONF_BUILD_STATIC ? "static" : "");
    printf("  modules:      %s%s%s%s%s%s%s\n",
           PQC_CONF_WITH_DILITHIUM ? "dilithium " : "",
           PQC_CONF_WITH_KYBER ? "kyber " : "",
           PQC_CONF_WITH_SNOVA_R2 ? "snova_r2 " : "",
           PQC_CONF_WITH_SNOVA_R3 ? "snova_r3 " : "",
           PQC_CONF_WITH_MAYO ? "mayo " : "",
           PQC_CONF_WITH_MQOM ? "mqom " : "",
           PQC_CONF_WITH_NGCC_DOVE ? "ngcc_dove" : "");
}
