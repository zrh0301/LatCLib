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
    printf("  modules:      %s%s%s%s%s%s%s%s%s%s\n",
           PQC_CONF_WITH_DILITHIUM ? "dilithium " : "",
           PQC_CONF_WITH_KYBER ? "kyber " : "",
           PQC_CONF_WITH_SNOVA_R2 ? "snova_r2 " : "",
           PQC_CONF_WITH_SNOVA_R3 ? "snova_r3 " : "",
           PQC_CONF_WITH_MAYO ? "mayo " : "",
           PQC_CONF_WITH_MQOM ? "mqom " : "",
           PQC_CONF_WITH_NGCC_DOVE ? "ngcc_dove " : "",
           PQC_CONF_WITH_NGCC_VDOO ? "ngcc_vdoo " : "",
           PQC_CONF_WITH_NGCC_FLEXTREE ? "ngcc_flextree " : "",
           PQC_CONF_WITH_NGCC_AIGIS ? "ngcc_aigis" : "");
    printf("  params:       dilithium=%d kyber=%d mayo=%d "
           "mqom=%d/%d/%d/%d ngcc_dove=%d ngcc_vdoo=%d "
           "ngcc_flextree=%s%s ngcc_aigis=%s\n",
           PQC_CONF_DILITHIUM_MODE, PQC_CONF_KYBER_MODE, PQC_CONF_MAYO_VARIANT,
           PQC_CONF_MQOM_SECURITY, PQC_CONF_MQOM_FIELD, PQC_CONF_MQOM_ROUNDS,
           PQC_CONF_MQOM_TRADEOFF, PQC_CONF_NGCC_DOVE_PARAM,
           PQC_CONF_NGCC_VDOO_PARAM, PQC_CONF_NGCC_FLEXTREE_PARAM,
           PQC_CONF_NGCC_FLEXTREE_STRICT ? "+strict" : "",
           PQC_CONF_NGCC_AIGIS_PARAM);
}
