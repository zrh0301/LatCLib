#ifndef PQC_NGCC_AIGIS_NAMESPACE_H
#define PQC_NGCC_AIGIS_NAMESPACE_H

/* Keep the reference KEM entry points distinct when linking other modules. */
#define mdkem_enc pqc_aigis_mdkem_enc
#define mkem_dec pqc_aigis_mkem_dec
#define mkem_enc pqc_aigis_mkem_enc
#define mkem_keygen pqc_aigis_mkem_keygen

#endif
