#ifndef PQC_NGCC_SCLOUD_NAMESPACE_H
#define PQC_NGCC_SCLOUD_NAMESPACE_H

/* Scloud+ 的参考实现使用通用的内部名字（pack_pk、sample_s、shake128、
 * randombytes 等），与其他模块或库自身的 Keccak/AES 层会发生符号冲突。
 * 这里统一加 pqc_ngcc_scloud_ 前缀，保证多个算法可以同时链接。 */

/* KEM / PKE 控制流 */
#define scloud_kemkeygen pqc_ngcc_scloud_kemkeygen
#define scloud_kemencaps pqc_ngcc_scloud_kemencaps
#define scloud_kemdecaps pqc_ngcc_scloud_kemdecaps
#define pke_keygen pqc_ngcc_scloud_pke_keygen
#define pke_enc pqc_ngcc_scloud_pke_enc
#define pke_dec pqc_ngcc_scloud_pke_dec

/* 采样与编码 */
#define sample_s pqc_ngcc_scloud_sample_s
#define sample_sp pqc_ngcc_scloud_sample_sp
#define sample_e pqc_ngcc_scloud_sample_e
#define sample_e12 pqc_ngcc_scloud_sample_e12
#define msg_encode pqc_ngcc_scloud_msg_encode
#define msg_decode pqc_ngcc_scloud_msg_decode

/* 序列化 */
#define pack_sk pqc_ngcc_scloud_pack_sk
#define unpack_sk pqc_ngcc_scloud_unpack_sk
#define pack_pk pqc_ngcc_scloud_pack_pk
#define unpack_pk pqc_ngcc_scloud_unpack_pk
#define reduce_c1 pqc_ngcc_scloud_reduce_c1
#define reduce_c2 pqc_ngcc_scloud_reduce_c2
#define pack_c1 pqc_ngcc_scloud_pack_c1
#define unpack_c1 pqc_ngcc_scloud_unpack_c1
#define pack_c2 pqc_ngcc_scloud_pack_c2
#define unpack_c2 pqc_ngcc_scloud_unpack_c2

/* 矩阵运算 */
#define mat_add pqc_ngcc_scloud_mat_add
#define mat_sub pqc_ngcc_scloud_mat_sub
#define mul_as_e pqc_ngcc_scloud_mul_as_e
#define mul_sa_e pqc_ngcc_scloud_mul_sa_e
#define mul_sb_e pqc_ngcc_scloud_mul_sb_e
#define mul_cs pqc_ngcc_scloud_mul_cs

/* 工具 */
#define scloudplus_verify pqc_ngcc_scloud_verify
#define scloudplus_cmov pqc_ngcc_scloud_cmov
#define scloudplus_secure_zeroize pqc_ngcc_scloud_secure_zeroize

/* 家族选择的 F/G/H/K 原语加前缀避免与其它候选冲突。
 * SHAKE 本身直接使用库内 sym/keccak.c，因此不再重命名 shake128/shake256。 */
#define scloudplus_F pqc_ngcc_scloud_hash_F
#define scloudplus_G pqc_ngcc_scloud_hash_G
#define scloudplus_H pqc_ngcc_scloud_hash_H
#define scloudplus_K pqc_ngcc_scloud_hash_K

/* 参考实现的随机入口由 ngcc_scloud/random.c 提供 */
#define randombytes pqc_ngcc_scloud_randombytes

#endif
