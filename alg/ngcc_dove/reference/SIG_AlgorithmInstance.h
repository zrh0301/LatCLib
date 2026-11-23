/*
The software is provided by the Institute of Commercial Cryptography Standards
(ICCS), and is used for algorithm submissions in the Next-generation Commercial
Cryptographic Algorithms Program (NGCC).

ICCS doesn't represent or warrant that the operation of the software will be
uninterrupted or error-free in all cases. ICCS will take no responsibility for
the use of the software or the results thereof, if the software is used for any
other purposes.
*/

#ifndef SIG_ALGORITHM_INSTANCE_H
#define SIG_ALGORITHM_INSTANCE_H

// Set "OUTPUT_BLANK_TEST_VECTORS" as 0 to generate test vector files
// Set "OUTPUT_BLANK_TEST_VECTORS" as 1 to generate blank template (default)
#define OUTPUT_BLANK_TEST_VECTORS 0

// Set "ALGORITHM_INSTANCE" as your algorithm instance name (no more than 64 bytes)
// Only letters, numbers, '-' or '_' are permitted
#define ALGORITHM_INSTANCE "dove"

#ifndef DOVE_PARAM_SET
#define DOVE_PARAM_SET DOVE512
#endif



#ifdef __cplusplus
extern "C"
{
#endif

#define DOVE128 1
#define DOVE256 2
#define DOVE512 3

	#if DOVE_PARAM_SET == DOVE128
		#define DOVE_O             44
		#define DOVE_M             44
		#define DOVE_N             112
		#define DOVE_V 		      68
		#define DOVE_L             7
		#define DOVE_SALT_BYTES    24
		#define DOVE_SK_SEED_BYTES 24
		#define DOVE_PK_SEED_BYTES 16
	#elif DOVE_PARAM_SET == DOVE256
		#define DOVE_O             96
		#define DOVE_M             96
		#define DOVE_N             256
		#define DOVE_V 			   160
		#define DOVE_L             10
		#define DOVE_SALT_BYTES    40
		#define DOVE_SK_SEED_BYTES 40
		#define DOVE_PK_SEED_BYTES 16
	#elif DOVE_PARAM_SET == DOVE512
		#define DOVE_O             216
		#define DOVE_M             216
		#define DOVE_N             576
		#define DOVE_V 			   360
		#define DOVE_L             15
		#define DOVE_SALT_BYTES    72
		#define DOVE_SK_SEED_BYTES 72
		#define DOVE_PK_SEED_BYTES 16
	#else
		#error "Invalid DOVE_PARAM_SET"
	#endif

#define DOVE_R 8   // 每个域元素比特数
#define DOVE_Q 256 // 域大小

#define UPPER_SIZE(n) ((n) * ((n) + 1) / 2)


	/// @brief Obtain the claimed byte length of the public key
	/// @return Claimed byte length of the public key
	unsigned long long sig_get_pk_len_bytes();

	/// @brief Obtain the claimed byte length of the private key
	/// @return Claimed byte length of the private key
	unsigned long long sig_get_sk_len_bytes();

	/// @brief Obtain the claimed byte length of the signature
	/// @return Claimed byte length of the signature
	unsigned long long sig_get_sn_len_bytes();

	/// @brief Key generate
	/// @param[out] pk Public key
	/// @param[out] pk_len_bytes Byte length of the public key
	/// @param[out] sk Private key
	/// @param[out] sk_len_bytes Byte length of the private key
	/// @return If run successfully, return 0; otherwise, return a self-defined negative (-1 to -99) error code
	int sig_keygen(
		unsigned char *pk, unsigned long long *pk_len_bytes,
		unsigned char *sk, unsigned long long *sk_len_bytes);

	/// @brief Sign
	/// @param[in] sk Private key
	/// @param[in] sk_len_bytes Byte length of the private key
	/// @param[in] m Message
	/// @param[in] m_len_bytes Byte length of the message
	/// @param[out] sn Signature
	/// @param[out] sn_len_bytes Byte length of the signature
	/// @return If run successfully, return 0; otherwise, return a self-defined negative (-1 to -99) error code
	int sig_sign(
		unsigned char *sk, unsigned long long sk_len_bytes,
		unsigned char *m, unsigned long long m_len_bytes,
		unsigned char *sn, unsigned long long *sn_len_bytes);

	/// @brief Verify
	/// @param[in] pk Public key
	/// @param[in] pk_len_bytes Byte length of the public key
	/// @param[in] sn Signature
	/// @param[in] sn_len_bytes Byte length of the signature
	/// @param[in] m Message
	/// @param[in] m_len_bytes Byte length of the message
	/// @return If the signature is valid, return 0; if the signature is invalid, return -1; otherwise, return a self-defined negative (-2 to -99) error code
	int sig_verify(
		unsigned char *pk, unsigned long long pk_len_bytes,
		unsigned char *sn, unsigned long long sn_len_bytes,
		unsigned char *m, unsigned long long m_len_bytes);

#ifdef __cplusplus
}
#endif
#endif