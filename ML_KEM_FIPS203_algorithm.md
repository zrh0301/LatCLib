# ML-KEM FIPS 203 Algorithm Notes

Source: `C:/Users/34621/Downloads/NIST.FIPS.203_ML-KEM.pdf`, FIPS 203, Module-Lattice-Based Key-Encapsulation Mechanism Standard.

This file is an implementation-oriented summary of ML-KEM. It is intended to guide the next coding step in this repository.

## 1. Constants and Parameter Sets

Global constants:

- `N = 256`
- `Q = 3329`
- `ZETA = 17`, a primitive 256-th root of unity modulo `Q`
- Shared secret length: `32` bytes

Approved parameter sets:

| Name | k | eta1 | eta2 | du | dv | ek bytes | dk bytes | ct bytes | RBG strength |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ML-KEM-512 | 2 | 3 | 2 | 10 | 4 | 800 | 1632 | 768 | 128 |
| ML-KEM-768 | 3 | 2 | 2 | 10 | 4 | 1184 | 2400 | 1088 | 192 |
| ML-KEM-1024 | 4 | 2 | 2 | 11 | 5 | 1568 | 3168 | 1568 | 256 |

Derived lengths:

- `ekPKE_len = 384*k + 32`
- `dkPKE_len = 384*k`
- `ct_len = 32 * (du*k + dv)`
- `dkKEM_len = 768*k + 96`
- `poly_bytes_d(d) = 32*d`

## 2. Core Data Model

ML-KEM works over:

- `Rq = Z_q[X] / (X^256 + 1)`
- NTT domain `Tq`, represented as 256 coefficients grouped into 128 base-case degree-one residues.

Implementation types:

- `poly`: 256 coefficients modulo `Q`
- `polyvec`: vector of `k` polynomials
- `matrix`: `k x k` polynomials, usually already in NTT domain

Important representation rule:

- K-PKE public vector `t_hat`, secret vector `s_hat`, and matrix `A_hat` are stored/used in NTT representation.
- Ciphertext component `u` and message polynomial `mu` are in normal coefficient representation.

## 3. Hash and XOF Wrappers

All inputs and outputs are byte arrays.

- `H(s) = SHA3-256(s)`, output 32 bytes
- `J(s) = SHAKE256(s, 32 bytes)`, output 32 bytes
- `G(s) = SHA3-512(s)`, output 64 bytes, parsed as `(a, b)` where each is 32 bytes
- `PRF_eta(s, b) = SHAKE256(s || b, 64*eta bytes)`, where `s` is 32 bytes and `b` is one byte
- `XOF = SHAKE128` incremental API, used by `SampleNTT`

Implementation notes:

- FIPS 203 uses output lengths in bytes in its wrapper functions.
- `G(d || byte(k))` includes the module dimension `k` for keygen domain separation.
- ML-KEM differs from older Kyber behavior: encapsulation does not hash `m` before using it, and the shared secret derivation does not include `H(c)`.

## 4. Bits, Bytes, Encoding, and Compression

### Bits and Bytes

Bits are little-endian inside each byte:

- byte value `B[i] = sum_{j=0..7} bit[8*i+j] * 2^j`

### ByteEncode_d(F)

Input: 256 integers.

- For `1 <= d <= 11`, entries are modulo `2^d`.
- For `d = 12`, entries are modulo `Q`.
- Output length is `32*d` bytes.
- For each coefficient, output its `d` low bits little-endian, then pack bits into bytes.

### ByteDecode_d(B)

Input length: `32*d` bytes.

- Convert bytes to little-endian bits.
- For each coefficient `i`, read `d` bits as `sum bit[i*d+j] * 2^j`.
- If `d < 12`, modulus is `2^d`.
- If `d == 12`, reduce the decoded 12-bit value modulo `Q`.

Important for encapsulation key check:

- Check `ByteEncode_12(ByteDecode_12(ek[0 : 384*k])) == ek[0 : 384*k]`.
- This rejects public-key coefficients that relied on the modulo reduction in `ByteDecode_12`.

### Compress_d and Decompress_d

For `d < 12`:

- `Compress_d(x) = round((2^d / Q) * x) mod 2^d`
- `Decompress_d(y) = round((Q / 2^d) * y)`

Do not use floating point. Suitable integer formulas:

```c
compress_d(x, d)   = (((x << d) + Q / 2) / Q) & ((1u << d) - 1)
decompress_d(y, d) = ((y * Q) + (1u << (d - 1))) >> d
```

Before compression, first normalize `x` into `[0, Q-1]`.

## 5. Sampling

### SampleNTT(B)

Input: 34 bytes `B = rho || byte(j) || byte(i)`.

Output: a polynomial in NTT representation.

Algorithm:

1. Initialize SHAKE128 XOF and absorb `B`.
2. Set output index `j = 0`.
3. While `j < 256`:
   - Squeeze 3 bytes `C[0..2]`.
   - `d1 = C[0] + 256 * (C[1] & 0x0f)`.
   - `d2 = (C[1] >> 4) + 16 * C[2]`.
   - If `d1 < Q`, set `a[j++] = d1`.
   - If `d2 < Q` and `j < 256`, set `a[j++] = d2`.
4. Return `a`.

Loop bound:

- Prefer no hard bound.
- If a bound is required, FIPS Appendix B says it must not be lower than 280 iterations.
- On bound failure, destroy intermediate data and return a uniform error.

### SamplePolyCBD_eta(B)

Input: `64*eta` bytes from `PRF_eta`, where `eta` is 2 or 3.

Output: a polynomial in normal representation.

For each coefficient `i`:

- Let `b = BytesToBits(B)`.
- `x = sum_{j=0..eta-1} b[2*i*eta + j]`
- `y = sum_{j=0..eta-1} b[2*i*eta + eta + j]`
- `f[i] = x - y mod Q`

Implementation note:

- Store small negative values either centered or as `Q - value`; be consistent with NTT and encode paths.
- For serialization, coefficients must be normalized to `[0, Q-1]`.

## 6. NTT and Multiplication

ML-KEM uses a 7-layer NTT:

- Forward loop lengths: `128, 64, 32, 16, 8, 4, 2`
- Inverse loop lengths: `2, 4, 8, 16, 32, 64, 128`
- Inverse final multiplier: `3303`, because `3303 == 128^-1 mod Q`

Forward NTT:

```text
i = 1
for len = 128; len >= 2; len /= 2:
  for start = 0; start < 256; start += 2*len:
    zeta = zetas[i] = ZETA^BitRev7(i) mod Q
    i += 1
    for j = start; j < start + len; j++:
      t = zeta * f[j + len]
      f[j + len] = f[j] - t
      f[j] = f[j] + t
return f
```

Inverse NTT:

```text
i = 127
for len = 2; len <= 128; len *= 2:
  for start = 0; start < 256; start += 2*len:
    zeta = zetas[i] = ZETA^BitRev7(i) mod Q
    i -= 1
    for j = start; j < start + len; j++:
      t = f[j]
      f[j] = t + f[j + len]
      f[j + len] = zeta * (f[j + len] - t)
f = f * 3303 mod Q
return f
```

Multiply two NTT-domain polynomials:

For `i = 0..127`, multiply coefficient pairs:

```text
gamma = ZETA^(2*BitRev7(i) + 1) mod Q
(h[2*i], h[2*i+1]) =
  BaseCaseMultiply(f[2*i], f[2*i+1], g[2*i], g[2*i+1], gamma)
```

Base-case multiplication:

```text
c0 = a0*b0 + a1*b1*gamma mod Q
c1 = a0*b1 + a1*b0 mod Q
```

Current repository note:

- This repository uses one shared 32-bit NTT implementation for Dilithium and Kyber through `lp_ctx`.
- Kyber/ML-KEM must use `INIT_NTT_KYBER(&ctx)` and the existing `poly_ntt`, `poly_invntt`, and Kyber polynomial multiplication path.
- Do not introduce a separate classic 16-bit Kyber NTT path.
- Kyber polynomial multiplication is treated as already validated in this repository; the ML-KEM implementation should focus on upper-layer correctness: sampling, byte encoding, compression, matrix index order, transpose use, and normal-domain versus NTT-domain transitions.
- If using Montgomery arithmetic, keep all zetas and pointwise products in the existing repository convention.
- The standard table in Appendix A is natural modulo `Q`; current dynamic table generation and Montgomery representation should remain centralized in `Poly/src/ntt_table.c`.

## 7. K-PKE Algorithms

K-PKE is the deterministic public-key encryption component used inside ML-KEM.

### K-PKE.KeyGen(d)

Input: 32-byte seed `d`.

Output:

- `ekPKE`, length `384*k + 32`
- `dkPKE`, length `384*k`

Algorithm:

1. `(rho, sigma) = G(d || byte(k))`.
2. Set nonce `N = 0`.
3. Generate matrix `A_hat[k][k]`:
   - For `i = 0..k-1`
   - For `j = 0..k-1`
   - `A_hat[i][j] = SampleNTT(rho || byte(j) || byte(i))`
4. Generate secret vector `s[k]`:
   - `s[i] = SamplePolyCBD_eta1(PRF_eta1(sigma, N))`
   - increment `N`
5. Generate error vector `e[k]`:
   - `e[i] = SamplePolyCBD_eta1(PRF_eta1(sigma, N))`
   - increment `N`
6. `s_hat = NTT(s)`, applied to every vector element.
7. `e_hat = NTT(e)`, applied to every vector element.
8. `t_hat = A_hat o s_hat + e_hat`, where matrix-vector products use `MultiplyNTTs`.
9. `ekPKE = ByteEncode_12(t_hat) || rho`.
10. `dkPKE = ByteEncode_12(s_hat)`.

Implementation notes:

- The matrix index order is exactly `rho || j || i`, not `rho || i || j`.
- `A_hat` is already in NTT form because `SampleNTT` outputs NTT-domain polynomials.
- `s_hat` is the secret key material stored by K-PKE.

### K-PKE.Encrypt(ekPKE, m, r)

Input:

- `ekPKE`, length `384*k + 32`
- plaintext `m`, 32 bytes
- randomness `r`, 32 bytes

Output:

- ciphertext `c`, length `32*(du*k + dv)`

Algorithm:

1. Set nonce `N = 0`.
2. Decode `t_hat = ByteDecode_12(ekPKE[0 : 384*k])`, producing `k` NTT-domain polynomials.
3. Extract `rho = ekPKE[384*k : 384*k + 32]`.
4. Regenerate `A_hat[i][j] = SampleNTT(rho || byte(j) || byte(i))`.
5. Generate `y[k]`:
   - `y[i] = SamplePolyCBD_eta1(PRF_eta1(r, N))`
   - increment `N`
6. Generate `e1[k]`:
   - `e1[i] = SamplePolyCBD_eta2(PRF_eta2(r, N))`
   - increment `N`
7. Generate `e2 = SamplePolyCBD_eta2(PRF_eta2(r, N))`.
8. `y_hat = NTT(y)`.
9. `u = NTT^-1(A_hat^T o y_hat) + e1`.
10. `mu = Decompress_1(ByteDecode_1(m))`.
11. `v = NTT^-1(t_hat^T o y_hat) + e2 + mu`.
12. `c1 = ByteEncode_du(Compress_du(u))`.
13. `c2 = ByteEncode_dv(Compress_dv(v))`.
14. Return `c = c1 || c2`.

Implementation notes:

- The transpose matters. If `A_hat` is row-major, then:
  - `u[i] = INTT(sum_j A_hat[j][i] * y_hat[j]) + e1[i]`
- `v` is a single polynomial:
  - `v = INTT(sum_i t_hat[i] * y_hat[i]) + e2 + mu`
- Compression is coefficient-wise.

### K-PKE.Decrypt(dkPKE, c)

Input:

- `dkPKE`, length `384*k`
- ciphertext `c`, length `32*(du*k + dv)`

Output:

- plaintext `m`, 32 bytes

Algorithm:

1. Split ciphertext:
   - `c1 = c[0 : 32*du*k]`
   - `c2 = c[32*du*k : 32*(du*k + dv)]`
2. `u_prime = Decompress_du(ByteDecode_du(c1))`, producing `k` normal-domain polynomials.
3. `v_prime = Decompress_dv(ByteDecode_dv(c2))`, producing one normal-domain polynomial.
4. Decode `s_hat = ByteDecode_12(dkPKE)`.
5. Compute `u_hat = NTT(u_prime)`.
6. `w = v_prime - NTT^-1(s_hat^T o u_hat)`.
7. `m = ByteEncode_1(Compress_1(w))`.
8. Return `m`.

## 8. ML-KEM Internal Algorithms

The internal algorithms are deterministic and are useful for tests/KATs. They should not be exposed as normal application APIs except for validation use.

### ML-KEM.KeyGen_internal(d, z)

Input:

- `d`: 32-byte seed
- `z`: 32-byte implicit-rejection seed

Algorithm:

1. `(ekPKE, dkPKE) = K-PKE.KeyGen(d)`.
2. `ek = ekPKE`.
3. `dk = dkPKE || ek || H(ek) || z`.
4. Return `(ek, dk)`.

### ML-KEM.Encaps_internal(ek, m)

Input:

- `ek`, length `384*k + 32`
- random message `m`, 32 bytes

Algorithm:

1. `(K, r) = G(m || H(ek))`.
2. `c = K-PKE.Encrypt(ek, m, r)`.
3. Return `(K, c)`.

Implementation note:

- Many C APIs return `(c, K)` instead of `(K, c)`. Be explicit in function names and tests.

### ML-KEM.Decaps_internal(dk, c)

Input:

- `dk`, length `768*k + 96`
- `c`, length `32*(du*k + dv)`

Algorithm:

1. Parse `dkPKE = dk[0 : 384*k]`.
2. Parse `ekPKE = dk[384*k : 768*k + 32]`.
3. Parse `h = dk[768*k + 32 : 768*k + 64]`.
4. Parse `z = dk[768*k + 64 : 768*k + 96]`.
5. `m_prime = K-PKE.Decrypt(dkPKE, c)`.
6. `(K_prime, r_prime) = G(m_prime || h)`.
7. `K_bar = J(z || c)`.
8. `c_prime = K-PKE.Encrypt(ekPKE, m_prime, r_prime)`.
9. If `c != c_prime`, set `K_prime = K_bar`.
10. Return `K_prime`.

Security notes:

- Compare `c` and `c_prime` in constant time.
- Select between `K_prime` and `K_bar` with constant-time conditional move.
- Destroy the implicit-rejection flag and sensitive intermediate data before returning.

## 9. Public ML-KEM APIs

### ML-KEM.KeyGen()

1. Generate random 32-byte `d`.
2. Generate random 32-byte `z`.
3. If random generation fails, return an error.
4. Return `ML-KEM.KeyGen_internal(d, z)`.

### ML-KEM.Encaps(ek)

Required input checks:

1. `ek` must be byte array of length `384*k + 32`.
2. Modulus check:
   - `test = ByteEncode_12(ByteDecode_12(ek[0 : 384*k]))`
   - reject unless `test == ek[0 : 384*k]`

Then:

1. Generate random 32-byte `m`.
2. If random generation fails, return an error.
3. Return `ML-KEM.Encaps_internal(ek, m)`.

### ML-KEM.Decaps(dk, c)

Required input checks:

1. `c` must be byte array of length `32*(du*k + dv)`.
2. `dk` must be byte array of length `768*k + 96`.
3. Hash check:
   - `test = H(dk[384*k : 768*k + 32])`
   - reject unless `test == dk[768*k + 32 : 768*k + 64]`

Then:

1. Return `ML-KEM.Decaps_internal(dk, c)`.

Implementation note:

- Ciphertext length checking must happen on every decapsulation.
- Decapsulation key checking can be assured elsewhere, but the implementation should provide the check path.

## 10. Implementation Plan for This Repository

Recommended file layout:

- `MLKEM/inc/mlkem.h`: parameter struct, public API declarations, sizes.
- `MLKEM/src/mlkem.c`: public and internal KEM algorithms.
- `MLKEM/src/kem_pke.c`: K-PKE KeyGen/Encrypt/Decrypt.
- `MLKEM/src/kem_encode.c`: ByteEncode/ByteDecode, Compress/Decompress.
- `MLKEM/src/kem_sample.c`: SampleNTT, SamplePolyCBD, PRF wrappers.
- Reuse `Poly/src/poly_ntt.c` and dynamic zeta generation for Kyber/ML-KEM mode.

Suggested build order:

1. Add `mlkem_param_t` with `k, eta1, eta2, du, dv, ek_len, dk_len, ct_len`.
2. Implement and test `ByteEncode/ByteDecode` for `d = 1, 4, 5, 10, 11, 12`.
3. Implement and test `Compress/Decompress` with no floating point.
4. Implement `SampleNTT` and `SamplePolyCBD_eta`.
5. Implement base multiplication and matrix/vector NTT-domain products.
6. Implement K-PKE KeyGen/Encrypt/Decrypt.
7. Implement internal ML-KEM algorithms.
8. Implement public input-checked ML-KEM APIs.
9. Add tests:
   - encode/decode round trips for all required `d`
   - self-test `KeyGen -> Encaps -> Decaps`
   - wrong ciphertext decapsulation returns a deterministic 32-byte rejection key, not an error
   - encapsulation-key modulus check rejects invalid 12-bit encodings
   - decapsulation hash check rejects inconsistent `dk`
   - CAVP/KAT vectors when available

## 11. Common Pitfalls

- Do not use old Kyber FO transform details in FIPS 203 ML-KEM.
- Do not hash `m` before `G(m || H(ek))`.
- Do not include `H(c)` in shared secret derivation.
- Do not swap matrix indices: FIPS uses `SampleNTT(rho || j || i)`.
- Do not forget transpose in encryption: `A_hat^T o y_hat`.
- Do not use floating point for compression/decompression.
- Do not expose internal deterministic APIs as normal application APIs.
- Do not return an error on decapsulation ciphertext mismatch; use implicit rejection.
- Keep NTT-domain and normal-domain polynomials separate in names.
- If using Montgomery arithmetic, convert constants and serialization boundaries consistently.

## 12. 与当前 Dilithium 工程风格对齐的设计

本节是下一步实现时的直接蓝图。目标不是先写外层公开 API，而是先把 ML-KEM/Kyber 的内部正确性跑通，并让代码组织、数据结构和基础函数命名尽量贴近当前 `Dilithium` 模块。

当前 Dilithium 风格要点：

- 头文件中用 `DILI_*` 宏定义最大尺寸和参数。
- 用 `dili_param` 在运行时选择模式。
- 用 `poly` 作为统一多项式类型，矩阵和向量用 `poly *`/`poly name[MAX]`。
- 每个主算法内部创建 `lp_ctx ctx; INIT_NTT_DILI(&ctx);`。
- 多项式向量操作集中在 `polyvec_*` 风格函数中。
- 采样函数集中在 `dili_sample.*`。
- pack/unpack 函数集中在 `pack.*`，但 Kyber 模块不建议再叫 `pack.h`，否则会和 `Dilithium/inc/pack.h` 在 include path 上冲突。
- 测试风格是 `init_xxx_param()`、`get_xxx_mode_name()`、`test_xxx_mode()`、`test_xxx()`、`test_perf()`。

建议新增模块目录：

```text
Kyber/inc/kyber.h
Kyber/inc/kyber_pack.h
Kyber/inc/kyber_sample.h
Kyber/inc/kyber_polyvec.h

Kyber/src/kyber.c
Kyber/src/kyber_pack.c
Kyber/src/kyber_sample.c
Kyber/src/kyber_polyvec.c

Test/src/test_kyber.c
```

说明：

- 文件夹名用 `Kyber`，但宏和函数命名建议统一用 `kyber_` 前缀；FIPS 203 名称 `ML-KEM` 可以在注释中说明。
- 不使用裸名 `pack.h`、`polyvec.h`，避免与 Dilithium 头文件重名。
- 外层 `kyber_keygen/kyber_encaps/kyber_decaps` 暂不设计。先实现 deterministic/internal 层与 K-PKE 层。
- Kyber 与 Dilithium 共用 `Poly` 模块中的 NTT 框架；Kyber 只通过 `INIT_NTT_KYBER(&ctx)` 选择参数。
- 不按经典 Kyber 参考实现的 16-bit 运算方式重写底层；本工程 Kyber 使用 32-bit `poly`/`lssi_t` 路径。
- 当前 Kyber 多项式乘法视为已确定正确，后续设计和实现只需要保证上层调用时的域转换、矩阵转置、pack/unpack、compress/decompress 和采样 nonce 顺序正确。

## 13. Kyber 数据结构和参数设计

### 13.1 宏定义

放在 `Kyber/inc/kyber.h`：

```c
#ifndef KYBER_H
#define KYBER_H

#include <stdint.h>
#include "poly.h"
#include "keccak.h"

#define KYBER_N 256
#define KYBER_Q 3329
#define KYBER_QBITS 12

#define KYBER_SYMBYTES 32
#define KYBER_MAX_K 4

#define KYBER_POLYBYTES 384              // ByteEncode_12(poly)
#define KYBER_MAX_EK_BYTES 1568
#define KYBER_MAX_DK_PKE_BYTES 1536
#define KYBER_MAX_CT_BYTES 1568

#endif
```

注意：

- `KYBER_MAX_K = 4` 对应 ML-KEM-1024。
- 运行时 `param.k` 只取 `2, 3, 4`。
- Kyber/ML-KEM 的 `q = 3329` 已由 `INIT_NTT_KYBER(&ctx)` 支持。

### 13.2 参数结构

风格对齐 `dili_param`：

```c
typedef struct {
    int mode;       // KYBER_MODE_2, KYBER_MODE_3, KYBER_MODE_5
    int k;
    int eta1;
    int eta2;
    int du;
    int dv;
    int ek_bytes;
    int dk_pke_bytes;
    int ct_bytes;
} kyber_param;
```

测试中先提供初始化函数，后续可以移动到 `kyber.c`：

```c
void kyber_init_param(kyber_param *param, int mode) {
    param->mode = mode;
    switch(mode) {
        case KYBER_MODE_2:
            param->k = 2;
            param->eta1 = 3;
            param->eta2 = 2;
            param->du = 10;
            param->dv = 4;
            param->ek_bytes = 800;
            param->dk_pke_bytes = 768;
            param->ct_bytes = 768;
            break;
        case KYBER_MODE_3:
            param->k = 3;
            param->eta1 = 2;
            param->eta2 = 2;
            param->du = 10;
            param->dv = 4;
            param->ek_bytes = 1184;
            param->dk_pke_bytes = 1152;
            param->ct_bytes = 1088;
            break;
        case KYBER_MODE_5:
            param->k = 4;
            param->eta1 = 2;
            param->eta2 = 2;
            param->du = 11;
            param->dv = 5;
            param->ek_bytes = 1568;
            param->dk_pke_bytes = 1536;
            param->ct_bytes = 1568;
            break;
        default:
            break;
    }
}
```

### 13.3 内部解包结构

FIPS 的真实接口是 byte-array，但当前 Dilithium 更偏向 unpacked struct。为了正确性和测试方便，建议同时保留 byte-array internal 函数和 unpacked 辅助结构。

```c
typedef struct {
    uint8_t rho[KYBER_SYMBYTES];
    poly t[KYBER_MAX_K];       // t_hat, NTT domain
} kyber_pke_pk;

typedef struct {
    poly s[KYBER_MAX_K];       // s_hat, NTT domain
} kyber_pke_sk;
```

建议正确性主路径仍以 FIPS byte layout 为准：

```c
int kyber_pke_keygen(uint8_t *ek_pke, uint8_t *dk_pke,
                     const uint8_t d[KYBER_SYMBYTES],
                     kyber_param param);

int kyber_pke_encrypt(uint8_t *ct, const uint8_t *ek_pke,
                      const uint8_t m[KYBER_SYMBYTES],
                      const uint8_t r[KYBER_SYMBYTES],
                      kyber_param param);

int kyber_pke_decrypt(uint8_t m[KYBER_SYMBYTES],
                      const uint8_t *dk_pke,
                      const uint8_t *ct,
                      kyber_param param);
```

Current implementation scope: only the bottom K-PKE layer is implemented and
tested. KEM/internal/public Encaps/Decaps APIs are intentionally left out.

这里虽然函数名是 `kyber_*`，语义按 FIPS 203 ML-KEM。

## 14. 基础函数命名设计

### 14.1 哈希和 PRF

放在 `kyber_sample.c` 或 `kyber.c` 的 static helper 中，命名对齐用途：

```c
void kyber_g(uint8_t out[64], const uint8_t *in, size_t inlen);
void kyber_prf(uint8_t *out, size_t outlen,
               const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);
```

伪代码：

```text
kyber_g(out, in, inlen):
    SHA3-512(in, inlen) -> out[64]

kyber_prf(out, outlen, seed, nonce):
    buf = seed || nonce
    SHAKE256(buf, 33, outlen bytes) -> out
```

### 14.2 pack/encode/compress

放在 `Kyber/inc/kyber_pack.h` 和 `Kyber/src/kyber_pack.c`。

```c
void kyber_byte_encode(uint8_t *out, const poly a, int d);
void kyber_byte_decode(poly a, const uint8_t *in, int d);

void kyber_poly_compress(poly r, const poly a, int d);
void kyber_poly_decompress(poly r, const poly a, int d);

void kyber_pack_ek(uint8_t *ek, const poly *t_hat,
                   const uint8_t rho[32], kyber_param param);
void kyber_unpack_ek(poly *t_hat, uint8_t rho[32],
                     const uint8_t *ek, kyber_param param);

void kyber_pack_dk_pke(uint8_t *dk_pke, const poly *s_hat, kyber_param param);
void kyber_unpack_dk_pke(poly *s_hat, const uint8_t *dk_pke, kyber_param param);

void kyber_pack_ct(uint8_t *ct, const poly *u, const poly v, kyber_param param);
void kyber_unpack_ct(poly *u, poly v, const uint8_t *ct, kyber_param param);
```

伪代码：

```text
kyber_byte_encode(out, a, d):
    bitpos = 0
    zero out[0 .. 32*d-1]
    for i in 0..255:
        x = normalize input coefficient
        if d < 12: x = x mod 2^d
        else: x = x mod KYBER_Q
        for j in 0..d-1:
            if (x >> j) & 1:
                out[bitpos >> 3] |= 1 << (bitpos & 7)
            bitpos++

kyber_byte_decode(a, in, d):
    bitpos = 0
    modulus = (d == 12) ? KYBER_Q : (1 << d)
    for i in 0..255:
        x = 0
        for j in 0..d-1:
            bit = (in[bitpos >> 3] >> (bitpos & 7)) & 1
            x |= bit << j
            bitpos++
        a[i] = x % modulus
```

Compression uses integer-only rounding:

```text
kyber_compress_one(x, d):
    x = reduce x into [0, KYBER_Q-1]
    return (((x << d) + KYBER_Q/2) / KYBER_Q) & ((1 << d) - 1)

kyber_decompress_one(y, d):
    return (y * KYBER_Q + (1 << (d-1))) >> d
```

Ciphertext pack:

```text
kyber_pack_ct(ct, u, v, param):
    ptr = ct
    for i in 0..param.k-1:
        tmp = Compress_du(u[i])
        ByteEncode_du(ptr, tmp)
        ptr += 32 * param.du
    tmpv = Compress_dv(v)
    ByteEncode_dv(ptr, tmpv)
```

### 14.3 sampling

放在 `kyber_sample.h/c`，不要复用 `dili_sample.c` 中的 `sample_eta`。FIPS 203 CBD 的 `eta=3` 与 Dili 的 `eta=4` 路径不同，单独写更清楚。

```c
void kyber_sample_ntt(poly a, const uint8_t seed[34], lp_ctx *ctx);
void kyber_sample_cbd(poly a, const uint8_t *buf, int eta);
void kyber_sample_noise(poly a, const uint8_t sigma[32], uint8_t nonce, int eta);
```

伪代码：

```text
kyber_sample_ntt(a, seed, ctx):
    init SHAKE128 ctx
    absorb seed[34]
    ctr = 0
    while ctr < 256:
        squeeze 3 bytes C
        d1 = C[0] + 256 * (C[1] & 0x0f)
        d2 = (C[1] >> 4) + 16 * C[2]
        if d1 < ctx->q: a[ctr++] = d1
        if d2 < ctx->q and ctr < 256: a[ctr++] = d2

kyber_sample_cbd(a, buf, eta):
    bits = BytesToBits(buf)
    for i in 0..255:
        x = sum bits[2*i*eta + 0 .. 2*i*eta + eta-1]
        y = sum bits[2*i*eta + eta .. 2*i*eta + 2*eta-1]
        a[i] = x - y
        if a[i] < 0: a[i] += KYBER_Q

kyber_sample_noise(a, sigma, nonce, eta):
    buf_len = 64 * eta
    buf = PRF_eta(sigma, nonce)
    kyber_sample_cbd(a, buf, eta)
```

### 14.4 polyvec/matrix

放在 `kyber_polyvec.h/c`。不要复用 Dilithium 的 `polyvec_matrix_expand()`，因为它会先采样再 `poly_ntt()`；ML-KEM 的 `SampleNTT` 已经直接输出 NTT 表示。

```c
void kyber_polyvec_matrix_expand(poly *mat, int k,
                                 const uint8_t rho[32], lp_ctx *ctx);
void kyber_polyvec_ntt(poly *v, int len, lp_ctx *ctx);
void kyber_polyvec_invntt(poly *v, int len, lp_ctx *ctx);
void kyber_polyvec_add(poly *w, const poly *u, const poly *v, int len);
void kyber_polyvec_sub(poly *w, const poly *u, const poly *v, int len);
void kyber_polyvec_reduce(poly *v, int len);
void kyber_polyvec_copy(poly *dst, const poly *src, int len);

void kyber_polyvec_pointwise_acc(poly r, const poly *a_hat,
                                 const poly *b_hat, int len, lp_ctx *ctx);
void kyber_matrix_vector_mul(poly *r, const poly *mat,
                             const poly *v_hat, int k, lp_ctx *ctx);
void kyber_matrix_transpose_vector_mul(poly *r, const poly *mat,
                                       const poly *v_hat, int k, lp_ctx *ctx);
```

伪代码：

```text
kyber_polyvec_matrix_expand(mat, k, rho, ctx):
    for i in 0..k-1:
        for j in 0..k-1:
            seed = rho || byte(j) || byte(i)
            kyber_sample_ntt(mat[i*k + j], seed, ctx)

kyber_polyvec_pointwise_acc(r, a_hat, b_hat, len, ctx):
    zero r
    for i in 0..len-1:
        tmp = existing Kyber NTT-domain polynomial multiply(a_hat[i], b_hat[i], ctx)
        r += tmp
    reduce r

kyber_matrix_vector_mul(r, mat, v_hat, k, ctx):
    for i in 0..k-1:
        kyber_polyvec_pointwise_acc(r[i], &mat[i*k], v_hat, k, ctx)

kyber_matrix_transpose_vector_mul(r, mat, v_hat, k, ctx):
    for i in 0..k-1:
        zero acc
        for j in 0..k-1:
            tmp = MultiplyNTTs(mat[j*k + i], v_hat[j])
            acc += tmp
        r[i] = acc
```

重要：

- ML-KEM 上层不要直接调用普通 `poly_pointwise()` 来做 Kyber 多项式乘法；应调用当前工程已经验证过的 Kyber NTT-domain 乘法路径。
- 如果底层函数名仍是 `poly_base_pointwise()`，Kyber 上层可以用 `kyber_poly_basemul()` 或 `kyber_polyvec_pointwise_acc()` 做一层薄封装，封装内部只调用现有正确实现。
- 不新增 16-bit Kyber 乘法和表；所有系数类型仍使用当前 `poly` / `lssi_t` 的 32-bit 路径。
- K-PKE KeyGen 的 `t_hat = A_hat * s_hat + e_hat` 要使用 `kyber_matrix_vector_mul()`。
- K-PKE Encrypt 的 `u = INTT(A_hat^T * y_hat) + e1` 要使用 transpose 版本。
- 上层重点检查：`SampleNTT` 是否已经在 NTT 域、`A_hat` 是否被错误地再次 NTT、`u/v` 是否在 pack 前回到普通系数域。

## 15. 内部算法伪代码（按工程命名）

### 15.1 `kyber_pke_keygen`

```text
kyber_pke_keygen(ek_pke, dk_pke, d, param):
    ctx = INIT_NTT_KYBER()
    k = param.k

    seedbuf = d || byte(k)
    buf64 = G(seedbuf)
    rho = buf64[0..31]
    sigma = buf64[32..63]

    kyber_polyvec_matrix_expand(A_hat, k, rho, &ctx)

    N = 0
    for i in 0..k-1:
        kyber_sample_noise(s[i], sigma, N, param.eta1)
        N++
    for i in 0..k-1:
        kyber_sample_noise(e[i], sigma, N, param.eta1)
        N++

    kyber_polyvec_copy(s_hat, s, k)
    kyber_polyvec_ntt(s_hat, k, &ctx)

    kyber_polyvec_copy(e_hat, e, k)
    kyber_polyvec_ntt(e_hat, k, &ctx)

    kyber_matrix_vector_mul(t_hat, A_hat, s_hat, k, &ctx)
    kyber_polyvec_add(t_hat, t_hat, e_hat, k)
    kyber_polyvec_reduce(t_hat, k)

    kyber_pack_ek(ek_pke, t_hat, rho, param)
    kyber_pack_dk_pke(dk_pke, s_hat, param)
    return 0
```

Correctness checks:

- `ek_pke` length is `384*k + 32`.
- `dk_pke` length is `384*k`.
- `ek_pke[384*k .. 384*k+31] == rho`.
- Re-decoding `ek_pke` gives the same `t_hat` modulo `Q`.

### 15.2 `kyber_pke_encrypt`

```text
kyber_pke_encrypt(ct, ek_pke, m, r, param):
    ctx = INIT_NTT_KYBER()
    k = param.k

    kyber_unpack_ek(t_hat, rho, ek_pke, param)
    kyber_polyvec_matrix_expand(A_hat, k, rho, &ctx)

    N = 0
    for i in 0..k-1:
        kyber_sample_noise(y[i], r, N, param.eta1)
        N++
    for i in 0..k-1:
        kyber_sample_noise(e1[i], r, N, param.eta2)
        N++
    kyber_sample_noise(e2, r, N, param.eta2)

    kyber_polyvec_copy(y_hat, y, k)
    kyber_polyvec_ntt(y_hat, k, &ctx)

    kyber_matrix_transpose_vector_mul(u, A_hat, y_hat, k, &ctx)
    kyber_polyvec_invntt(u, k, &ctx)
    kyber_polyvec_add(u, u, e1, k)
    kyber_polyvec_reduce(u, k)

    kyber_polyvec_pointwise_acc(v, t_hat, y_hat, k, &ctx)
    poly_invntt(v, &ctx)
    POLY_ADD(v, v, e2, KYBER_N)

    kyber_byte_decode(mu_bits, m, 1)
    kyber_poly_decompress(mu, mu_bits, 1)
    POLY_ADD(v, v, mu, KYBER_N)
    kyber_reduce_poly(v)

    kyber_pack_ct(ct, u, v, param)
    return 0
```

Correctness notes:

- `kyber_byte_decode(mu_bits, m, 1)` treats 32 bytes as 256 one-bit coefficients.
- `kyber_poly_decompress(mu, mu_bits, 1)` maps bits to `0` or approximately `Q/2`.
- `t_hat` and `y_hat` are NTT domain; `u`, `v`, `e1`, `e2`, `mu` are normal domain after inverse NTT.

### 15.3 `kyber_pke_decrypt`

```text
kyber_pke_decrypt(m, dk_pke, ct, param):
    ctx = INIT_NTT_KYBER()
    k = param.k

    kyber_unpack_ct(u, v, ct, param)
    kyber_unpack_dk_pke(s_hat, dk_pke, param)

    kyber_polyvec_ntt(u, k, &ctx)
    kyber_polyvec_pointwise_acc(w, s_hat, u, k, &ctx)
    poly_invntt(w, &ctx)

    POLY_SUB(w, v, w, KYBER_N)
    kyber_reduce_poly(w)

    kyber_poly_compress(m_bits, w, 1)
    kyber_byte_encode(m, m_bits, 1)
    return 0
```

Correctness notes:

- `s_hat` is stored in NTT domain in `dk_pke`.
- `u` is decoded/decompressed into normal domain, then transformed by NTT before multiplication.

### 15.4 KEM scope

Current code intentionally stops at K-PKE:

```text
implemented:
    kyber_pke_keygen
    kyber_pke_encrypt
    kyber_pke_decrypt

not implemented in this phase:
    kyber_keygen_internal
    kyber_encaps_internal
    kyber_decaps_internal
    public Kyber KeyGen/Encaps/Decaps wrappers
```

The omitted layer is the CCA/KEM wrapper around K-PKE. It needs `H(ek)`,
`G(m || H(ek))`, `J(z || ct)`, re-encryption, and constant-time selection.
Those details should only be added after the bottom PKE layer remains stable.

- 第一版为了功能正确，可以先实现普通 compare/select，但文档和代码注释必须标出最终需要 constant-time。
- 不要在 ciphertext mismatch 时返回错误；必须输出 `K_bar`。

## 16. 测试函数设计

新增 `Test/src/test_kyber.c`，整体风格对齐 `test_dili.c`。

### 16.1 测试辅助函数

```c
void init_kyber_param(kyber_param *param, int mode);
const char* get_kyber_mode_name(int mode);
static int bytes_equal(const uint8_t *a, const uint8_t *b, int len);
static int poly_in_range(const poly a, int n, int min, int max);
static void fill_test_bytes(uint8_t *buf, int len, uint8_t seed);
```

模式命名：

```text
mode 2 -> "Kyber2"
mode 3 -> "Kyber3"
mode 5 -> "Kyber5"
```

### 16.2 `test_kyber_param_mode`

目标：确认参数长度没有写错。

```text
test_kyber_param_mode(mode):
    init param
    assert ek_bytes == 384*k + 32
    assert dk_pke_bytes == 384*k
    assert ct_bytes == 32*(du*k + dv)
    assert eta1 in {2,3}
    assert eta2 == 2
    assert du in {10,11}
    assert dv in {4,5}
```

### 16.3 `test_kyber_encode_decode`

目标：覆盖 `ByteEncode/ByteDecode` 的所有必要 `d`。

```text
for d in {1,4,5,10,11,12}:
    construct poly a:
        if d < 12: a[i] = deterministic value mod 2^d
        if d == 12: a[i] = deterministic value mod KYBER_Q
    kyber_byte_encode(buf, a, d)
    kyber_byte_decode(b, buf, d)
    assert b == a
```

额外 modulus check 测试：

```text
make a 384-byte public-key poly block containing one 12-bit value >= KYBER_Q
decode then encode
assert encoded block != original block
```

### 16.4 `test_kyber_compress_decompress`

目标：保证压缩/解压缩范围正确，并验证 FIPS 要求的 `Compress(Decompress(y)) == y`。

```text
for d in {1,4,5,10,11}:
    for y in 0..(1<<d)-1:
        x = kyber_decompress_one(y, d)
        y2 = kyber_compress_one(x, d)
        assert y2 == y

    construct poly a in [0,Q)
    kyber_poly_compress(c, a, d)
    assert c[i] in [0, 2^d-1]
    kyber_poly_decompress(b, c, d)
    assert b[i] in [0, Q-1]
```

### 16.5 `test_kyber_sample`

目标：检查采样范围、确定性和 nonce 区分。

```text
test_sample_ntt:
    seed = deterministic 34 bytes
    kyber_sample_ntt(a, seed, &ctx)
    kyber_sample_ntt(b, seed, &ctx)
    assert a == b
    assert all a[i] in [0,Q)

test_sample_cbd_eta2_eta3:
    buf = deterministic bytes
    kyber_sample_cbd(a, buf, 2)
    assert all a[i] in [0,2] or [Q-2,Q-1] representation
    kyber_sample_cbd(a, buf, 3)
    assert all a[i] in [0,3] or [Q-3,Q-1] representation

test_sample_noise_nonce:
    sigma = deterministic 32 bytes
    sample nonce 0 twice -> equal
    sample nonce 0 and 1 -> usually different; assert not equal for test vector seed
```

### 16.6 `test_kyber_ntt_mul`

目标：对已经确定正确的共享 32-bit Kyber NTT/乘法路径做防回归验证，避免 ML-KEM 上层误用普通点乘或重复 NTT。

```text
test_kyber_ntt_roundtrip:
    init ctx with INIT_NTT_KYBER
    generate random small poly a
    b = a
    poly_ntt(b)
    poly_invntt(b)
    convert/from Montgomery as current Poly code requires
    assert b == a mod Q

test_kyber_basemul_vs_naive:
    generate random small a,b
    a_hat = NTT(a)
    b_hat = NTT(b)
    c_hat = kyber polynomial NTT-domain multiply
    c = INTT(c_hat)
    c_naive = schoolbook negacyclic a*b mod (X^256+1, Q)
    assert c == c_naive mod Q
```

这个测试不是为了重写底层乘法，而是为了确认上层封装调用的是现有正确的 Kyber 乘法路径。测试失败时应优先检查上层域转换、Montgomery 表示、是否重复 NTT、是否误用 `poly_pointwise()`。

### 16.7 `test_kyber_pke_mode`

目标：先跑通 deterministic PKE 正确性。

```text
test_kyber_pke_mode(mode):
    init param
    d = deterministic 32 bytes
    r = deterministic 32 bytes
    m = deterministic 32 bytes

    kyber_pke_keygen(ek_pke, dk_pke, d, param)
    assert ek_pke not all zero
    assert dk_pke not all zero

    kyber_pke_encrypt(ct, ek_pke, m, r, param)
    assert ct not all zero

    kyber_pke_decrypt(m2, dk_pke, ct, param)
    assert m2 == m
```

### 16.8 `test_kyber_perf_mode`

第一阶段只做粗性能，不作为正确性阻塞：

```text
test_kyber_perf_mode(mode):
    bench kyber_pke_keygen
    bench kyber_pke_encrypt
    bench kyber_pke_decrypt
```

### 16.9 `main`

对齐 `test_dili.c`：

```text
main:
    print banner
    test_init()
    test_kyber_params()
    test_kyber_encode_decode()
    test_kyber_compress_decompress()
    test_kyber_sample()
    test_kyber_ntt_mul()
    test_kyber_pke()
    test_kyber_perf()
    print complete
```

## 17. 下一步实现顺序

为了降低调试难度，建议按这个顺序实现：

1. `kyber.h` 参数结构和 `init_kyber_param` 测试版本。
2. `kyber_pack.c`: encode/decode/compress/decompress。
3. `kyber_sample.c`: hash wrappers、PRF、SampleNTT、CBD。
4. `kyber_polyvec.c`: matrix expand、matrix/vector multiply、transpose multiply。
5. `kyber_pke_keygen`。
6. `kyber_pke_encrypt` 和 `kyber_pke_decrypt`。
7. Stop here for the current phase. KEM/internal/public Encaps/Decaps APIs are future work.
8. 最后再考虑公开 `ML-KEM.KeyGen/Encaps/Decaps` 输入检查 API。

这一顺序可以让每一步都有独立测试，并且能尽早发现 NTT/Montgomery 表示错误。
