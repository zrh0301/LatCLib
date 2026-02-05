# Dilithium 算法详细规范

## 参数定义

### Dilithium2 参数集
- `k = 4`, `l = 4` (矩阵维度)
- `η = 2` (密钥系数范围 [-2, 2])
- `d = 13` (t 的分解位数)
- `τ = 39` (挑战多项式中 ±1 的个数)
- `γ₁ = 2^17 = 131072` (掩码 y 的范围)
- `γ₂ = (q-1)/88 = 95232` (舍入参数)
- `ω = 80` (hint 的最大权重)
- `β = 78` (拒绝采样界)
- `q = 8380417` (模数)
- `n = 256` (多项式次数)

### Dilithium3 参数集
- `k = 6`, `l = 5`
- `η = 4`, `d = 13`, `τ = 49`
- `γ₁ = 2^19 = 524288`
- `γ₂ = (q-1)/32 = 261888`
- `ω = 55`, `β = 196`

### Dilithium5 参数集
- `k = 8`, `l = 7`
- `η = 2`, `d = 13`, `τ = 60`
- `γ₁ = 2^19 = 524288`
- `γ₂ = (q-1)/32 = 261888`
- `ω = 75`, `β = 120`

---

## 1. 密钥生成 (KeyGen)

### **步骤 1.1: 生成随机种子**

**输入**: 无（或系统熵源）

**处理**:
```
ζ ← {0,1}^256                // 256位随机种子
```

**输出**:
- `ζ`: 256位 (32字节) 随机比特串

---

### **步骤 1.2: 派生种子**

**输入**: `ζ` (256位)

**处理**:
```
(ρ, ρ', K) ← SHAKE256(ζ, 96)    // 扩展到96字节
ρ = [0:32]                        // 取前32字节
ρ' = [32:64]                      // 取中间32字节
K = [64:96]                       // 取最后32字节
```

**输出**:
- `ρ`: 256位 (32字节) - 用于展开矩阵 A
- `ρ'`: 256位 (32字节) - 用于采样 s1, s2
- `K`: 256位 (32字节) - 签名时的密钥

---

### **步骤 1.3: 展开矩阵 A**

**输入**: `ρ` (32字节)

**处理**:
```
for i = 0 to k-1:
    for j = 0 to l-1:
        seed = ρ || byte(i) || byte(j)    // 34字节
        stream = SHAKE128(seed, ∞)         // 无限流
        A[i][j] = SampleUniform(stream)    // 采样256个系数
```

**SampleUniform 详细处理**:
```
输入: stream (SHAKE128 输出流)
输出: poly ∈ R_q (256个系数)

ctr = 0
pos = 0
while ctr < 256:
    t = stream[pos] | (stream[pos+1] << 8) | (stream[pos+2] << 16)
    t = t & 0x7FFFFF                       // 取23位
    pos += 3
    if t < q:                               // 拒绝采样
        poly[ctr] = t
        ctr++

每个系数编码: 23位 → 3字节/系数 (约)
总计: 256系数 × 3字节 ≈ 768字节/多项式
矩阵 A: k×l 个多项式
```

**输出**:
- `A`: k×l 矩阵，每个元素是 R_q 中的多项式
- 存储大小: `k × l × 256系数 × 4字节 = k × l × 1KB`
- Dilithium2: 4×4×1KB = 16KB (运行时)

---

### **步骤 1.4: 采样密钥 s1, s2**

**输入**: `ρ'` (32字节), `η` (参数)

**处理**:
```
// 采样 s1 ∈ R_q^l
for i = 0 to l-1:
    seed = ρ' || byte(i) || byte(0)
    stream = SHAKE256(seed, ∞)
    s1[i] = SampleEta(stream, η)

// 采样 s2 ∈ R_q^k  
for i = 0 to k-1:
    seed = ρ' || byte(l+i) || byte(0)
    stream = SHAKE256(seed, ∞)
    s2[i] = SampleEta(stream, η)
```

**SampleEta 详细处理 (η=2)**:
```
输入: stream, η=2
输出: poly ∈ R_q, 系数 ∈ [-2, 2]

ctr = 0
pos = 0
while ctr < 256:
    t0 = stream[pos] & 0x0F           // 低4位
    t1 = stream[pos] >> 4              // 高4位
    pos++
    
    if t0 < 15:
        t0 = t0 - (205*t0 >> 10)*5     // 模5运算
        poly[ctr++] = 2 - t0           // 范围 [-2,2]
    
    if t1 < 15 and ctr < 256:
        t1 = t1 - (205*t1 >> 10)*5
        poly[ctr++] = 2 - t1

编码效率: 4位/系数 (约)
总需字节: 256系数 × 0.5字节 = 128字节/多项式
```

**SampleEta 详细处理 (η=4)**:
```
if t0 < 9:
    poly[ctr++] = 4 - t0               // 范围 [-4,4]
if t1 < 9:
    poly[ctr++] = 4 - t1
```

**输出**:
- `s1`: l 个多项式，每个系数 ∈ [-η, η]
- `s2`: k 个多项式，每个系数 ∈ [-η, η]
- 存储: `(k+l) × 256系数 × 1字节 ≈ (k+l) × 256B`
- Dilithium2: 8 × 256B = 2KB

---

### **步骤 1.5: 计算 t = A·s1 + s2**

**输入**: `A` (k×l矩阵), `s1` (l向量), `s2` (k向量)

**处理**:
```
// 转到 NTT 域
for i = 0 to l-1:
    s1_ntt[i] = NTT(s1[i])

for i = 0 to k-1:
    for j = 0 to l-1:
        A_ntt[i][j] = NTT(A[i][j])

// 矩阵乘法 (NTT 域中逐系数相乘)
for i = 0 to k-1:
    t_ntt[i] = 0
    for j = 0 to l-1:
        for coef = 0 to 255:
            t_ntt[i][coef] += Montgomery(A_ntt[i][j][coef] * s1_ntt[j][coef], q)

// 转回时域并加 s2
for i = 0 to k-1:
    t[i] = INTT(t_ntt[i])
    for coef = 0 to 255:
        t[i][coef] = (t[i][coef] + s2[i][coef]) mod q
```

**NTT/INTT 规格**:
- NTT: 8层蝶形运算，O(n log n)
- Montgomery 约减: `(a*b) mod q` 用 32位×32位 → 64位 → 32位
- 每系数范围: [0, q-1]

**输出**:
- `t`: k 个多项式，每个系数 ∈ [0, q-1]
- 存储: `k × 256 × 4字节 = k × 1KB`

---

### **步骤 1.6: Power2Round 分解 t**

**输入**: `t` (k向量), `d = 13`

**处理**:
```
for i = 0 to k-1:
    for coef = 0 to 255:
        r = t[i][coef]
        t1[i][coef] = (r + 2^(d-1)) >> d        // 高位
        t0[i][coef] = r - (t1[i][coef] << d)    // 低位
```

**位分解说明**:
- `t = t1·2^d + t0`
- `t1`: 高 (32-d) 位，范围 [0, 2^(32-d)-1]
- `t0`: 低 d 位，范围 [-2^(d-1), 2^(d-1))
- d=13: `t1` 19位, `t0` 13位

**输出**:
- `t1`: k 个多项式，每个系数 19位 (Dilithium2)
- `t0`: k 个多项式，每个系数 13位
- 打包大小:
  - `t1`: 256系数 × 19位 / 8 = 608字节/多项式
  - `t0`: 256系数 × 13位 / 8 = 416字节/多项式

---

### **步骤 1.7: 计算公钥哈希 tr**

**输入**: `ρ` (32字节), `t1` (k向量)

**处理**:
```
pk_bytes = Pack(ρ, t1)        // 打包公钥
tr = SHAKE256(pk_bytes, 64)   // 512位哈希
```

**Pack(ρ, t1) 格式**:
```
输出: 32 + k×320 字节 (Dilithium2)

[ρ: 32字节]
[t1[0]: 320字节 = 256系数×10位/8]
[t1[1]: 320字节]
...
[t1[k-1]: 320字节]

注: 实际 t1 每系数10位压缩编码 (去掉高位0)
```

**输出**:
- `tr`: 512位 (64字节) - 公钥的哈希摘要

---

### **步骤 1.8: 输出密钥对**

**公钥 (pk)**:
```
pk = (ρ, t1)
大小: 32 + k×320 字节
Dilithium2: 32 + 4×320 = 1312字节
```

**私钥 (sk)**:
```
sk = (ρ, K, tr, s1, s2, t0)

打包格式:
[ρ: 32字节]
[K: 32字节]
[tr: 64字节]
[s1: l×128字节, η位编码]
[s2: k×128字节, η位编码]
[t0: k×416字节, 13位编码]

Dilithium2: 32+32+64 + 4×128 + 4×128 + 4×416 = 2560字节
```

---

## 2. 签名 (Sign)

### **步骤 2.1: 计算消息哈希 μ**

**输入**: `M` (任意长度消息), `tr` (64字节)

**处理**:
```
μ = SHAKE256(tr || M, 64)     // 512位哈希
```

**输出**:
- `μ`: 64字节 - 绑定公钥和消息的哈希

---

### **步骤 2.2: 初始化拒绝采样计数器**

**输入**: `K` (32字节), `μ` (64字节)

**处理**:
```
κ = 0                          // 拒绝采样计数器
ρ' = K                         // 初始种子
```

---

### **步骤 2.3: 拒绝采样循环开始**

**标签**: `REJECTION_LOOP`

---

### **步骤 2.3.1: 更新种子**

**输入**: `K`, `ρ'`, `μ`, `κ`

**处理**:
```
seed = K || ρ' || μ
ρ' = SHAKE256(seed, 64)        // 512位新种子
```

**输出**:
- `ρ'`: 64字节 - 用于展开掩码 y

---

### **步骤 2.3.2: 展开掩码向量 y**

**输入**: `ρ'` (64字节), `κ`, `γ₁`

**处理**:
```
for i = 0 to l-1:
    seed = ρ' || uint16(κ) || byte(i)
    stream = SHAKE256(seed, ∞)
    y[i] = SampleGamma1(stream, γ₁)
```

**SampleGamma1 详细处理 (γ₁ = 2^17)**:
```
输入: stream, γ₁ = 131072
输出: poly, 系数 ∈ [-(γ₁-1), γ₁]

ctr = 0
pos = 0
while ctr < 256:
    // 18位编码
    t = stream[pos] | (stream[pos+1]<<8) | (stream[pos+2]<<16)
    t = t & 0x3FFFF              // 18位掩码
    pos += 2.25字节 (实际按字节对齐)
    
    poly[ctr++] = t - γ₁         // 转换到 [-(γ₁-1), γ₁]

每系数: 18位
总需: 256×18位/8 = 576字节/多项式
```

**SampleGamma1 详细处理 (γ₁ = 2^19)**:
```
// 20位编码
t = stream[pos] | (stream[pos+1]<<8) | (stream[pos+2]<<16)
t = t & 0xFFFFF                 // 20位掩码
poly[ctr++] = t - γ₁

每系数: 20位
总需: 256×20位/8 = 640字节/多项式
```

**输出**:
- `y`: l 个多项式，系数 ∈ [-(γ₁-1), γ₁]
- 存储: `l × 640字节` (γ₁=2^19)

---

### **步骤 2.3.3: 计算 w = A·y**

**输入**: `A` (k×l), `y` (l向量)

**处理**:
```
// NTT 域乘法
for i = 0 to l-1:
    y_ntt[i] = NTT(y[i])

for i = 0 to k-1:
    w_ntt[i] = 0
    for j = 0 to l-1:
        for coef = 0 to 255:
            w_ntt[i][coef] += Montgomery(A_ntt[i][j][coef] * y_ntt[j][coef], q)

for i = 0 to k-1:
    w[i] = INTT(w_ntt[i])
    Reduce(w[i])                 // 约减到 [0, q-1]
```

**输出**:
- `w`: k 个多项式，系数 ∈ [0, q-1]

---

### **步骤 2.3.4: 提取高位 w1 = HighBits(w, 2γ₂)**

**输入**: `w` (k向量), `α = 2γ₂`

**处理**:
```
for i = 0 to k-1:
    for coef = 0 to 255:
        r = w[i][coef] mod^± q     // 中心化到 [-(q-1)/2, (q-1)/2]
        
        r1 = r / α                  // 除法取整
        r0 = r - r1 * α
        
        if r0 > α/2:
            r1 = r1 + 1
            r0 = r0 - α
        
        w1[i][coef] = r1 mod^+ (q/α)  // 约减到 [0, q/α)
```

**位宽说明 (Dilithium2, γ₂=95232)**:
- `α = 2γ₂ = 190464`
- `q/α ≈ 44` → w1 每系数 6位
- 打包: 256系数 × 6位 / 8 = 192字节/多项式

**位宽说明 (Dilithium3/5, γ₂=261888)**:
- `α = 2γ₂ = 523776`
- `q/α = 16` → w1 每系数 4位
- 打包: 256系数 × 4位 / 8 = 128字节/多项式

**输出**:
- `w1`: k 个多项式，每系数 4或6位
- `w0`: k 个多项式 (后续用) ，系数 ∈ [-(α/2), α/2]

---

### **步骤 2.3.5: 计算挑战 c̃**

**输入**: `μ` (64字节), `w1` (k向量)

**处理**:
```
w1_bytes = Pack(w1)              // 打包 w1
c_seed = μ || w1_bytes
c_hash = SHAKE256(c_seed, 32)    // 256位哈希
c̃ = c_hash                       // 保存原始哈希
c = SampleInBall(c_hash, τ)     // 展开挑战多项式
```

**SampleInBall 详细处理 (τ=39)**:
```
输入: seed (32字节), τ
输出: c ∈ R_q, 恰好 τ 个系数为 ±1, 其余为 0

初始化: c[0..255] = 0
signs = uint64(seed[0..7])       // 前8字节作为符号位

pos = 8
for i = 256-τ to 255:
    // 生成随机位置 j ∈ [0, i]
    j = seed[pos++] mod (i+1)
    while j already used:
        j = seed[pos++] mod (i+1)
    
    // 交换 c[i] 和 c[j]
    c[i] = c[j]
    
    // 设置 c[j] = ±1
    sign_bit = (signs >> (i-256+τ)) & 1
    c[j] = 1 - 2*sign_bit           // +1 或 -1

输出: 稀疏多项式, τ 个 ±1
编码: 256位位图 + τ 位符号 ≈ 32+5 = 37字节
```

**输出**:
- `c̃`: 32字节 - 挑战哈希（签名的一部分）
- `c`: 多项式，τ 个系数为 ±1

---

### **步骤 2.3.6: 计算 z = y + c·s1**

**输入**: `y` (l向量), `c` (多项式), `s1` (l向量)

**处理**:
```
c_ntt = NTT(c)

for i = 0 to l-1:
    s1_ntt[i] = NTT(s1[i])
    cs1_ntt[i] = PointwiseMul(c_ntt, s1_ntt[i])
    cs1[i] = INTT(cs1_ntt[i])
    
    for coef = 0 to 255:
        z[i][coef] = y[i][coef] + cs1[i][coef]
        z[i][coef] = z[i][coef] mod^± q      // 中心化到 [-(q-1)/2, (q-1)/2]
```

**输出**:
- `z`: l 个多项式

---

### **步骤 2.3.7: 检查 ||z||_∞ < γ₁ - β**

**输入**: `z` (l向量), `γ₁`, `β`

**处理**:
```
bound = γ₁ - β

for i = 0 to l-1:
    for coef = 0 to 255:
        val = z[i][coef] mod^± q
        if |val| >= bound:
            κ++
            goto REJECTION_LOOP      // 拒绝，重新采样
```

**检查阈值** (Dilithium2):
- `γ₁ - β = 131072 - 78 = 130994`
- 失败概率: 约 1-2%

---

### **步骤 2.3.8: 计算 r0 = LowBits(w - c·s2, 2γ₂)**

**输入**: `w` (k向量), `c`, `s2` (k向量), `α = 2γ₂`

**处理**:
```
c_ntt = NTT(c)

for i = 0 to k-1:
    s2_ntt[i] = NTT(s2[i])
    cs2_ntt[i] = PointwiseMul(c_ntt, s2_ntt[i])
    cs2[i] = INTT(cs2_ntt[i])
    
    for coef = 0 to 255:
        tmp = w[i][coef] - cs2[i][coef]
        tmp = tmp mod^± q
        
        // 提取低位
        r1 = tmp / α
        r0[i][coef] = tmp - r1 * α
        
        if r0[i][coef] > α/2:
            r0[i][coef] -= α
```

**输出**:
- `r0`: k 个多项式，系数 ∈ [-(α/2), α/2]

---

### **步骤 2.3.9: 检查 ||r0||_∞ < γ₂ - β**

**输入**: `r0` (k向量), `γ₂`, `β`

**处理**:
```
bound = γ₂ - β

for i = 0 to k-1:
    for coef = 0 to 255:
        if |r0[i][coef]| >= bound:
            κ++
            goto REJECTION_LOOP
```

**检查阈值** (Dilithium2):
- `γ₂ - β = 95232 - 78 = 95154`

---

### **步骤 2.3.10: 计算 c·t0**

**输入**: `c`, `t0` (k向量)

**处理**:
```
c_ntt = NTT(c)

for i = 0 to k-1:
    t0_ntt[i] = NTT(t0[i])
    ct0_ntt[i] = PointwiseMul(c_ntt, t0_ntt[i])
    ct0[i] = INTT(ct0_ntt[i])
    Reduce(ct0[i])
```

**输出**:
- `ct0`: k 个多项式

---

### **步骤 2.3.11: 检查 ||c·t0||_∞ < γ₂**

**输入**: `ct0` (k向量), `γ₂`

**处理**:
```
for i = 0 to k-1:
    for coef = 0 to 255:
        val = ct0[i][coef] mod^± q
        if |val| >= γ₂:
            κ++
            goto REJECTION_LOOP
```

---

### **步骤 2.3.12: 生成 hint h**

**输入**: `ct0` (k向量), `w` (k向量), `cs2` (k向量), `α = 2γ₂`

**处理**:
```
for i = 0 to k-1:
    hint_count = 0
    for coef = 0 to 255:
        z_val = -ct0[i][coef]
        r_val = w[i][coef] - cs2[i][coef] + ct0[i][coef]
        r_val = r_val mod^± q
        
        // MakeHint 逻辑
        r1 = HighBits(r_val, α)
        v1 = HighBits(r_val + z_val, α)
        
        if r1 != v1:
            h[i][coef] = 1
            hint_count++
        else:
            h[i][coef] = 0
    
    if hint_count > ω:
        κ++
        goto REJECTION_LOOP
```

**Hint 编码**:
```
每个多项式: 位图 256位 + 权重信息
紧凑编码: 记录 1 的位置索引
最大 ω 个索引 × 8位 = ω 字节
Dilithium2: 最多 80字节/多项式
```

**输出**:
- `h`: k 个 hint 多项式（稀疏）
- 总 hint 权重 ≤ ω

---

### **步骤 2.4: 签名输出**

**签名格式**:
```
σ = (c̃, z, h)

[c̃: 32字节]
[z: l 个多项式, 每个 (γ₁-β) 编码]
[h: k 个 hint, 紧凑编码]

编码细节:
- z 系数范围: [-(γ₁-β), γ₁-β]
- 编码位宽: 18位 (Dilithium2) 或 20位 (Dilithium3/5)
- z 打包: l × 256 × 18位/8 = l × 576字节

Dilithium2 签名大小:
32 + 4×576 + 640(hint) ≈ 2420字节
```

---

## 3. 验证 (Verify)

### **步骤 3.1: 解析签名**

**输入**: `σ` (签名), `pk` (公钥), `M` (消息)

**处理**:
```
(c̃, z, h) = Unpack(σ)
(ρ, t1) = Unpack(pk)
```

**输出**:
- `c̃`: 32字节
- `z`: l 个多项式
- `h`: k 个 hint
- `ρ`: 32字节
- `t1`: k 个多项式

---

### **步骤 3.2: 检查 ||z||_∞ < γ₁ - β**

**输入**: `z` (l向量)

**处理**:
```
bound = γ₁ - β
for i = 0 to l-1:
    for coef = 0 to 255:
        if |z[i][coef]| >= bound:
            return REJECT
```

---

### **步骤 3.3: 检查 hint 权重**

**输入**: `h` (k向量)

**处理**:
```
total_weight = 0
for i = 0 to k-1:
    for coef = 0 to 255:
        total_weight += h[i][coef]

if total_weight > ω:
    return REJECT
```

---

### **步骤 3.4: 计算 tr 和 μ**

**输入**: `ρ`, `t1`, `M`

**处理**:
```
pk_bytes = Pack(ρ, t1)
tr = SHAKE256(pk_bytes, 64)
μ = SHAKE256(tr || M, 64)
```

---

### **步骤 3.5: 展开 A 和 c**

**输入**: `ρ`, `c̃`

**处理**:
```
A = ExpandA(ρ)                   // 同 KeyGen 步骤
c = SampleInBall(c̃, τ)          // 同 Sign 步骤
```

---

### **步骤 3.6: 计算 w' = A·z - c·t1·2^d**

**输入**: `A` (k×l), `z` (l), `c`, `t1` (k), `d`

**处理**:
```
// Part 1: A·z
z_ntt = NTT(z)
Az_ntt = MatMul(A_ntt, z_ntt)
Az = INTT(Az_ntt)

// Part 2: c·t1·2^d
c_ntt = NTT(c)
for i = 0 to k-1:
    t1_shifted[i] = t1[i] << d       // 乘 2^d
    t1_ntt[i] = NTT(t1_shifted[i])
    ct1_ntt[i] = PointwiseMul(c_ntt, t1_ntt[i])
    ct1[i] = INTT(ct1_ntt[i])

// 相减
for i = 0 to k-1:
    for coef = 0 to 255:
        w_prime[i][coef] = Az[i][coef] - ct1[i][coef]
        w_prime[i][coef] = w_prime[i][coef] mod q
```

**输出**:
- `w'`: k 个多项式

---

### **步骤 3.7: 使用 hint 恢复 w1'**

**输入**: `h` (k向量), `w'` (k向量), `α = 2γ₂`

**处理**:
```
for i = 0 to k-1:
    for coef = 0 to 255:
        r = w_prime[i][coef] mod^± q
        
        // 提取高位
        r1 = r / α
        r0 = r - r1 * α
        if r0 > α/2:
            r1++
            r0 -= α
        
        // 应用 hint
        if h[i][coef] == 1:
            if r0 > 0:
                r1 = (r1 + 1) mod (q/α)
            else:
                r1 = (r1 - 1) mod (q/α)
        
        w1_prime[i][coef] = r1
```

**输出**:
- `w1'`: k 个多项式（高位部分）

---

### **步骤 3.8: 重新计算挑战 c̃'**

**输入**: `μ`, `w1'`

**处理**:
```
w1_prime_bytes = Pack(w1_prime)
c_seed = μ || w1_prime_bytes
c_tilde_prime = SHAKE256(c_seed, 32)
```

**输出**:
- `c̃'`: 32字节

---

### **步骤 3.9: 比较挑战**

**输入**: `c̃`, `c̃'`

**处理**:
```
if c̃ == c̃':
    return ACCEPT
else:
    return REJECT
```

---

## 附录：数据大小汇总

### Dilithium2

| 项目 | 大小 | 说明 |
|------|------|------|
| 公钥 (pk) | 1312字节 | ρ(32) + t1(4×320) |
| 私钥 (sk) | 2560字节 | ρ(32) + K(32) + tr(64) + s1(512) + s2(512) + t0(1664) |
| 签名 (σ) | 2420字节 | c̃(32) + z(2304) + h(640) |
| 中间矩阵 A | 16KB (运行时) | 4×4×256×4字节 |

### Dilithium3

| 项目 | 大小 | 说明 |
|------|------|------|
| 公钥 (pk) | 1952字节 | ρ(32) + t1(6×320) |
| 私钥 (sk) | 4032字节 | 更大的 k,l,η |
| 签名 (σ) | 3293字节 | 更大的 γ₁ |

### Dilithium5

| 项目 | 大小 | 说明 |
|------|------|------|
| 公钥 (pk) | 2592字节 | ρ(32) + t1(8×320) |
| 私钥 (sk) | 4864字节 | 最大的 k,l |
| 签名 (σ) | 4595字节 | 最大的签名 |

---

## 性能特征

### 拒绝采样预期次数

| 参数集 | 平均重试次数 | 说明 |
|--------|-------------|------|
| Dilithium2 | ~1.5 次 | β 较小，拒绝概率高 |
| Dilithium3 | ~1.8 次 | γ₁ 更大，成功率提高 |
| Dilithium5 | ~1.3 次 | 参数优化 |

### 计算复杂度

| 操作 | 复杂度 | 次数 |
|------|--------|------|
| NTT | O(n log n) | KeyGen: (k×l + k+l), Sign: (3l+4k), Verify: (2l+3k) |
| PointwiseMul | O(n) | 每次矩阵乘: k×l 次 |
| SHAKE | 可变 | 多次（种子展开、哈希） |
