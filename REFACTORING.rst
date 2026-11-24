第一轮重构说明
==============

本轮建立独立的素域、商环、NTT 计划和随机源接口，保留现有 Dilithium
及 Kyber K-PKE 的算法流程。当前算法代码仍是原有实现；回归一致不等于
已经完成 FIPS 203/204 符合性验证。Kyber 的完整 KEM 不在本轮范围内。

构建与验证
----------

在仓库根目录执行::

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure

也可以使用 ``cmake --build build --target run_test`` 运行全部正确性测试。
``test_poly --bench``、``test_ntt --bench``、``test_sig --bench`` 和
``test_enc --bench`` 显式开启性能测试；默认不运行耗时基准。

``PQC_NATIVE_OPTIMIZATION=ON`` 选择本机指令集；默认构建不强制启用。
GCC/Clang 可以通过 ``PQC_ENABLE_SANITIZERS=ON`` 开启 ASan/UBSan。
Linux 测试使用链接器包装 getrandom，注入 EINTR、短读及部分写入后失败。
Windows/POSIX 备用随机源需要在对应平台另行验证。

构建目标与依赖
--------------

* ``field``：fp32.c，仅依赖标准整数类型。
* ``algebra``：ring、poly、NTT、polyvec，仅链接 field。
* ``common``：系统随机源及 Keccak。
* ``core``：算法、兼容 NTT 参数预设和多项式采样，链接 algebra/common。

旧的 CMake 消费者可以继续链接 core。现在 core 是组合目标，手工链接
静态库时需同时提供其依赖库，不能只拷贝 libcore.a。
数学层可以独立链接，不会引入随机源、哈希或算法实现。

数学接口
--------

``fp32_ctx``（inc/poly/fp32.h）描述 3 到 INT32_MAX 的奇素数域。
``fp32_init`` 校验公开模数；失败返回 -1，原上下文不变。
``fp32_reduce_raw`` 等原始取模工具还允许正的合数/偶数模数，但并不因此
将这些对象当作域。一般规范化接口使用整数余数，是可移植参考后端；
它们未承诺所有处理器上的常数时间执行，也未完成侧信道审计。

``fp16_ctx``（inc/poly/fp16.h）是同一约定的 16 位版本，描述 3 到
INT16_MAX 的奇素数域，供系数以 int16_t 存储的算法复用。Montgomery
约简使用 R=2^16，``barrett`` 为 floor(2^15/q)，``mont_r2`` 为 R^2 mod q。
Aigis-Enc+ 的 ``reduce.c`` 只保留算法侧入口，约简运算全部转发到该算子。

``inc/poly/poly16.h`` 提供 int16_t 系数向量的通用算子：逐系数加/减、
Barrett 归约、``caddq``/``caddq2``、转入 Montgomery 域，以及
``fp16_poly_cbd`` 中心二项分布采样。采样器把 eta 拆成不超过 eta 的最大 2
的幂 hi 与 lo：前 n*hi/4 字节编码一个标准 CBD_hi，其后编码 CBD_lo 并与前者
逐系数相加。eta 为 1、2、4 时即普通连续 CBD；3、5、6 使用 ICCS NGCC 参考
实现的拆分布局，因此与参考输出逐位一致。

``inc/poly/ntt16.h`` 提供 128 点负循环 NTT（q=3329，官方 Kyber 旋转因子表，
逆变换自带 1/128 缩放，故 ``pqc_invntt128(pqc_ntt128(a)) == a``）。它与
``ntt_plan`` 的差别是系数为 int16_t、按 128 点分块，Aigis-Enc+ 的 512 项
多项式正好是 4 个块。

Aigis-Enc+ 的 ``ntt.c``、``cbd.c`` 与 ``poly.c`` 中对应的函数现在只转发到
上述算子；仍留在算法侧的是与规格绑定的位打包（压缩/解压缩、frommsg/tomsg、
tobytes/frombytes）、``gen_a`` 的拒绝采样，以及 NTT 域的 4 块 Karatsuba
乘法 ``poly_mont_mul``（含 NTT_Y 旋转因子表，随参数集与环布局变化）。

``inc/util/ct.h`` 提供常数时间等值比较与条件选择。Aigis-Enc+ 的 verify/cmov
只把 0/1 约定映射到这两个算子，因此 KEM 的密文校验与隐式拒绝共用同一实现。
两个函数都遍历全部字节，运行时间只与长度有关；它们不是编译器屏障，调用方
仍需避免把比较结果用于数据相关分支。

Montgomery 常量统一为正逆元 qinv = q^(-1) mod 2^32，使用减法形式。
``fp32_mont_reduce`` 是热路径内联内核，要求 |a| <= 2^31*q。
函数的输入范围、规范表示与延迟归约范围详见头文件。

``ring_ctx``（inc/poly/ring.h）描述 F_q[X]/(X^n+1)，包含域的值副本；
不包含 NTT 表，也不依赖 NTT 是否存在。初始化允许任意正维度。

``ntt_plan``（inc/poly/ntt.h）描述具体变换后端。本轮显式支持：

* q=8380417、n=256、8 层，标量乘法；
* q=8380417、n=256、7 层，二次块乘法；
* q=3329、n=256、7 层，二次块乘法。

``ntt_plan_init_256`` 拒绝其他组合。成功后将上下文视为只读；禁止直接修改
层数或模数。计划按值拥有 ring，表指向静态 const 数据，没有指向自身的
存储指针，因而可以安全复制，并在线程间只读共享。

例如::

    fp32_ctx field;
    ring_ctx ring;
    ntt_plan plan;
    if (fp32_init(&field, 3329) != 0 ||
        ring_init_negacyclic(&ring, &field, 256) != 0 ||
        ntt_plan_init_256(&plan, &ring, 7) != 0) {
        /* configuration error */
    }

NTT 保持 Montgomery 指数；变换域乘法使指数减一；逆变换使指数加一。
因此 INTT(NTT(a)) = R*a，而 INTT(MulNTT(NTT(a),NTT(b))) = a*b，
这里 R=2^32，等式均在商环中成立。

``poly_mul_ntt`` 根据计划选择标量/二次块乘法，返回规范余数并支持原地输出。
向量/矩阵乘法也从计划派生模式，不再接受 FULL/BASE 参数。
累加使用 int64，逆变换在蝶形运算前中心化输入，避免累加或表示差异导致溢出。

``poly256_i32`` 明确表示固定 256 项；``poly`` 和 ``lp_ctx`` 分别是
poly256_i32、ntt_plan 的过渡别名。标量 Montgomery 宏已由 fp32 API 取代。
旧 ``ctx.q/n/qc/layer`` 访问分别迁移为
``ctx.ring.field.q``、``ctx.ring.n``、``ctx.ring.field.qinv``、``ctx.layers``。
更改层数必须重新调用构造函数，不能只改字段。

``INIT_NTT_DILI`` / ``INIT_NTT_KYBER`` 作为 core 中的兼容预设保留。
普通多项式宏仍要求调用方保证 int32 算术范围和有效移位；所有 poly/vector
操作仍以 256 项为后端约束，ring 的维度可配置不代表这些内核支持任意维度。

随机源与采样迁移
----------------

生产 rand_seed 和全局伪随机状态已移除。默认使用 OS 安全随机源，
不会回退到时间或固定种子。``pqc_rng_func`` 由调用方提供状态，要求填满
所请求字节并返回 0；非零表示失败。``pqc_randombytes`` 在失败后清零
随机输出缓冲区，将错误统一为 -1，NULL 回调选择系统随机源。

新增 ``dili_keygen_with_seed``、``dili_keygen_with_rng``、
``kyber_pke_keygen_with_rng``；现有 ``dili_sign_with_rng`` 采用共同回调类型。
RNG 失败不会改写算法密钥或签名输出。dili_sign 的默认随机化行为保留。
旧 void 接口 dili_keygen、rand_bytes、rand_bits 在无法获得安全随机数时
会 abort；新调用者应采用可返回错误的接口。

数学层不再主动取得随机数。旧多项式采样声明迁入 inc/poly/poly_sample.h，
实现迁入 core；新增 poly_rand_with_rng 和 poly_uniform_CBD_with_rng。
确定性伪随机辅助函数只位于 Test/，不提供生产环境全局切换。

测试边界与后续工作
------------------

GF(16) 算术现位于 ``inc/gf/gf16.h`` / ``src/gf/gf16.c``，采用
``GF(2)[x]/(x^4+x+1)`` 的低四位多项式表示。矩形矩阵乘法、异或累加及
1 到 5 阶行列式位于 ``inc/gf/gf16_matrix.h`` / ``src/gf/gf16_matrix.c``。
GF(31) 的域运算、逆元表、向量算子与 ``GF(31^2)=GF(31)[x]/(x^2-3)`` 的
平面布局矩阵乘法位于 ``inc/gf/gf31.h``、``inc/gf/gf31_matrix.h`` 和
``src/gf/gf31.c``、``src/gf/gf31_matrix.c``，对应
``operator_gf31_correctness`` 与 ``gf31_*`` 性能用例。
2 的幂模数算术位于 ``inc/gf/mod_pow2.h`` / ``src/gf/mod_pow2.c``
（``pqc_mod_pow2``、``pqc_round_pow2``），对应
``operator_mod_pow2_correctness`` 与 ``mod_pow2_*`` 性能用例。
``q = 4171777``、``n = 512`` 的环算术位于 ``inc/poly/ring4171.h`` /
``src/poly/ring4171.c``：Montgomery/Barrett/通用约减、规范与中心表示、
512 点正反 NTT、Montgomery 逐点乘、范数检查及系数级加减/缩放。官方
zeta 表随之迁入库内。Aigis-Sig+ 的 ``poly.c`` 只剩算法侧的打包、采样与
类型适配，``reduce.c``、``ntt.c`` 已删除；对应
``operator_ring4171_correctness`` 与 ``q4171_*`` 性能用例。
各算法模块的底层依赖现状：

* DOVE、VDOO、Aigis-Enc+、Aigis-Sig+、Scloud+、TSUOV、FlexTree 均已
  不携带自己的 GF/AES/SHA 实现，统一调用 ``inc/gf``、``inc/sym``、
  ``inc/ngcc``、``inc/util`` 与 ``inc/tree``。
* Scloud+ 的采样与 BW/RBW 编解码仍留在 ``reference``，属于参数专用编码；
  其 608 维矩阵加减已改用库内 2 的幂模算子。
* FlexTree 的 ``reference/params/pors_precomputed-*.h`` 是每个参数集的
  预计算 PORS 表，属于算法数据而非算子实现。
* SNOVA R2/R3、MAYO、MQOM 使用上游提交的 AES/Keccak/Rijndael 优化内核
  （含 AVX2/NEON 汇编）；这些内核与其 SIMD 数据布局绑定，替换会改变
  性能特征，暂按“不破坏现有实现”保留，后续可先替换其通用 C 分支。
  SNOVA 的算法实现、签名 API、对称层及平台变体位于 ``alg/snova_r3``；
  CMake 兼容目标 ``pqc_snova`` 指向 ``pqc_snova_r3`` 并链接 ``pqc``。
  当前仅接入用户优化实现：x86-64 默认 AVX2，AArch64 默认 NEON。公共
  标量算子按需复用；扩展表示及 SIMD 热路径仍留在 SNOVA 后端。GF(16)
  求解器和其他多变量专用运算尚未拆分。

NGCC 算法位于 ``alg/ngcc_*``，包含参数、密钥和密文或签名格式、算法流程及
公共 API 适配。哈希和确定性随机数统一调用 ``inc/ngcc/auxfunc.h`` 和
``inc/ngcc/drng.h``。底层算子按表示和用途放在共享模块中。

- DOVE 复用多项式基 GF(256) 内核及矩阵乘法，矩阵向量积、求逆、求解和
  上三角布局转换位于 ``inc/gf`` 和 ``src/gf/common``，签名实现直接调用这些接口。
- VDOO 的 GF(16)/GF(256) 字节使用 GF(4) 塔域坐标。算术、向量乘加和批量矩阵
  运算位于 ``inc/gf/composite*.h`` 和 ``src/gf/tower``，以兼容提交版 KAT。
- Aigis-Enc+ 的 KEM 和公钥加密流程直接调用 ``poly16``、``ntt16``、
  ``ngcc`` 和 ``util/ct.h``。系数编解码和重复位消息编解码统一由
  ``src/poly/poly16_pack.c`` 实现。密钥与密文布局、种子派生和噪声参数
  在算法层定义。
- Aigis-Sig+ 的哈希与扩展复用 ``inc/ngcc`` 的 pseudoXOF/pseudohash/sm3hash
  和 DRNG。其环参数为 q=4171777、n=512，与 Kyber/Dilithium 的
  q=3329 不同，因此 NTT、多项式算术和约减保留在
  ``alg/ngcc_aigis_sig/reference``，未并入 ``poly16``/``ntt16``。
- Scloud+ 的公开矩阵 A 用 AES-128-CTR 逐行展开，接口保留在
  ``alg/ngcc_scloud/aes_shim.c``，实际分组加密复用 ``src/sym/aes.c``
  （纯 C 或 AES-NI 后端），输出与参考实现逐字节一致。F/G/H/K 是带单字节
  标签的 SHAKE256/SHA3-512，直接构建在 ``inc/sym/keccak.h`` 上：增量
  absorb/finalize/squeeze/squeezeblocks 只是对库内海绵接口的薄封装，参考
  实现自带的 TweetFips202 副本已删除。常量时间比较、条件选择和安全清零
  分别复用 ``util/ct.h`` 与 ``util/secure_clear.h``。采样分布为 BD4/BD2
  （128 级），矩阵乘法、BW/RBW 消息编解码和 Pack10/ShortPack 布局按参考
  实现保留，因为其维度为 608、模数 q=1024，与现有 ``poly16``/``ntt16``
  不兼容。评审报告 kem-35-1（重加密拒绝采样泄漏解密消息类别）只影响
  256/384/512 级的 BD6/BD12 拒绝采样；128 级是非拒绝采样对照实现，不受影响。
- TSUOV 的哈希、伪哈希与 XOF 调用 ``inc/ngcc/auxfunc.h`` 的
  sm3hash/pseudohash/pseudoXOF，随机源使用 ``inc/ngcc/drng.h``，安全清零
  复用 ``util/secure_clear.h``。域元素与签名的 5 比特打包、Scloud+ 的
  消息块打包共用 ``inc/util/bits.h`` 的 LSB-first 位打包算子。域运算、
  逆元、向量算子与 L=2 平面塔的矩阵乘法和 phi 坐标映射都已移入库内
  ``inc/gf/gf31*.h`` / ``src/gf/gf31*.c``，``alg/ngcc_tsuov/reference``
  只保留 E 算子、扩展与拒绝采样、UOV 消元求解和签名编码。评审报告
  sign-31 的处理：
  报告缺陷 sign-31-3/5（分配失败复用旧消息代表、长度加法溢出）、sign-31-4
  （忽略声明长度）、sign-31-6（非规范 5 比特编码）和 sign-31-8（忽略 XOF
  失败）已在参考源码中修复，签名格式不变，官方 KAT 仍逐字节一致；
  sign-31-1、sign-31-7 属于规范/证明问题，sign-31-2 是参考实现的局部
  计时侧信道，均在 ``alg/ngcc_tsuov/adapter.c`` 顶部记录。
- FlexTree 的完整树与截短底层树共用 ``inc/tree/merkle.h`` 和
  ``src/tree/merkle.c``，覆盖树构建、单叶认证路径、多叶认证节点位置与根重建。
  PORS 的组合数查表、字典序组合解码和按字节查表模约减调用
  ``util/combinatorics.h``。算法层保留域分离地址、哈希输入格式、WOTS 链、
  签名计数器搜索和 Flextree-160f 参数表。WOTS 叶生成与验签位于同一个
  ``wots.c``，计数器编解码直接调用 ``util/bytes.h``。

NGCC 正确性测试覆盖提交版 KAT、公共 API 及多算法链接。算子测试覆盖塔域参考
乘法、矩阵求解及逆矩阵、分块 NTT 与朴素卷积的对照、位编码及大整数边界。
Merkle 算子与递归构造的树对照，对 1 到 9 个叶子的全部非空叶集合验证多叶认证；
组合解码与小参数的全部字典序组合枚举对照。对应
基准由 CMake 注册，打开 ``PQC_PERF_TESTS`` 后可通过 ``run_performance``
或 CTest 的 ``performance`` 标签运行。

test_fp32 使用独立取模/Euclid 参考和极值检查；test_fp16 使用相同的参考
策略覆盖 16 位 Montgomery/Barrett 与条件加法。test_poly16 用逐位参考验证
16 位多项式算术与 eta 1 到 6 的 CBD；test_ntt16 用直接按定义求值的 NTT
（psi^(2*bitrev(j)+1)）验证 128 点变换、往返一致性与分块包装。
test_contexts 只链接 algebra，验证构造失败、复制、每个旋转因子、缩放、
朴素负循环乘法、原地输出、矩形矩阵及长向量累加。
test_rng 检查注入、系统接口故障传播及输出不变契约。

Test/regression/test_regression.c 中的 18 个 SHA3-256 摘要从重构前的代码
单独构建捕获：覆盖三个 Kyber 参数集的公钥/私钥/密文，以及三个 Dilithium
参数集的公钥/私钥有效部分/签名。种子和签名随机量固定，哈希不包含
结构体填充或未使用槽位。这是行为回归数据，不是官方标准 KAT。

后续再处理：任意维度缓冲区视图、算法参数集封装、完整字节串 API、
矩阵扩展域分离从 polyvec 迁回算法、GF(16) 消元/求解器及其他多变量专用
算子、经验证的常数时间优化后端，以及标准向量/侧信道验证。

GF 与多算法接入（2026-11-05）
----------------------------

当前构建以 CMakeLists.txt、cmake/options.cmake、cmake/arch.cmake、
cmake/gf.cmake、cmake/algorithms.cmake、cmake/tests.cmake 和
cmake/install.cmake 为准，配置头由 inc/core/pqc_conf.h.in 生成。
GF 运算独立为 pqc_gf；签名目标为 pqc_snova_r2、pqc_snova_r3、
pqc_mayo、pqc_mqom，pqc_snova 保留为 R3 兼容目标。
共享库与静态库由同一 OBJECT 库生成，默认同时构建，并支持安装后
find_package(pqc)。
四算法公共 signed-message 测试位于 Test/sig/test_signature_api.c。
源码来源、后端选择、参数范围与验证边界见 docs/gf-integration/README.md
及 VALIDATION.md。本文前述 field/algebra/common/core 是早期重构说明，
并非当前实际链接目标。

构建入口已移除 ref/portable 基线选择；GF 和四算法仅选择用户的优化后端。
公共标量运算和优化实现的内部依赖保留，所有四算法公共测试继续位于 Test/sig。
