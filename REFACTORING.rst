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

``fp32_ctx``（inc/fp32.h）描述 3 到 INT32_MAX 的奇素数域。
``fp32_init`` 校验公开模数；失败返回 -1，原上下文不变。
``fp32_reduce_raw`` 等原始取模工具还允许正的合数/偶数模数，但并不因此
将这些对象当作域。一般规范化接口使用整数余数，是可移植参考后端；
它们未承诺所有处理器上的常数时间执行，也未完成侧信道审计。

Montgomery 常量统一为正逆元 qinv = q^(-1) mod 2^32，使用减法形式。
``fp32_mont_reduce`` 是热路径内联内核，要求 |a| <= 2^31*q。
函数的输入范围、规范表示与延迟归约范围详见头文件。

``ring_ctx``（inc/ring.h）描述 F_q[X]/(X^n+1)，包含域的值副本；
不包含 NTT 表，也不依赖 NTT 是否存在。初始化允许任意正维度。

``ntt_plan``（inc/ntt.h）描述具体变换后端。本轮显式支持：

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

数学层不再主动取得随机数。旧多项式采样声明迁入 inc/poly_sample.h，
实现迁入 core；新增 poly_rand_with_rng 和 poly_uniform_CBD_with_rng。
确定性伪随机辅助函数只位于 Test/，不提供生产环境全局切换。

测试边界与后续工作
------------------

GF(16) 算术现位于 ``inc/gf16.h`` / ``src/gf16.c``，采用
``GF(2)[x]/(x^4+x+1)`` 的低四位多项式表示。矩形矩阵乘法、异或累加及
1 到 5 阶行列式位于 ``inc/gf16_matrix.h`` / ``src/gf16_matrix.c``。
SNOVA 的算法实现、签名 API、对称层及平台变体位于 ``alg/snova``；CMake
默认目标 ``pqc_snova`` 构建 ``SNOVA_OPT=10`` 优化实现，并链接 ``pqc``。
参考实现和默认优化实现均复用这些标量表示算子；扩展表示及 SIMD 热路径
仍留在 SNOVA 后端。GF(16) 求解器和其他多变量专用运算尚未拆分。

test_fp32 使用独立取模/Euclid 参考和极值检查。
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
