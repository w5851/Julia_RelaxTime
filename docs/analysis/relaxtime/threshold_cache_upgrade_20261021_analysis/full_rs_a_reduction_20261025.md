# 完整 RS 单线 A 的解析角积分降维

任务：`transport-full-rs-a-reduction`。本轮已完成可选内部算法入口、数值稳定性验证和单线内核计时；生产默认仍为二维 `A_aniso`。
作者于 2026-10-25 要求评估完整 A^xi 口径的解析降维；本轮只优化算法并验证，不重做完整输运计算，
也不把降维入口晋升为正式数据口径。

## 公式与适用域

现有 `A_aniso` 的定义为

\[
A^\xi=2\int_0^{P}dp\,\frac{p^2}{\sqrt{p^2+M^2}}
\int_{-1}^{1}du\,[f^+(E_\xi)+f^-(E_\xi)]-4V_\Lambda(M),
\qquad E_\xi=\sqrt{M^2+p^2(1+\xi u^2)}.
\]

其中 \(u=\cos\theta\)，\(\xi>-1\)，生产候选的热项取 \(P=\infty\)。
真空项 \(V_\Lambda\) 继续使用现有函数及其固定 \(\Lambda\)；
原分母是普通在壳能量，不能改为分布函数中的 \(E_\xi\)。

令 \(q=p\sqrt{1+\xi u^2}\)，则

\[
\frac{p^2dp}{\sqrt{p^2+M^2}}
=\frac{q^2dq}{(1+\xi u^2)\sqrt{q^2+M^2(1+\xi u^2)}}.
\]

在半无穷热域上，分布函数只依赖 \(E_q=\sqrt{q^2+M^2}\)，从而

\[
A^\xi_{\rm th}=2\int_0^\infty dq\,
W(q,M,\xi)[f^+(E_q)+f^-(E_q)],
\]

\[
W(q,M,\xi)=
\begin{cases}
\displaystyle\frac{2q}{\sqrt{\xi}}\arctan
\frac{q\sqrt{\xi}}{\sqrt{q^2+M^2(1+\xi)}},&\xi>0,\\[6pt]
\displaystyle\frac{2q^2}{\sqrt{q^2+M^2}},&\xi=0,\\[6pt]
\displaystyle\frac{2q}{\sqrt{-\xi}}\operatorname{artanh}
\frac{q\sqrt{-\xi}}{\sqrt{q^2+M^2(1+\xi)}},&-1<\xi<0.
\end{cases}
\]

负 xi 分支在实现中使用等价形式

\[
W=\frac{2q}{\sqrt{-\xi}}\operatorname{arsinh}
\frac{q\sqrt{-\xi}}{\sqrt{q^2+M^2}\sqrt{1+\xi}},
\]

避免 xi 接近 -1 时 artanh 的浮点自变量被舍入为 1。
q=0 时直接采用 W=0；小自变量使用 atan(z)/z 或 asinh(z)/z 的连续极限。
这保留完整 RS 分布，并非对 xi 的线性展开。由于普通在壳分母保留，
A 不服从内部数密度所用的简单 \(1/\sqrt{1+\xi}\) 缩放。

## 有限热截止

原球形限制 \(p\le P\) 变为 \(q\le P\sqrt{1+\xi u^2}\)。
原型显式分段处理，不把 P 直接当作变换后的统一 q 上限。

- xi>0：q≤P 时角域完整；P<q≤P√(1+xi) 时仅保留
  \(|u|\ge\sqrt{((q/P)^2-1)/\xi}\)。
- -1<xi<0：q≤P√(1+xi) 时角域完整；P√(1+xi)<q≤P 时仅保留
  \(|u|\le\sqrt{(1-(q/P)^2)/(-\xi)}\)。
- xi=0：回到原径向区间 [0,P]。

部分角域用上述角原函数在相应端点求差；径向分段的平方根端点另作平方变量代换。
这里只改变积分变量，不改变真空截止、背景热积分上限或散射的有限 s 支持。

## 验证范围与来源

[独立原型与计时脚本](../../../../benchmark/relaxtime/bench_full_rs_a_reduction.jl)
检查 120 个角权重极限组合、4 个非法 xi 输入、27 组半无穷积分及 45 组有限球形截止积分；
当前 benchmark 的等价性循环同时调用 `A_aniso_reduced` API 与独立参考表达式。
角权重用 256-bit 另一等价表达式复核；二维/一维积分各自升阶，比较热项以避免真空项抵消掩盖误差。
积分比较阈值固定为 rtol=1e-9、atol=1e-11；它不替代完整输运的 1% 验收。
8 组内核计时同时记录现有 A 和一维候选相对高阶参考的误差。

源码来自当前工作树；实现入口为 `OneLoopIntegralsCorrection.A_aniso_reduced`，
`AFieldBuilder.build_A_triplet(...; a_method=:reduced)` 与 workflow 的 `a_builder_config` 只在显式选择时启用。
[v1 来源清单](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_a_reduction_probe_v1/package_manifest.json)
对应脚本头文档字符串的加载失败，尚未进行数值评估。
[v2 来源清单](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_a_reduction_probe_v2/package_manifest.json)
仅修正该加载问题，公式与容差不变。

## 原型证据

W 上 Julia 1.12.5、单线程的 v2 运行于 2026-10-25 完成，退出码为 0。
[结果摘要](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_a_reduction_probe_v2/returned/results/summary.json)
记录 120 个角权重组合和 72 组积分比较全部通过，4 个非法 xi 输入均被拒绝。
二维/一维热项最大相对差为 `3.9968028886505635e-15`，两种积分各自升阶也通过原设门槛。
质量为零的检查针对热项与角权重；没有据此修改现有真空函数。

[8 组计时](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_a_reduction_probe_v2/returned/results/timings.json)
在 p32/cos16 的现有 A 与 q32 的一维原型之间比较，内核加速中位数为 8.397。
现有单次 A 中位耗时为 25.14–33.02 微秒。这是固定输入、预热后的 A 内核计时，
没有测量整个输运流程的加速，也没有包括加载/JIT、背景求解和通道缓存构建。
源码包的 1598 个文件在远端运行前后身份一致；返回压缩包、来源清单和结果均已核验。
后续[真实默认入口计时](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_baseline_acceptance_setup_v1/returned/reviewed_result.json)
在 W 上得到两个代表点的稳态中位时间 143.087 秒、147.331 秒，五次实际运行均复现 v3 基线。
这些运行继续使用二维 A，尚不构成一维 A 接入后完整输运收益的比较。

当前 `TransportSamplingProfiles.point_state` 在通道循环之前构造并传递 A 与有效耦合，
随后 17 道缓存消费该点状态。现有证据不支持把 A 视为每次截面求值都要重算的热点。

## 本轮实现与验证（2026-11-06）

- `A_aniso_reduced` 已进入 `OneLoopIntegralsCorrection`，支持半无穷热域和原球形有限 `p_max`；
  负 `ξ` 采用 `asinh` 稳定形式，输入显式拒绝 `ξ≤-1`、非有限参数、非正温度和非法截止。
- `AFieldBuilder` 增加 `a_method=:quadrature|:reduced`、`reduced_p_max`、`reduced_radial_nodes` 和
  `reduced_scale`。默认配置和 `TransportWorkflow` 的默认值保持 `:quadrature`，并将实际选择写入
  workflow reproducibility metadata。
- 聚焦单元测试：`tests/unit/runtests.jl` 选择
  `models/test_transport_workflow.jl` 与 `relaxtime/test_afieldbuilder.jl`，`116/116` 通过（其中 workflow 子集 `78/78`）。
- benchmark：Julia 1.12.5、单线程；120 个角权重边界组合、72 组积分对照和 8 个计时案例全部通过，
  最大半无穷热项相对差 `3.9968028886505635e-15`；重复运行的单线 A 内核中位加速约 `8.6x`
  （本轮 `8.56x`，前一轮 `8.89x`，计时受运行时波动影响）。
  计时只覆盖固定背景下的 A 内核，不代表完整输运流程加速；本轮没有改动正式输运结果。

本项的工程目标至此完成：算法已可选、默认口径未变、数值边界和性能证据已覆盖。完整输运入口的收益测量、
更广泛扫描和是否把 `:reduced` 晋升为默认，仍应作为后续独立任务处理，不阻塞本项归档。
