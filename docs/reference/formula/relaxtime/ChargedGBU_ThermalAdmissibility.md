# Charged GBU：热处方的谱解释与连续核验收

## 范围与当前决定

本文件是 [PNJL 稳定性审查](ChargedGBU_PNJLStability.md) 的后续 required follow-up，
不选择新的 regulator，也不更改 [研究方法 v1](ChargedGBU_ResearchMethod_v1.md)
或生产默认。固定 quark-only BQS 背景、质量、耦合、通道与外频率坐标。
执行证据见 [strict-audit 第22.18节](../../../dev/active/2026-08-30_charged-rpa-bu-strict-audit.md)。

**判定对象必须分开**：

1. 完整物理算符/共轭算符的热关联函数应满足的谱符号；
2. 已减法、近似或辅助场关联函数的可解释范围；
3. GBU 相移导数给出的有符号连续产额；
4. UHP 零点、Mott/Levinson 与积分收敛。

第1项不通过不等于已经找到 UHP 不稳定极点，也不由此推出介子总密度一定为负。
反之，第4项的一部分数值通过不能自动恢复第1项的物理解释。不得仅给产额改名
为“减法关联贡献”就继续宣称可与实验粒子数直接比较；改变观察量需要作者评审。

## 1. 文献事实

- [Pereira 等，Phys. Rev. C 109, 025206](https://doi.org/10.1103/PhysRevC.109.025206)，
  Sec.II C，第9页：讨论只截 A 的真空项、A/B 热项混用等不同策略；本文选择
  把 regulator 作为模型的一部分，施加于发散和收敛的全部项。其后引入两线球交集。
  因此本项目 thermal-extended 不能称为该文原样处方。详见
  [Pereira 审查](Pereira2024_RegularizationAudit.md)。
- [Laine 与 Vuorinen，Basics of thermal field theory](https://arxiv.org/abs/1701.01554)，
  所保留版本 Eq.(8.2)--(8.21)，印刷页115--117：
  Wightman 能量本征态展开、KMS 关系、retarded 上半平面解析性及谱表示。
  其 retarded 定义为 i theta(t) times commutator，且 rho=Im PiR，
  与此处使用的 Cauchy 分母 v-z 和谱符号匹配。
  来源页并非针对本项目 charged PNJL-BQS 模型，以下带化学势应用是项目推导。

文本公式已核对，PDF 原件及哈希随诊断 manifest 保留；本轮页面渲染后的图像工具
不支持读取，未声称完成版式/图像视觉验收。

## 2. 数学推导：KMS 频率与必要谱条件

对具有正热密度矩阵的正常平衡体系，以
H_gc=H-sum(mu_i N_i) 生成时间演化。对于同一算符 O_q 及其共轭，
Wightman 谱具有形式

```math
W^>(k_0,q)=\frac{2\pi}{Z}\sum_{m,n}e^{-\beta\kappa_m}
|\langle m|O_q|n\rangle|^2
\delta(k_0+\kappa_m-\kappa_n)\ge0 .
```

KMS 给出 W^<=exp(-beta*k0) W^>，于是

```math
\rho(k_0,q)=\operatorname{Im}G^R(k_0,q)
=\frac{1-e^{-\beta k_0}}{2}W^>(k_0,q),\qquad
k_0\rho(k_0,q)\ge0 .
```

正性是分布/正测度意义的；简单实极点另外按 delta 分布计入。
这里 k0 是 grand-canonical 频率，本项目
lambda=k0+mu1-mu2。不能在该不等式中直接把 k0 换成 lambda。
若使用 H 演化的物理能量，必须同步平移 KMS/Bose 参数。
这不是上游最低自由能分支的选择，也不是要求 GBU 的 d(delta)/dk0 为正。

上述是**完整物理关联函数解释的必要条件**，而非任意减法函数的公理。
PNJL 的实 Phi/PhiBar 平均场占据不自动构造正 Hilbert-space 热迹；
其推广范围须单独说明。本轮另保留普通 Fermi 控制，避免把失败全归给 PNJL。
上一阶段的“谱符号 + F(0)>0 可排除 UHP 零点”仍只是充分条件，两个命题不混同。

## 3. 不依赖相位的 q=0 反例

取正 lambda pair，公共壳上动量 p，Ei=sqrt(p^2+mi^2)，lambda=E1+E2。
在当前真空/热分别截断的同域两线核中，直接 Jacobian 给出

```math
\operatorname{Im}B_0^R(\lambda,0)
=\frac{2\pi p}{\lambda}
\left[\mathbf1_{p<\Lambda}
-\mathbf1_{p<L_{\rm th}}(n_1(E_1)+\bar n_2(E_2))\right].
```

当 Lambda<p<L_th 且 T>0、占据非零时，方括号严格为负。
这是开放频段，不是硬截断端点的赋值问题。在这一 timelike pair 区，

```math
\operatorname{Im}\Pi_P^R=
\frac{N_c}{8\pi^2}\big[\lambda^2-(m_1-m_2)^2\big]\operatorname{Im}B_0^R<0 .
```

若 k0>0，则不满足上一节的完整关联函数谱条件。有限 q 下，真空 pair 支撑
之外但热 pair 支撑以内也可构造同类证据；相位、unwrap、插值均未参与。
并且对当前实正耦合

```math
D^R=\frac{2K_a}{1-4K_a\Pi^R},\qquad
\operatorname{Im}D^R=
\frac{8K_a^2\operatorname{Im}\Pi^R}{|1-4K_a\Pi^R|^2}.
```

所以实 contact/subtraction 的改变不能纠正连续 cut 的符号。它可以移动实部和
极点，但不把该负谱变成正谱。不能用高能相位 anchor、fold 或改符号掩盖。

这个结论只针对**当前这一真空谱截断加未配对热尾的组合及其物理解释**，
不宣称所有热项不截断方法都不可用。若将其解释为低能有效/已减法响应，需要
补充适用域、匹配及高能补全，并重推与 GBU 观察量的关系；不能直接删除负尾。

## 4. PNJL 占据代数控制

当 Phi/PhiBar 非负，PNJL 占据可写为 f=<N>/3，N=0,1,2,3，
权重正比于 (1,3Phi exp(-x),3PhiBar exp(-2x),exp(-3x))。因此

```math
\frac{df}{dx}=-\frac{\operatorname{Var}(N)}{3}\le0,\qquad
f(x;\Phi,\bar\Phi)+f(-x;\bar\Phi,\Phi)=1.
```

同硬域 pair 权重 1-f(x;Phi,PhiBar)-f(y;PhiBar,Phi) 与 x+y 同号；
其中 x+y=k0/T。只延伸负占据修正而不延伸“1”破坏这个配对关系。
这里只证明占据与 pair 的代数，不冒充完整 PNJL Ward 恒等式或所有态的稳定性。

## 5. 连续核与近轴的数值检查

`causal_gbu_thermal_admissibility.jl` 提供独立 q0 壳上表达式、
Pereira (E,epsilon) 原 cut 和不经过谱插值的 Cauchy 求积。
积分使用

```math
\Pi(z)=\frac1\pi\int_a^b
\frac{\rho(v)-\rho(c)}{v-z}\,dv
+\frac{\rho(c)}{\pi}\,[\log(b-z)-\log(a-z)],\qquad c=\operatorname{Re}z .
```

c 在支撑内部时作此减法；外部取减法常数零。实轴采用 PV 对数和 +i*pi，
近轴按 eta 尺度添加积分分区；探针不得位于硬跳变端点。
分区继承当前几何支撑及 cutoff 断点提示，但被积谱值来自独立参考 cut，
不是插值节点。该分区不是全几何支撑独立证明。

自适应求积比较 n/2n 阶 Gauss 积分，以所有叶区间误差估计之和控制全局 atol；
优先细分估计误差最大的可细分区间。达到深度限制且总估计仍超标时返回失败。
估计差不是严格误差界；同时比较两套 cut/求积阶数与目标误差。
合成测试覆盖解析线性谱、内部硬跳变、近轴值、局部窄区间及预算耗尽。

独立原 cut 与现核虚部可在远离不连续点处逐点比较。PV 与有限 eta 在此只作为
各自解析函数的求值探针，不要求两者强行逐点相等。非线性 GBU 仍须另做分布/
积分极限，不能直接沿用本节的线性 Cauchy 结果。

即使样本处误差满足 4K*abs(Pi_interp-Pi_direct)<abs(F_interp)，也不代表
围道上的 uniform bound。Rouche 计数转移需要整条围道的统一上界与边界裕量；
本轮 `continuous_UHP_count_certified=false`，不能称近轴条带已完全清除。

## 6. 运行入口和作者决策点

`audit_causal_gbu_thermal_admissibility.jl` 读取保留 T170/muB240 BQS 背景，
不调用 solver。`GBU_THERMAL_OUTPUT` 必须指向不存在的 diagnostic 目录。
设置在 `settings()` 中固定且记录在 manifest：
四通道、Fermi/PNJL、四 q 的负尾探针；pi+/K+、q0/3.2、两热处方的 PV/近轴
探针；谱128/256/512与直接cut32/64及Gauss8/12阶。
输出 `spectral_witnesses.csv`、`continuous_probes.csv`，异常另存 failures；
源码、参数、冻结输入、文献及输出哈希分别绑定。

当前停止条件是完整物理介子传播子/部分粒子产额解释需要评审；不是生产
配置已被程序修改。此解释未获支持前，不继续扩大 Mott/GBU 生产验收。
可提交作者比较的候选为：

| 候选 | 已知优点 | 仍需解决 |
|---|---|---|
| Pereira 同域同硬截断 | 文献处方明确，消除本节特定的未配对热尾 | 与现有热延伸平衡背景不是同一驻点；需明确仅作冻结背景探针，或另行授权一致模型背景 |
| 重新推导带匹配/减法的热延伸 | 可探索保留高温热行为 | 需要完整公式、谱解释、适用域和 GBU 观察量推导，不能只改 contact |
| 保留当前函数作有符号诊断 | 现有数据/比较可追溯 | 不自动具有完整粒子谱或实验产额资格 |

这些是待评审选项，不是本轮自动择模。当前正/负密度或 ratio 趋势都不用于
选择 regulator。旧默认、PNJLCore、baseline、PR310 合并权限均保持不变。

作者随后授权优先研究B，见[方向B静态基础与匹配](ChargedGBU_DirectionB.md)。
本节对完整物理谱的限制不否定上游热延伸势；regulator辅助响应若要通过
GBU与产额连接，需要审查其推导及可能的额外围道/减法项，不可静默改名后采用。
