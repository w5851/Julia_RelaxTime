# 介子数密度与同位旋不对称 PNJL/NJL：计算流程定向文献调研

> **状态**：研究背景调研（diagnostic background survey），不是数值结果、回归基线或 production 授权。
> **日期**：2026-08-28
> **工作区**：codex/charged-kpi-conserved-chemicals 隔离 worktree
> **研究对象**：当前项目关注的 \(K^+/\pi^+\)、\(K^-/\pi^-\)、\(2+1\) 味 PNJL/NJL、同位旋不对称、完整 \(\mathcal K_{ab}\) 相互作用、RPA/相移/Beth–Uhlenbeck（BU）介子密度以及化学冻结线。

本文记录与当前问题最相关的已发表论文和公开预印本中所采用的**计算流程**，并将它们映射到本项目现有实现。文献中的符号、\(\mu_I\) 归一化、KMT 耦合常数符号可能不同；下文保留各文献自己的约定，并在需要处显式标注，不能把不同约定的数值直接互换。

## 1. 先给结论

### 1.1 “完整同位旋交叉项”是一个明确的模型层问题

师兄所说的同位旋交叉项，与附件中出现的

\[
\mathcal K_{30}^{\pm},\quad \mathcal K_{03}^{\pm},\quad
\mathcal K_{38}^{\pm},\quad \mathcal K_{83}^{\pm}
\]

相符。它们来自三味 KMT 六费米相互作用在收缩一个夸克对后形成的有效四费米耦合矩阵，而不是单纯把 \(\mu_u\) 和 \(\mu_d\) 设成不同。Rehberg、Klevansky 和 Hüfner 的 SU(3) NJL 推导直接给出了这些项；在 \(m_u=m_d\)、\(\mu_u=\mu_d\) 的同位旋对称极限，\(\mathcal K_{03}\) 与 \(\mathcal K_{38}\) 消失，\(\pi^0\) 与 \(0/8\) 中性扇区解耦。因而“只保留 \(0/8\) 混合”是一个有明确适用条件的简化，而不是完整三味拉氏量在任意同位旋不对称下的结果。

### 1.2 文献中至少有三条不同的计算链路

1. **比值/介子谱的现象学链路**：对角夸克平均场（MFA）
   \(\rightarrow\) 带电介子 RPA 极点或相移
   \(\rightarrow\) 稳定粒子、BW 或 BU 密度
   \(\rightarrow\) 沿手征交叉线或化学冻结线画 \(K/\pi\)。这正是 Friesen 2019、Blaschke 2020 系列工作的主要路线。
2. **完整有效耦合矩阵链路**：从 KMT 相互作用得到所有 \(\mathcal K_{ab}\)
   \(\rightarrow\) 求解含全部允许凝聚的 MFA
   \(\rightarrow\) 用完整 RPA 矩阵求中性介子混合。Rehberg 1996 给出一般耦合公式；Mei 2023 在外磁场导致 \(\sigma_u\neq\sigma_d\) 时使用了 \(0/3/8\) 中性混合。
3. **显式同位旋化学势链路**：以 \(\mu_B,\mu_I,\mu_S\) 为输入，允许 \(\pi\)（必要时 \(K\)）凝聚，将凝聚作为 MFA 变量，再以扩展的 RPA 矩阵求谱。Xia 2013、Lopes 2021 代表三味模型，Wu 2017 代表两味带电/中性 pion 的明确同位旋破缺处理。

这三条链路可以组合，但文献通常不会一次性完成“完整 \(\mathcal K_{ab}\) + BQS 约束 + 全部 BU 介子热力学反馈 + 冻结线生产”。因此当前项目应把它们视为正交的模型选择，而不是一个开关。

本轮补充的 Tian 等工作（Phys. Rev. D 114, 034012 (2026)，DOI：[10.1103/d7nm-y2vp](https://doi.org/10.1103/d7nm-y2vp)）提供了可直接复核的 Eq. (2)-(3) 有效四费米耦合展开，以及 Eq. (20)-(22)、(26) 的中性 `(0,3,8)` RPA 矩阵组织。该文采用 NJL、外磁场、Ritus 基底和 Pauli--Villars 正则化，化学势取味道相等，目标是极点质量而非 BU 密度或冻结线；因此在本项目中作为“公式/矩阵结构信源”，不作为 PNJL/BU 数值基线。

### 1.3 \(\mu_s=0.55\mu_u\) 是现象学消融设定，不是普适平衡条件

Friesen 2019 和 Blaschke 2020 的相关工作确实使用了 \(\mu_u=\mu_d\)、\(\mu_s=0.55\mu_u\) 的案例；同时它们明确把这类 Case I/Case II 视作模型情景或现象学输入。它不是由 \(\rho_Q/\rho_B=0.4\)、\(\rho_S=0\) 自动推出的平衡解，也不能替代固定 \((T,\mu_B)\) 后求解 \(\mu_Q,\mu_S\) 的守恒荷路线。

对于本项目，\(\mu_s=0.2\mu_u\) 和 \(0.55\mu_u\) 最适合作为**同点消融**，用于回答“上游奇异化学势设定改变了多少”，而不是作为完整同位旋不对称方程组的最终定义。

### 1.4 完整 \(\mathcal K_{ab}\) 与介子反馈属于两个层次

完整 \(\mathcal K_{ab}\) 首先改变的是：

- 中性 \(\pi^0-\eta-\eta'\) 混合及其谱函数；
- 由 \(\sigma_u-\sigma_d\) 产生的耦合矩阵非对角元；
- 通过自洽 MFA/RPA 重新得到的夸克质量、极化函数、带电介子耦合和相移。

它**不等价于**把 \(\rho_s=0\) 改写为 \(\rho_s+\rho_K=0\)，也不自动要求将 \(\Omega_M\) 放入平均场驻点方程。后者是独立的热力学闭合问题。当前分支中的 partial feedback 只把 \(\pi^\pm,K^\pm\) 的守恒荷贡献放进外层约束，并未把介子贡献加入 Omega PNJL 的 stationarity；这一点与文献中常见的“夸克平衡解后处理介子谱/密度”路线更接近，而不是完整 \(\Omega_q+\Omega_M\) 自洽求解。

### 1.5 当前 \(K^\pm/\pi^\pm\) 接近 \(O(1)\) 不能单独归因于缺少 \(K_{30}\) 等项

在带电通道中，化学势差、介子质量、相移/谱权重、简并度、积分窗口和 \(x_{\min}\) 截断都会共同影响密度。完整 \(K_{30}\)、\(K_{38}\) 的最直接效应位于中性混合扇区；它们对 \(K^\pm/\pi^\pm\) 的影响通常是通过重新求解凝聚、质量和耦合的**间接效应**。因此“加入完整拉氏量后比值一定变成 \(0.3\)”或“比值接近 1 一定是 \(\mu_s\) 被拉回来了”都超出了目前证据。应先做同一数值设定下的 full-\(\mathcal K_{ab}\)/diagonal-\(\mathcal K_{ab}\) 对照，再讨论归因。

## 2. 研究问题、检索和筛选方法

### 2.1 研究问题

本次调研只回答以下四类问题：

1. 已发表 PNJL/NJL 工作如何选择 \(\mu_u,\mu_d,\mu_s\)，尤其如何区分 \(\mu_s=0.55\mu_u\)、显式 \(\mu_I\) 和 BQS 守恒荷约束？
2. 完整 KMT 有效四费米耦合中，\(\mathcal K_{30},\mathcal K_{03},\mathcal K_{38},\mathcal K_{83}\) 在什么条件下出现、何时消失，进入哪一个 RPA 扇区？
3. \(K/\pi\) 介子密度采用稳定粒子、BW 极点还是 BU 相移？介子密度是否反馈回平衡态？
4. 计算结果如何沿手征线、伪临界线或化学冻结线组织？这些路线和当前项目的 \(K^+/\pi^+\)、\(K^-/\pi^-\) 目标如何对应？

不把本次调研当作穷尽所有 PNJL 文献的系统综述，也不把文献中的现象学参数化当作本项目的物理标准答案。

### 2.2 数据库、时间和检索日志

- **检索日期**：2026-08-28；优先公开全文和 DOI/期刊元数据。
- **主要来源**：arXiv HTML/full text、APS/DOI 页面、SCOAP3 开放全文片段、IHEP/CPC 开放页面。
- **检索策略**：先以模型和问题词定位记录，再回到全文提取“化学势—MFA—RPA/相移—密度—路径—反馈”链路。
- 动态搜索页面的总结果数随索引和过滤器变化，未把不稳定的网页计数伪装成精确统计；以下记录保留可复查的查询族、纳入记录数和排除原因。

| 查询族（2026-08-28） | 目的 | 处理结果 |
|---|---|---|
| site:arxiv.org PNJL isospin asymmetric meson K30 K03 lambda 0 3 8 | 找完整 \(0/3/8\) 耦合与中性混合 | 纳入 Rehberg、Mei；补充 Xia |
| site:arxiv.org Beth Uhlenbeck kaon pion ratio 2+1 flavor PNJL Blaschke | 找 BU \(K/\pi\) 密度 | 纳入 Blaschke 2020 BU |
| site:arxiv.org entanglement SU(3) PNJL kaon to pion ratio Blaschke | 找 \(\mu_s=0.55\mu_u\)、冻结线/伪临界线 | 纳入 Blaschke 2020 EPNJL |
| PNJL isospin chemical potential meson properties neutral pi eta mixing | 找显式 \(\mu_I\) 和 RPA 矩阵 | 纳入 Xia、Lopes、Wu |
| K30 NJL meson | 追踪 KMT 有效耦合记号 | 纳入 Rehberg、Mei |
| K03 sigma_u sigma_d NJL | 核对 \(\sigma_u-\sigma_d\) 依赖 | 纳入 Mei；与 Rehberg 交叉核对 |
| Three-flavor Nambu-Jona-Lasinio model at finite isospin chemical potential arXiv | 三味 \(\mu_I\) MFA/RPA 流程 | 纳入 Xia |
| Hot QCD at finite isospin density SU(3) NJL lattice data | 三味有限同位旋热力学路线 | 纳入 Lopes |
| Rehberg Klevansky Hüfner 1996 SU(3) NJL K30 K38 effective couplings | 获取早期完整耦合矩阵基准 | 纳入 Rehberg |

### 2.3 纳入、去重和排除

**纳入条件**：至少满足以下一项，并能从全文或方法摘要重建流程：

- 三味/2+1 味 NJL/PNJL 的 MFA 与介子 RPA；
- 显式同位旋化学势或 \(u/d\) 不对称导致的中性混合；
- KMT 有效耦合矩阵中的 \(0,3,8\) 非对角项；
- \(K/\pi\) 介子密度、相移/BU 或冻结线扫描。

**去重规则**：同一工作在 arXiv、期刊和会议页面的版本按一个记录计；例如 Rehberg 的预印本与 PRC 53 论文、Blaschke 的 arXiv 与期刊页面不重复计数。

**排除规则**：

- 只有实验数据、没有模型计算流程的记录不作为方法证据；
- 综述、重复会议摘要和与介子/同位旋无关的 PNJL 工作不进入核心表；
- 纯磁场研究若没有明确的 \(\mathcal K_{ab}\) 或混合矩阵构造，不作为完整交叉项证据；Mei 2023 保留的原因是它确实给出了完整 \(0/3/8\) 方法，磁场只是在其算例中产生 \(\sigma_u\neq\sigma_d\)。

问题核心记录共 8 项：Rehberg 1996、Xia 2013、Friesen 2019、Blaschke 2020（两篇工作）、Lopes 2021、Mei 2023、Wu 2017。另列 2 篇专门说明 BU 介子热力学反馈的 supporting literature。它们的目标不同，因此下面按“流程”而不是按年份评价。

## 3. 统一符号和当前项目的对照基线

### 3.1 三种化学势语言不能混用

本项目的 BQS 映射（见 src/models/solver/spec/ConstraintModes.jl）是

\[
\mu_u=\frac{\mu_B}{3}+\frac{2\mu_Q}{3},\qquad
\mu_d=\frac{\mu_B}{3}-\frac{\mu_Q}{3},\qquad
\mu_s=\frac{\mu_B}{3}-\frac{\mu_Q}{3}-\mu_S.
\]

对应的密度定义为

\[
\rho_B=\frac{\rho_u+\rho_d+\rho_s}{3},\qquad
\rho_Q=\frac{2\rho_u-\rho_d-\rho_s}{3},\qquad
\rho_S=-\rho_s.
\]

显式 \(\mu_I\) 文献有至少两种归一化：

- Xia et al. 2013 使用
  \[
  \mu_u=\mu_B/3+\mu_I/2,\quad
  \mu_d=\mu_B/3-\mu_I/2,
  \]
  因而这里 \(\mu_I=\mu_u-\mu_d\)。
- Lopes et al. 2021 使用
  \[
  \mu_u=\mu_B/3+\mu_I,\quad
  \mu_d=\mu_B/3-\mu_I,
  \]
  因而这里 \(\mu_u-\mu_d=2\mu_I\)。

所以“指定 \(\mu_I\)”必须先说明归一化；它与给定 \(\rho_Q/\rho_B\) 的约束求解也不是同一个入口。二者在小不对称、无介子凝聚的局部线性近似下可以建立响应关系，但不能直接把一个参数替换成另一个参数。

### 3.2 附件中的 \(\mathcal K_{30}\) 等项来自哪里

Rehberg et al. 在收缩 KMT 六费米项后给出的有效耦合满足（这里采用其记号，整体正负号依赖标量/赝标量通道约定）

\[
\begin{aligned}
K_{03}^{\pm}=K_{30}^{\pm}
&=\mp\frac{1}{2\sqrt 6}N_cK
   \left(i\operatorname{tr}_\gamma S^u-i\operatorname{tr}_\gamma S^d\right),\\
K_{38}^{\pm}=K_{83}^{\pm}
&=\pm\frac{1}{2\sqrt 3}N_cK
   \left(i\operatorname{tr}_\gamma S^u-i\operatorname{tr}_\gamma S^d\right),\\
K_{08}^{\pm}=K_{80}^{\pm}
&=\pm\frac{\sqrt2}{12}N_cK
   \left(i\operatorname{tr}_\gamma S^u+i\operatorname{tr}_\gamma S^d-2i\operatorname{tr}_\gamma S^s\right).
\end{aligned}
\]

由于 \(i\operatorname{tr}_\gamma S^f\) 与味凝聚等价，后两组可理解为

- \(K_{03},K_{38}\propto \sigma_u-\sigma_d\)：同位旋破缺的直接诊断；
- \(K_{08}\propto \sigma_u+\sigma_d-2\sigma_s\)：即使 \(u/d\) 对称也通常不为零，负责 \(0/8\) 混合。

在 \(\sigma_u=\sigma_d\) 且没有 pion condensate 的极限，\(K_{03}=K_{30}=K_{38}=K_{83}=0\)，中性矩阵退化为 \(0/8\) 块。Mei et al. 2023 用等价的 \(\sigma_u-\sigma_d\) 结构构造了中性 \(\pi^0-\eta-\eta'\) 的 \(3\times3\) 混合矩阵；其符号与 Rehberg 的凝聚定义可能相反，不能逐项比较正负号而不先统一约定。

### 3.3 当前代码已经实现和没有实现的部分

| 当前入口 | 已有能力 | 与完整同位旋路线的差距 |
|---|---|---|
| src/models/pnjl_physics/PNJLCore.jl | 对角 \(\phi_u,\phi_d,\phi_s\) 质量方程；\(\Omega\) 含凝聚、真空、热项和 Polyakov 项 | 没有非对角夸克自能或 pion/kaon 凝聚变量 |
| src/models/omega.jl | 只组装 PNJL 夸克/平均场巨势 | 没有 \(\Omega_M\) 介子热力学项 |
| src/models/solver/spec/ConstraintModes.jl | BQS 与 flavor 化学势映射及 \(\rho_Q/\rho_B,\rho_S\) 约束 | 约束是夸克层闭合；尚未把介子、反介子荷作为完整热力学自由度 |
| src/relaxtime/EffectiveCouplings.jl | \(K_0,K_{123},K_{4567},K_8,K_{08}\)；\(u/d\) 使用合并的 \(G^\mu\) | 没有独立 \(G^u,G^d\)，也没有 \(K_{03},K_{30},K_{38},K_{83}\) |
| src/relaxtime/MesonPropagator.jl | \(0/8\) 两维混合矩阵 | 没有 \(0/3/8\) 三维中性矩阵 |
| src/relaxtime/PolarizationAniso.jl | 带电极化函数使用 \(\lambda=k_0+\mu_1-\mu_2\) | 带电通道虽可使用 \(\mu_u\neq\mu_d\)，但耦合矩阵仍是简化版本 |
| src/relaxtime/MesonDensity.jl | \(\pi^\pm,K^\pm\) 组分和 BU 密度；当前诊断为 x_min_cut | 目前没有 \(K^0,\bar K^0\)、\(\pi^0,\eta,\eta'\) 密度的完整生产链 |
| scripts/analysis/relaxtime/meson_conserved_charge_feedback_runtime.jl | 外层 partial feedback：重算夸克 gap，再加四个带电介子守恒荷 | 不是 \(\Omega_q+\Omega_M\) 的驻点求解，也不是完整 hadronic strangeness constraint |

因此，当前项目最接近下面第一个文献族的“带电介子后处理”版本，并正在试探其外层守恒荷修正；它还不是第二、第三族的完整拉氏量/凝聚实现。

## 4. 核心文献的计算流程对照

下表中的“证据强度”含义为：**直接全文公式**（可逐式核对）、**全文方法段**（流程明确但细节分散）或**摘要/项目映射**（只能支持较弱结论）。

| 文献 | 化学势/约束 | MFA 与介子谱 | 密度/路径 | 介子反馈 | 与本项目的关系 |
|---|---|---|---|---|---|
| Rehberg, Klevansky & Hüfner, PRC 53 (1996) | 给出任意 \(\mu_f\) 的夸克传播子公式；实际数值主要在 \(\mu_f=0\) | Hartree + RPA；显式给出所有 \(K_{ab}\)，同位旋对称时只剩 \(0/8\) 中性混合 | 介子质量和 Mott 温度；非 BU 冻结线 | 未做 \(\Omega_M\) 反馈 | **完整 \(\mathcal K_{ab}\) 的基准推导**，但不是 PNJL/Kπ 冻结线计算 |
| Xia, He & Zhuang, PRD 88 (2013) | \(\mu_u=\mu_B/3+\mu_I/2\)、\(\mu_d=\mu_B/3-\mu_I/2\)、\(\mu_s=\mu_B/3-\mu_S\)；允许 pion condensate | 三味 NJL MFA；以 \(\sigma_f,\phi_{ud},\phi_{us},\phi_{ds}\) 为变量；扩展耦合矩阵并用 \(\det[1-2G\Pi]=0\) 求谱 | 有限 \(\mu_I\) 相图和介子质量；不做 K/π 冻结线 | 凝聚通过 MFA 驻点反馈，未做 BU 热介子 \(\Omega_M\) | **显式 \(\mu_I\)+扩展 RPA** 的主要方法参考 |
| Friesen, Kalinovsky & Toneev, PRC 99 (2019) | Case I：\(\mu_u=\mu_d=\mu_s\)；Case II：\(\mu_u=\mu_d,\ \mu_s=0.55\mu_u\) | SU(3) PNJL 对角 MFA；Bethe–Salpeter/RPA 介子极点，可出现复极点 | 稳定 Bose 密度；\(\mu_K=\mu_u-\mu_s\)，\(\mu_\pi=0.135\) GeV；沿手征转变线画比值；另给 Cleymans 型冻结线参数化 | 介子作为后处理，不把介子密度加入 MFA 巨势 | **0.55 消融、K/π 量级和冻结线参数化**的直接先例；作者也指出 Case I/II 是人为情景 |
| Blaschke et al., EPJST 229 (2020) | \(\mu_s=0.55\mu_u,\ \mu_u=\mu_d\)；另讨论守恒律给出的 x 依赖 \(\mu_s,\mu_\pi\) | EPNJL 对角 MFA；RPA 极点、散射切和相移；比较 BW 与 BU | 稳定粒子、BW 和广义 BU；\(\mu_{K^+}=\mu_u-\mu_s\)、\(\mu_{K^-}=\mu_s-\mu_u\)；沿伪临界线作为冻结线代理，并试验固定 \(\mu_\pi\) 的等温缝合 | 介子密度用于比值后处理/统计模型，不是完整 \(\Omega_M\) stationarity | **当前 K±/π± BU 路线最接近的公开先例**；其自身结论是冻结线代理与 BW/BU 选择会显著改变 horn |
| Blaschke, Friesen, Kalinovsky & Radzhabov, Particles 3 (2020) | 2+1 PNJL，参数/守恒律情景与上文工作相关 | Beth–Uhlenbeck 相移，强调介质中的 anomalous K+ mode | 讨论 kaon-to-pion 比值的温度/密度行为 | 后处理 | 作为短篇会议版背景记录；细节以 EPNJL 全文为准 |
| Lopes et al., PRD 103 (2021) | \(\mu_u=\mu_B/3+\mu_I,\ \mu_d=\mu_B/3-\mu_I,\ \mu_s=\mu_B/3-\mu_S\)；明确讨论 \(\mu_I\) | 三味 NJL MFA；包含 pion condensate \(\Delta\)，比较 TRS/MSS 正则化；重点是 EOS/热力学 | 与晶格数据比较；不是 K/π BU 冻结线 | \(\Delta\) 通过 MFA 反馈，未加入热介子 \(\Omega_M\) | **显式 \(\mu_I\) 与平衡热力学**的参考，但不直接给带电 K/π 密度 |
| Wu & collaborators, CPC 41 (2017) | 两味 \(m_u\neq m_d\)，有限 \(T,\mu_I,\mu_B\) | SU(2) NJL MFA + RPA；计算带电/中性 pion 质量和衰变常数分裂 | 不含 strange/K，也不做冻结线 | 仅夸克 MFA 自洽 | 说明即使没有三味 KMT，也可直接把 \(u/d\) 不对称带入带电 pion RPA；是两味方法支撑 |
| Mei et al., PRD 107 (2023) | 外磁场、\(T,\mu_B\) 造成 \(\sigma_u\neq\sigma_d\)，不采用本项目的 \(\rho_Q/\rho_B\) | SU(3) NJL；完整 \(K_{30},K_{03},K_{38},K_{83}\) 等；中性 \(\pi^0-\eta-\eta'\) 用 \(3\times3\) RPA 混合 | 中性介子质量谱，研究磁场/温度/密度，不做 K/π BU | 没有热介子 \(\Omega_M\) | **完整 \(0/3/8\) 实现的现代方法证据**；外场来源不同，不能直接当作重离子冻结线数值基准 |
| Blaschke et al., Ann. Phys. 348 (2014) | SU(2) NJL，有限 \(T,\mu\)，不以 BQS 冻结线为目标 | 玻色化后在 meson/diquark 场的 Gaussian 近似中构造相移 | BU 形式的热力学势；研究 Mott 解离、连续谱和 pionic pressure | **介子项进入 \(\Omega\)**，并讨论自洽准粒子图景 | supporting：说明“BU 密度/相移”可以上升为热力学闭合，但不是三味 K/π 计算 |
| Maslov & Blaschke, PRD 107 (2023) | 2 味 PNJL，研究 EOS，不采用本项目的 \(\rho_Q/\rho_B\) 冻结线约束 | 全动量依赖的介子极化；保留 timelike 与 spacelike（Landau damping）相关 | BU 介子压力、trace anomaly，并用 Feynman–Hellmann 评估 pion/Landau 对凝聚的贡献 | **介子 off-shell 贡献反馈到热力学量/凝聚分析**，但不是统一 BQS+K/π 求解 | supporting：说明 \(\Omega_M\) 反馈需要全动量谱与 cutoff 依赖审计 |

### 4.1 Rehberg 1996：完整耦合矩阵的早期基准

Rehberg et al. 的工作最适合回答“附件中的交叉耦合在理论上是否应该出现”。其流程是：

\[
\mathcal L_{\rm NJL}^{(4q+6q)}
\longrightarrow
\text{收缩一条夸克线得到 }K_{ab}^{\pm}
\longrightarrow
\text{Hartree gap equations}
\longrightarrow
\text{RPA 极化函数与介子传播子}.
\]

他们列出了

\[
K_{03}=K_{30}\propto i\operatorname{tr}S^u-i\operatorname{tr}S^d,
\qquad
K_{38}=K_{83}\propto i\operatorname{tr}S^u-i\operatorname{tr}S^d,
\]

并指出这些非对角耦合会造成 \(\pi^0-\eta-\eta'\) 混合；若恢复 SU(2) 同位旋对称，则它们为零，\(\pi^0\) 解耦。需要特别保留的限制是：论文的公式允许一般 \(\mu_f\)，但其数值结果主要取零化学势，所以它证明的是**模型结构和 RPA 组织方式**，不是本项目的冻结线数值验证。

### 4.2 Friesen 2019：\(\mu_s=0.55\mu_u\) 与 K/π 现象学路线

这篇 SU(3) PNJL 工作从对角 MFA 出发，先解夸克质量和 Polyakov 背景，再通过 Bethe–Salpeter 方程得到介子极点，最后使用

\[
n_M=d_M\int\frac{d^3q}{(2\pi)^3}
\frac{1}{\exp[(E_M(q)-\mu_M)/T]-1}
\]

计算稳定介子密度。带电通道使用 \(\mu_{K^+}=\mu_u-\mu_s\)、\(\mu_{K^-}=\mu_s-\mu_u\)；pion 处使用固定的 \(\mu_\pi=0.135\) GeV。论文给出的 Case II 正是 \(\mu_u=\mu_d\)、\(\mu_s=0.55\mu_u\)，但同时承认这些 Case 是为研究有限密度效应而构造的情景，不是对真实碰撞守恒律的唯一推导。

其“冻结线”更多是手征转变线上的比值扫描或 Cleymans 型参数化的现象学比较，并非在每个 \((T,\mu_B)\) 上同时解完整 BQS 加介子约束。因此该工作支持“0.55 是文献中的可复现消融值”，不支持把它当作本项目的最终平衡条件。

### 4.3 Blaschke 2020：RPA 极点、BW、BU 和路径选择的分层比较

EPNJL 工作把介子极化函数分解为夸克-反夸克 pair 部分和 scattering 部分，先构造 RPA 传播子，再由相移计算：

\[
n_M=d_M\int\frac{d^3q}{(2\pi)^3}
\int\frac{d\omega}{2\pi}
g_M(\omega,q)\,\frac{\partial\delta_M(\omega,q)}{\partial\omega}.
\]

它比较了：

- 只保留窄极点的稳定/BW 近似；
- 计入散射切和连续谱的广义 BU；
- 沿伪临界线把 \(x=T/\mu_B\) 作为冻结线代理；
- 固定 \(\mu_\pi\) 的等温试验与更具守恒律动机的 x 依赖 \(\mu_s,\mu_\pi\)。

论文的一个关键方法结论是：BW 和 BU、固定参数和路径相关的 \(\mu_s,\mu_\pi\) 会改变 \(K^+/\pi^+\) 的峰形；\(K^-/\pi^-\) 在其设定下更倾向于单调上升。因而“沿冻结线比较”不是单独的后处理绘图步骤，而是化学势参数化、谱函数方案和路径三者的联合选择。

### 4.4 Xia 2013 与 Lopes 2021：显式 \(\mu_I\) 的平衡态路线

Xia et al. 从

\[
\mu_u=\mu_B/3+\mu_I/2,\quad
\mu_d=\mu_B/3-\mu_I/2,\quad
\mu_s=\mu_B/3-\mu_S
\]

出发，在 MFA 中同时最小化 \(\sigma_l,\sigma_s,\phi_{ud}\)（以及可能的 \(us,ds\) 凝聚），再将有效耦合和极化函数组织成矩阵，以

\[
\det[1-2G\Pi]=0
\]

求介子质量。有限 \(\mu_I\) 不只是把带电极化函数中的 \(\mu_u-\mu_d\) 换成一个数；当 \(\phi_{ud}\neq0\) 时，标量—赝标量和更多中性/带电通道发生混合，平衡方程本身也扩大。

Lopes et al. 采用另一种 \(\mu_I\) 归一化，重点是三味 NJL 的 EOS、压力和晶格比较，并比较 TRS 与 MSS 正则化。它说明显式 \(\mu_I\) 是一条成熟的热力学路线，但其目标不是 K/π 介子密度。对本项目而言，这两篇文章最重要的启示是：如果真的要把 \(\mu_I\) 作为上游输入，就应连同可能的凝聚变量、正则化和扩展 RPA 一起定义，而不是只在下游把 \(\mu_I\) 塞进 \(\mu_M\)。

### 4.5 Mei 2023：\(0/3/8\) 中性混合的现代实现

Mei et al. 在磁场下求解三味 NJL 的 \(\sigma_f\)，随后使用包含 \(K_{30},K_{03},K_{38},K_{83}\) 的耦合矩阵，构造中性 \(\pi^0,\eta,\eta'\) 的三维传播子。其流程可以概括为：

\[
\text{外场/}T/\mu_B
\longrightarrow
\sigma_u,\sigma_d,\sigma_s
\longrightarrow
K_{ab}(\sigma_u,\sigma_d,\sigma_s)
\longrightarrow
\Pi_{ab}
\longrightarrow
3\times3\text{ 中性 RPA 谱}.
\]

该文不是本项目的冻结线/BU 先例，因为同位旋破缺来自磁场而不是 \(\rho_Q/\rho_B=0.4\)。但它证明了把当前代码的 \(0/8\) 矩阵扩展到 \(0/3/8\) 在数值上是可执行的，且应从同一个 \(K_{ab}\) 构造器同时生成耦合、极化和极点方程，不能只在下游手工添加一个混合项。

### 4.6 supporting：从相移后处理到 \(\Omega_M\) 的热力学闭合

Blaschke et al.（Annals of Physics 348, 2014）从玻色化的 NJL 出发，在介子和二夸克场的 Gaussian 近似下得到 BU 形式的热力学势，并用介质相移描述束缚态到散射连续谱的 Mott 解离。这里的关键不是某个固定 \(K/\pi\) 数值，而是把相移积分放进压力等热力学量；这与当前只计算 \(n_{\pi^\pm},n_{K^\pm}\) 并将其放入外层守恒荷残差的 partial feedback 有明确层次差异。

Maslov & Blaschke（PRD 107, 094010, 2023）进一步在 2 味 PNJL 中保留介子极化的完整动量依赖，指出 spacelike 区域的 Landau damping 会改变介子压力，并通过 Feynman–Hellmann 关系分析 pion 气体和 Landau damping 对手征凝聚的贡献。该结果提示：如果未来把 \(\Omega_M\) 真正反馈到 PNJL 平均场，不能只把一个 on-shell Bose 密度加到约束方程；至少要明确相移窗口、spacelike 连续谱、动量 cutoff 以及 \(\partial\Omega_M/\partial\sigma_f\) 的求导口径。

这两篇 supporting 文献并不提供本项目的三味 \(K^\pm/\pi^\pm\) 冻结线答案，也不意味着当前应立即实现完整反馈；它们只用于界定“完整反馈”所需的物理和数值工作量。

## 5. 按“计算流程”归类的证据地图

### 5.1 路线 A：对角 PNJL MFA + 带电介子后处理

这是当前 \(K^\pm/\pi^\pm\) 比值最直接的公开先例：

\[
\begin{aligned}
&\text{给定 }(T,\mu_u,\mu_d,\mu_s,\Phi,\bar\Phi)\\
&\quad\downarrow\\
&\text{解对角 PNJL gap equations}\ (M_u,M_d,M_s)\\
&\quad\downarrow\\
&\text{带电 RPA 极化}\ \Pi_{u\bar d},\Pi_{d\bar u},\Pi_{u\bar s},\Pi_{s\bar u}\\
&\quad\downarrow\\
&\text{极点/BW 或相移/BU}\ n_{\pi^\pm},n_{K^\pm}\\
&\quad\downarrow\\
&\text{沿手征线/冻结线输出 }K^+/\pi^+,K^-/\pi^- .
\end{aligned}
\]

优点是成本低、能直接分析带电通道；缺点是若上游只解夸克 \(\rho_S=0\)，介子携带的奇异荷并没有进入平衡条件。当前分支的 partial feedback 正是在这条链路末端加了一个外层守恒荷修正，但仍保留 quark-only \(\Omega\) 的驻点定义。

### 5.2 路线 B：完整 \(\mathcal K_{ab}\) + 中性 RPA

\[
\text{KMT}
\rightarrow K_{00},K_{33},K_{88},K_{03},K_{08},K_{38}
\rightarrow \Pi_{0,3,8}
\rightarrow \det D^{-1}_{0,3,8}=0.
\]

在同位旋对称极限，这个矩阵退化为 \(0/8\) 块；在 \(\sigma_u\neq\sigma_d\) 或显式 \(m_u\neq m_d\) 时，\(\pi^0\) 会与 \(\eta,\eta'\) 混合。此路线的验证首要是“对称极限回归”，不是先看 \(K/\pi\) 比值是否变成实验数值。

### 5.3 路线 C：显式 \(\mu_I\) + 凝聚 + 扩展 RPA

\[
\text{给定 }(\mu_B,\mu_I,\mu_S)
\rightarrow
\text{加入 }\sigma,\phi_{ud},\phi_{us},\phi_{ds}
\rightarrow
\partial\Omega_{\rm MFA}=0
\rightarrow
\text{矩阵 RPA/介子谱}.
\]

这条路线把同位旋不对称放在上游拉氏量和平均场中。如果 \(|\mu_I|\) 进入 pion 凝聚区，只在下游修改 \(\mu_{\pi^\pm}\) 会漏掉相变和混合，因此不能把它看作当前 partial feedback 的小改动。

### 5.4 路线 D：完整介子热力学反馈（尚非本项目当前实现）

严格的热力学反馈应从

\[
\Omega_{\rm tot}(T,\mu;X)
=\Omega_{\rm PNJL}(T,\mu;X_q)
+\Omega_M[T,\mu_M;\delta_M(X_q)]
\]

出发，对所有独立的夸克凝聚、Polyakov 变量和必要的介子序参量统一求驻点，并在守恒荷方程中同时包含夸克和介子贡献。当前 partial feedback 做的是：

\[
\text{quark-only baseline}
\rightarrow
\text{在候选 }(\mu_Q,\mu_S)\text{ 下重解夸克 gap}
\rightarrow
\text{计算 }\pi^\pm,K^\pm\text{ BU 密度}
\rightarrow
\text{把介子荷放进外层残差}.
\]

两者都可作为研究工具，但后者不能被称为完整 \(\Omega_M\) 自洽解；文献中 K/π 比值工作多采用路线 A，而不是路线 D。

## 6. 对当前项目问题的直接解释

### 6.1 为什么完整交叉项值得做，但不能先把它当作 ratio 的唯一原因

当前实现把 \(u/d\) 的凝聚和耦合压缩为 \(G^\mu\)，并只保留 \(0/8\) 中性矩阵。若目标是物理上自洽地研究 \(\rho_Q/\rho_B=0.4\)，至少存在两个可分离的缺口：

1. **上游模型缺口**：非零 \(\sigma_u-\sigma_d\) 时缺少 \(K_{03},K_{30},K_{38},K_{83}\) 和 \(0/3/8\) 中性矩阵；
2. **守恒荷闭合缺口**：quark-only 的 \(\rho_S=0\) 没有包含 \(K^\pm\) 等介子携带的 strange charge。

第一个缺口会改变质量、耦合和中性谱，并可能通过 condensate 反馈间接影响带电比值；第二个缺口会改变求解得到的 \(\mu_Q,\mu_S\)。两者的方向和大小不能从“\(\mu_s\) 是否变大”这一项单独判断。

### 6.2 对师兄“完整拉氏量后应直接用 quark-only 结果”的理解

这句话在一个清晰的近似层次下是合理的：

- 先把完整 \(K_{ab}\) 和 \(\mu_I\) 纳入 PNJL/NJL 的**夸克层平衡态**；
- 用该平衡态计算介子传播子、相移和 \(K^\pm/\pi^\pm\) 密度；
- 不把介子密度再反馈到同一个平衡方程，作为一个受控的“quark-equilibrium + meson postprocessing”近似。

它并不意味着介子对守恒荷在所有问题中都可以忽略。若研究问题是“给定重离子守恒荷的化学势”，介子也携带 \(Q,S\)，则应把以下两条结果并列报告：

\[
\text{(A) full-kernel quark-only equilibrium + meson postprocessing},
\]

\[
\text{(B) full-kernel quark-only seed + partial conserved-charge correction}.
\]

只有当两者差异在目标误差内，才有理由把 (A) 作为生产路线；若差异很大，再决定是否投入完整 \(\Omega_M\) 反馈。这样既尊重师兄对“先修正拉氏量”的优先级，也不把当前 partial feedback 的诊断价值丢掉。

### 6.3 推荐的最小可验证顺序

1. **耦合矩阵单元测试**：输入 \(\sigma_u=\sigma_d\) 时验证 \(K_{03}=K_{30}=K_{38}=K_{83}=0\)；输入小的 \(\sigma_u-\sigma_d\) 时验证线性响应和 \(K_{03}=K_{30}\)、\(K_{38}=K_{83}\) 的结构。
2. **中性谱回归**：完整 \(0/3/8\) 实现退回同位旋对称后，与现有 \(0/8\) 结果相同；检查极点、相移和介子质量的连续性。
3. **带电通道 A/B**：同一 \((T,\mu_B)\)、同一 x_min_cut/节点设置下，对比当前简化耦合与 full-\(\mathcal K_{ab}\) 的 \(\pi^\pm,K^\pm\) 密度及两种 ratio。
4. **守恒荷 A/B**：固定 \(T,\mu_B\) 后分别使用 quark-only \((\rho_Q/\rho_B=0.4,\rho_S=0)\) 和 partial conserved-charge correction，记录 \(\mu_Q,\mu_S\) 漂移与残差。
5. **只有在 1–4 通过后**，再评估是否需要 \(K^0,\bar K^0,\pi^0,\eta,\eta'\) 和完整 \(\Omega_M\) 反馈。否则很难区分是拉氏量、谱函数、积分门禁还是反馈造成的。

### 6.4 实验比较应怎样使用

本项目当前只需要实验的趋势和数量级，不需要逐点拟合。文献方法显示，以下因素本身就足以改变曲线形状：

- \(T,\mu_B\) 路径是手征线、伪临界线还是冻结线；
- \(\mu_s,\mu_\pi\) 是固定比例还是沿路径由守恒律给出；
- 使用稳定/BW 还是包含散射连续谱的 BU；
- 是否包含介子守恒荷或热力学反馈。

因此最稳妥的第一阶段产出是：在同一套模型口径下给出 \(K^+/\pi^+\)、\(K^-/\pi^-\) 的趋势、峰/平台/单调性和 \(O(10^{-1})\) 或 \(O(1)\) 数量级，并把实验的典型 \(O(10^{-1})\) 仅作为外部量级参照，不把它直接当成当前 x_min_cut 诊断结果的回归基线。

## 7. 文献到实现的缺口与后续任务边界

| 缺口 | 文献证据 | 最小实现 | 暂不应宣称 |
|---|---|---|---|
| \(u/d\) 不对称下的 \(K_{03},K_{38}\) | Rehberg 1996；Mei 2023 | 从独立 \(G^u,G^d,G^s\) 生成完整 \(K_{ab}\)，先接中性 \(3\times3\) RPA | “已经完成完整 PNJL 拉氏量” |
| 显式 \(\mu_I\) 上游耦合 | Xia 2013；Lopes 2021；Wu 2017 | 明确 \(\mu_I\) 归一化，决定是否允许 pion condensate，并扩展 MFA 变量 | “给带电极化函数传入 \(\mu_u-\mu_d\) 就等于完整 \(\mu_I\) 模型” |
| \(K^0,\bar K^0\) 等中性/非带电密度 | Friesen/Blaschke 的带电比值工作主要只需 \(K^\pm,\pi^\pm\) | 在确认中性通道的谱函数和守恒荷用途后再加入 | “当前四通道已代表所有介子荷” |
| 介子进入守恒荷 | 热统计/BQS 物理需要，当前分支仅 partial | 先保留 A/B 对照；若差异大，再设计稳定 \(\Omega_M\) 方案 | “partial feedback 就是完整反馈” |
| BU 支持、节点收敛和凝聚 | 当前项目已有 x_min_cut 诊断口径；公开 K/π 工作也常用有限窗口/相移近似 | 先把 full kernel 和 A/B 差异隔离，再单独做 strict-support gate | “本轮 ratio 已具备 production 物理完备性” |

## 8. 证据强度、限制和待补调研

### 8.1 证据强度

- **强**：Rehberg 全文直接给出 \(K_{03}=K_{30}\)、\(K_{38}=K_{83}\) 与 \(S^u-S^d\) 的关系，以及同位旋对称时的解耦条件；Mei 全文/公式片段给出同样的 \(\sigma_u-\sigma_d\) 依赖和 \(3\times3\) 中性混合。
- **中强**：Friesen、Blaschke 全文明确给出 \(\mu_s=0.55\mu_u\)、带电介子化学势差、极点/BW/BU 和路径比较，但这些论文不是完整 BQS+介子反馈计算。
- **中**：Xia、Lopes、Wu 证明显式 \(\mu_I\) 的 MFA/RPA 路线及可能的凝聚/混合，但目标分别是有限同位旋相图、EOS 或两味 pion 性质，不是本项目的 K/π 冻结线。
- **项目推论**：当前 EffectiveCouplings.jl、MesonPropagator.jl 和 MesonDensity.jl 的差距，是根据本仓库代码直接检查得到的实现事实；“完整交叉项对带电 ratio 的具体数值影响”仍需新一轮同点数值消融，不能从文献直接移植。

### 8.2 主要限制

1. 本次是问题导向的定向综述，不是覆盖所有年份和所有数据库的系统综述。
2. Rehberg 的一般 \(\mu_f\) 公式与其 \(\mu_f=0\) 数值结果必须区分；不能把其零密度谱当成冻结线基准。
3. Mei 的同位旋破缺来源是外磁场；它支撑的是矩阵构造，不支撑本项目具体的 \(\rho_Q/\rho_B=0.4\) 数值。
4. 不同作者对 \(\mu_I\) 有因子 2 的约定差异；任何代码接口都应显式命名 convention。
5. BU 的 \(x_{\min}\)、能量窗口、相移分支和凝聚处理会改变密度；因此文献间的 ratio 不能只按曲线高低比较。

### 8.3 待补的文献问题

- 需要进一步寻找在**三味 PNJL + BQS 守恒荷 + 显式带电/中性介子 BU + 冻结线**四者同时出现的论文；当前核心记录中没有一篇完全重合。
- 已补充两篇带介子热力学 \(\Omega_M\) 的 supporting 文献，但仍需要进一步核对不同广义 BU 工作采用的相移积分、Beth–Uhlenbeck virial 近似或 HRG 匹配，以及是否真的把 \(\Omega_M\) 对每个凝聚求导。
- 需要在本项目内做 full-\(\mathcal K_{ab}\) 与当前简化耦合的同点消融，再判断 ratio 接近 1 的主导来源是上游化学势、质量/谱权重，还是 partial feedback 的外层约束。

## 9. 核心与 supporting 文献清单

1. **Rehberg, Klevansky & Hüfner**, “Hadronization in the SU(3) NJL Model,” *Phys. Rev. C* **53**, 410 (1996). [DOI](https://doi.org/10.1103/PhysRevC.53.410) · [arXiv:hep-ph/9506436](https://arxiv.org/html/hep-ph/9506436).
2. **Xia, He & Zhuang**, “Three-flavor NJL model at finite isospin chemical potential,” *Phys. Rev. D* **88**, 056013 (2013). [DOI](https://doi.org/10.1103/PhysRevD.88.056013) · [arXiv:1307.4622](https://arxiv.org/html/1307.4622).
3. **Friesen, Kalinovsky & Toneev**, “Strange matter and kaon to pion ratio in SU(3) PNJL model,” *Phys. Rev. C* **99**, 045201 (2019). [DOI](https://doi.org/10.1103/PhysRevC.99.045201) · [arXiv:1808.04179](https://arxiv.org/html/1808.04179).
4. **Blaschke, Friesen, Kalinovsky & Radzhabov**, “Chiral phase transition and kaon-to-pion ratios in the entanglement SU(3) PNJL model,” *Eur. Phys. J. Spec. Top.* **229**, 3517–3536 (2020). [DOI](https://doi.org/10.1140/epjst/e2020-000218-1) · [arXiv:2009.03284](https://arxiv.org/html/2009.03284).
5. **Blaschke, Friesen, Kalinovsky & Radzhabov**, “Using the Beth-Uhlenbeck Approach to Describe the Kaon to Pion Ratio in a 2 + 1 Flavor PNJL Model,” *Particles* **3**, 169–177 (2020). [DOI](https://doi.org/10.3390/particles3010014) · [arXiv:1912.13162](https://arxiv.org/abs/1912.13162).
6. **Lopes et al.**, “Hot QCD at finite isospin density: confronting SU(3) NJL model with recent lattice data,” *Phys. Rev. D* **103**, 076023 (2021). [DOI](https://doi.org/10.1103/PhysRevD.103.076023) · [arXiv:2102.02844](https://arxiv.org/html/2102.02844).
7. **Wu et al.**, “Pion properties at finite isospin chemical potential with isospin symmetry breaking,” *Chinese Phys. C* **41**, 124106 (2017). [DOI](https://doi.org/10.1088/1674-1137/41/12/124106).
8. **Mei et al.**, “Mass spectra of neutral mesons \(K_0,\pi_0,\eta,\eta'\) at finite magnetic field, temperature and baryon chemical potential,” *Phys. Rev. D* **107**, 074018 (2023). [DOI](https://doi.org/10.1103/PhysRevD.107.074018) · [arXiv:2212.04778](https://arxiv.org/abs/2212.04778).
9. **Tian et al.**, “Mass spectra of charged kaons and pions at finite magnetic field, temperature, and density,” *Phys. Rev. D* **114**, 034012 (2026). [DOI](https://doi.org/10.1103/d7nm-y2vp).
10. **Blaschke et al.**, “Generalized Beth–Uhlenbeck approach to mesons and diquarks in hot, dense quark matter,” *Annals of Physics* **348**, 228–255 (2014). [DOI](https://doi.org/10.1016/j.aop.2014.06.002) · [arXiv:1305.3907](https://arxiv.org/abs/1305.3907).
11. **Maslov & Blaschke**, “Effect of mesonic off-shell correlations in the PNJL equation of state,” *Phys. Rev. D* **107**, 094010 (2023). [DOI](https://doi.org/10.1103/PhysRevD.107.094010) · [arXiv:2301.09882](https://arxiv.org/abs/2301.09882).

## 10. 本项目内的交叉引用

- [从拉氏量到巨热力学势：平均场推导统一流程](../../reference/formula/models/shared/Lagrangian_to_GrandPotential_MFA_Workflow.md)：本仓库的 MFA/巨势统一语义。
- [介子密度稳定粒子与 K/π 比值公式](../../reference/formula/relaxtime/meson_density/MesonDensity_稳定粒子与KPi比值.md)：稳定粒子密度、化学势差和比值定义。
- [介子密度 BU 相移公式](../../reference/formula/relaxtime/meson_density/MesonDensity_BU相移公式.md)：相移导数、窗口和积分约定。
- [部分反馈与冻结线粗扫任务记录](../../dev/active/2026-07-27_partial-feedback性能与冻结线粗扫.md)：当前 x_min_cut partial-feedback 诊断的成本、比值和生产边界。
- [带电 K/π BU kernel 门禁与同点消融](../../dev/active/2026-07-26_带电KPi_BU_kernel门禁与同点消融.md)：当前带电通道 kernel 门禁和 \(\mu_s\) 消融记录。
