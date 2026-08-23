## 【提取报告：外磁场下的PNJL模型相关公式】

### 📍 源信息
- **文献**：高雪艳博士论文《强相互作用物质相变与重子数涨落的研究》
- **位置**：第五章“磁场对QCD一阶相变的影响”，第5.1节“外磁场下的PNJL模型”
- **提取时间**：2026-2-20

---

## 🔢 核心公式组

### **公式1：外磁场下的PNJL模型拉格朗日密度**
**【标签】** 式(5-1)
**【类型】** 拉格朗日密度
**【内容】**
$$
\begin{aligned}
\mathcal{L}=&\bar{q}\left(i\gamma^{\mu}D_{\mu}+\gamma_{0}\hat{\mu}-\hat{m}_{0}\right)q \\
&+G\sum_{k=0}^{8}\left[\left(\bar{q}\lambda_{k}q\right)^{2}+\left(\bar{q}i\gamma_{5}\lambda_{k}q\right)^{2}\right] \\
&-K\left[\operatorname{det}_{f}\left(\bar{q}\left(1+\gamma_{5}\right)q\right)+\operatorname{det}_{f}\left(\bar{q}\left(1-\gamma_{5}\right)q\right)\right] \\
&-U(\Phi[A],\bar{\Phi}[A],T)-\frac{1}{4}F^{\mu\nu}F_{\mu\nu}
\end{aligned}
$$

**【上下文关联】**
- **引入文本**：“在外部磁场存在的情况下，夸克物质的拉格朗日密度具有如下形式”
- **关键修改**：协变导数增加电矢势项 $D^{\mu} = \partial^{\mu} - iA^{\mu} - iq_{i}A_{EM}^{\mu}$
- **磁场设定**：沿z方向静态恒定磁场 $A_{EM}^{\mu} = \delta^{\mu 2}x_{1}B$

---

### **公式2：磁场下PNJL模型的巨热力学势**
**【标签】** 式(5-2)
**【类型】** 巨热力学势表达式
**【内容】**
$$
\Omega = \sum_{f = u,d,s}\left(\Omega_{f}^{0} + \Omega_{f}^{T}\right) + 2G\left(\phi_{u}^{2} + \phi_{d}^{2} + \phi_{s}^{2}\right) -4K\phi_{u}\phi_{d}\phi_{s} + U(\Phi ,\bar{\Phi},T)
$$

---

### **公式3：真空项 $\Omega_{f}^{0}$**
**【标签】** 式(5-3)、(5-7)
**【类型】** 真空贡献的积分表达式
**【内容】**
$$
\Omega_{f}^{0} = -N_{c}\frac{|q_{f}|eB}{2\pi}\sum_{n = 0}^{\infty}\alpha_{n}\int_{-\infty}^{\infty}\frac{dp_{z}}{2\pi} E_{f,n}
$$

**平滑截断正则化形式**（式5-7）：
$$
\Omega_{f}^{0} = -N_{c}\frac{|q_{f}|eB}{2\pi}\sum_{n = 0}^{\infty}\alpha_{n}\int_{-\infty}^{\infty}\frac{dp_{z}}{2\pi} f_{\Lambda}^{2}(p_{f})E_{f,n}
$$
$$
f_{\Lambda}(p) = \sqrt{\frac{\Lambda^{2N}}{\Lambda^{2N} + p^{2N}}}
$$
其中 $N=10$ 在数值计算中采用。

实现中将截断函数的横向动量写为 $p_{f,n}=\sqrt{p_z^2+2n|q_f|eB}$；论文式(5-7)只记作 $f_\Lambda(p_f)$，这里的展开是代码与 Landau 能谱的变量映射，不是额外的物理假设。

---

### **公式4：温度相关项 $\Omega_{f}^{T}$**
**【标签】** 式(5-4)
**【类型】** 温度贡献的积分表达式
**【内容】**
$$
\Omega_{f}^{T} = -T\frac{|q_{f}|eB}{2\pi}\sum_{n = 0}^{\infty}\alpha_{n}\int_{-\infty}^{+\infty}\frac{dp_{z}}{2\pi}\left(Z_{f}^{+} + Z_{f}^{-}\right)
$$

**其中：**
$$
\begin{aligned}
Z_{f}^{+} &= \ln \left(1 + 3\Phi e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 3\bar{\Phi} e^{-2\frac{E_{f,n} - \mu_{f}}{T}} + e^{-3\frac{E_{f,n} - \mu_{f}}{T}}\right) \\
Z_{f}^{-} &= \ln \left(1 + 3\bar{\Phi} e^{-\frac{E_{f,n} + \mu_{f}}{T}} + 3\Phi e^{-2\frac{E_{f,n} + \mu_{f}}{T}} + e^{-3\frac{E_{f,n} + \mu_{f}}{T}}\right)
\end{aligned}
$$

---

### **公式5：准粒子能量（朗道能级）**
**【标签】** 式(5-6)
**【类型】** 能量本征值
**【内容】**
$$
E_{f,n} = \left(2n|q_{f}|B + p_{z}^{2} + M_{f}^{2}\right)^{1 / 2}
$$
其中 $n = 0,1,2,\ldots$ 表示朗道能级指标。

---

### **公式6：数密度公式**
**【标签】** 式(5-10)、(5-11)
**【类型】** 夸克数密度表达式
**【内容】**
$$
\rho_{f} = \frac{\partial\Omega}{\partial\mu_{f}} = \sum_{n = 0}^{\infty}\rho_{f,n}
$$
$$
\rho_{f,n} = 3\frac{q_{f}B}{2\pi}\alpha_{n}\int_{-\infty}^{\infty}\frac{\mathrm{d}p_{z}}{2\pi}\left[f_{f,n}^{+} - f_{f,n}^{-}\right]
$$

**其中分布函数：**
$$
\begin{aligned}
f_{f,n}^{+} &= \frac{\Phi e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 2\bar{\Phi} e^{-2\frac{E_{f,n} - \mu_{f}}{T}} + e^{-3\frac{E_{f,n} - \mu_{f}}{T}}}{\left(1 + 3\Phi e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 3\bar{\Phi} e^{-2\frac{E_{f,n} - \mu_{f}}{T}} + e^{-3\frac{E_{f,n} - \mu_{f}}{T}}\right)} \\
f_{f,n}^{-} &= \frac{\bar{\Phi}e^{-\frac{E_{f,n} + \mu_{f}}{T}} + 2\Phi e^{-2\frac{E_{f,n} + \mu_{f}}{T}} + e^{-3\frac{E_{f,n} + \mu_{f}}{T}}}{\left(1 + 3\bar{\Phi}e^{-\frac{E_{f,n} + \mu_{f}}{T}} + 3\Phi e^{-2\frac{E_{f,n} + \mu_{f}}{T}} + e^{-3\frac{E_{f,n} + \mu_{f}}{T}}\right)}
\end{aligned}
$$

---

### **公式7：低温近似下的数密度**
**【标签】** 式(5-13)
**【类型】** 近似表达式
**【内容】**
$$
\rho_{f,n}\approx 3\frac{q_{f}B}{2\pi^{2}}\alpha_{n}\int_{0}^{\infty}\frac{\mathrm{d}p_{z}}{(2\pi)}\left[\frac{1}{1 + e^{3\left(\frac{E_{f,n} - \mu_{f}}{T}\right)}}\right]
$$

---

### **公式8：磁场依赖的耦合常数**
**【标签】** 式(5-9)
**【类型】** 经验公式（用于描述IMC效应）
**【内容】**
$$
G(B) = G_{0}\frac{1 + a\zeta^{2} + b\zeta^{3}}{1 + c\zeta^{2} + d\zeta^{4}}
$$
其中 $\zeta = eB / \Lambda_{QCD}^{2}$，$\Lambda_{QCD} = 300\mathrm{MeV}$

---

## 开发者审核表：公式到实现的对应关系

本表是实现审计表，不是外部物理验证表。状态分为：

- **已接入（内核）**：公式已经进入候选分支的实际计算函数，并有 focused unit/regression 证据。
- **已接入（诊断）**：代码路径存在且可短样本运行，但尚未完成默认精度、全域收敛或分支覆盖验收。
- **仅源公式/未实现**：论文中有公式，但当前实现没有把它作为可用计算路径。
- **明确不包含**：当前合同主动排除，不能从现有输出推断出该能力。

审核对象分为两层：

1. **主线现状**：`main`/当前主工作树仍只有固定 `x_state` 的磁场内核与专题脚本；它没有
   `MagneticGapSolver.jl`，专题脚本也没有调用完整五维 equilibrium solver。
2. **隔离候选**：以下候选状态对应 `codex/core-algorithm-perf-ab@6697feddd1c6382955727b93881b68569c165779`
   及其未提交工作树改动；它们尚未进入 `main`，不构成当前生产能力。

| 物理/数值对象 | 论文或公式来源 | 候选实现入口 | 主线现状 | 隔离候选状态 | 开发者审核要点 |
|---|---|---|---|---|---|
| 外场协变导数与 Maxwell 项 | 式(5-1) | `PNJLMagneticModel.magnetic.eB_fm2` | 仅以标量 `eB` 传入背景场 | **明确不包含** Maxwell 自能 | 没有 gauge-field 求解，也没有 `B^2/2` 输出；不能解释为总 EOS |
| 巨热力学势分解 | 式(5-2) | `calculate_magnetic_omega_components` | 固定态内核已组装 | **已接入（诊断）** | 组装 `chi + poly + vac + therm`；尚无高节点 equilibrium 生产证据 |
| Landau 能谱 | 式(5-6) | `MagneticIntegrals.energy_landau` | **已接入（内核）** | **已接入（内核）** | `sqrt(pz^2 + M^2 + 2n|q|eB)`；有能谱单测 |
| Landau 简并度 | `alpha_n=2-delta_{n0}` | `MagneticIntegrals.alpha_n` | **已接入（内核）** | **已接入（内核）** | `n=0` 为 1，其余为 2；有单测 |
| 真空 Landau 项 | 式(5-3)、(5-7) | `omega0_flavor_landau` | **已接入（内核）** | **已接入（内核）** | 使用有限 `n_max`、pz Gauss 节点和平滑截断；数值截断需独立收敛 |
| 平滑截断 | 式(5-7) | `MagneticIntegrals.smooth_cutoff` | **已接入（内核）** | **已接入（内核）** | 默认 `cutoff_N=10`；论文公式不决定 pz 节点和截断上限 |
| Polyakov 热项 | 式(5-4) 及 `Z_f^+ + Z_f^-` | `omegat_flavor_landau`、`_log_polyakov_pair` | **已接入（内核）** | **已接入（内核）** | 使用一般温度双对数；低温指数稳定缩放已有短测试 |
| 动力学质量与手征项 | 式(2-55)、式(5-2)中的 `G(eB)` | `_calculate_mass_vec_with_GB`、`_chiral_with_GB` | 固定态内核已接入 | **已接入（诊断）** | 质量核支持泛型实数；完整 Omega 仍以 Float64 finite-difference solver 为主 |
| 对数 Polyakov 势 | 式(2-50)--(2-51) | `PNJLCore.polyakov_potential` | 固定态内核已接入 | **已接入（诊断）** | 由 PNJL core 复用；未单独证明磁场全域参数适用性 |
| 一般温度密度 | 式(5-10)--(5-12) | `density_flavor_landau`、`calculate_magnetic_rho` | **已接入（内核）** | **已接入（内核）** | 按用户确认采用 `q-qbar` 净密度和 `|q_f eB|` 相空间；独立差分交叉检查最大误差约 `6.1e-13` |
| 低温三夸克近似 | 式(5-13) | 无独立默认实现 | **仅源公式/未实现** | **仅源公式/未实现** | 不能把论文近似当作默认密度；磁场 equilibrium 明确要求 `T_fm>0` |
| 五维驻点条件 | 式(2-63) | `magnetic_gap_residual`、`solve_magnetic_gap` | **未接入主线** | **已接入（诊断）** | 候选 solver 用五维有限差分 residual；有低节点 stationarity probe，未完成生产精度验证 |
| 多分支候选 | 论文低温/高场多解讨论；五变量驻点条件 | `MagneticGapCandidate`、`MagneticGapResult` | **未接入主线** | **已接入（诊断）** | 多 seed、去重、Omega 排序可用；不能证明 seed 覆盖全部分支 |
| 局部稳定性分类 | 由 `Omega` Hessian 判定局部极小 | `_magnetic_hessian`、`classify_stability=true` | **未接入主线** | **已接入（诊断）** | 当前为有限差分 Hessian；低节点出现小 residual 但 `saddle_or_maximum`；该标签不作为 PNJL 默认生产过滤条件 |
| `G(eB)` IMC 参数化 | 式(5-9) | `coupling_GB` | **已接入（内核）** | **已接入（内核）** | 参数和单位已固定；不等于磁场 EOS 或磁化响应完整实现 |
| `n_max` 数值起点 | 论文低温占据估计 + 实现 cutoff 约定 | `resolve_nmax_from_cutoff` | **已接入（诊断）** | **已接入（诊断）** | 当前是 cutoff-based 求和起点，不是低温物理占据上限 |
| Landau 截断收敛 | 数值治理，不是论文单一公式 | `magnetic_nmax_convergence_report` | **已接入（诊断）** | **已接入（诊断）** | 只比较 `n_base` 与增量截断；不能替代全域收敛研究 |
| 标量物质压力 | `P=-Omega` | `calculate_magnetic_pressure` | **已接入（内核）** | **已接入（内核）** | 不含 Maxwell、磁化强度、纵向/横向压力张量 |
| 独立夸克/反夸克密度接口 | 通用 `Models.number_densities` 合同 | `PNJLMagneticModel.number_densities` | 旧主线宣称支持 | 非零 `eB` capability 已收窄为不支持 | 非零磁场只返回 `net`；不能直接接普通 transport workflow |
| RS 联合磁场路线 | 当前项目作者合同 | `xi` 参数校验 | 未实现 | **明确不包含** | 非零 `xi` 显式拒绝；磁场不复用 RS reduced adaptive 路线 |
| 磁场专题工作流 | 项目脚本合同 | `run_magnetic_point.jl`、`run_magnetic_eb_scan.jl`、`run_magnetic_stability_scan.jl` | 固定 `x_state` 诊断 | 候选 equilibrium API 另行存在 | 脚本输出不能解释为完整五维生产扫描 |

### 表格结论

从公式覆盖看，隔离候选已经覆盖“固定外部磁场背景下的物质巨势、一般温度净密度和五维驻点接口”；从生产资格看，仍不能称为完整磁场计算链路，原因是：

1. 主线尚未接入五维 solver；
2. 候选 solver 的多 seed 不是全分支证明；
3. Hessian 分类只是可选诊断，尚无高节点全域稳定性证据；
4. 专题脚本仍使用固定 `x_state`；
5. Maxwell、磁化响应、方向性压力、`T=0` 近似和 RS 联合路线均不在当前合同内。

因此“公式已实现”只能用于描述内核覆盖，不能替代“已完成生产验证”。Hessian 标签应作为分支审计和显式研究信息保留，不应在没有额外物理约定时升级为默认生产资格门槛。

---

## 🔧 变量说明表

| 符号 | 物理意义 | 单位/量纲 | 备注 |
|------|----------|-----------|------|
| $\Omega$ | 巨热力学势 | MeV⁴ |  |
| $\Omega_{f}^{0}$ | 真空贡献项 | MeV⁴ | 每味夸克 |
| $\Omega_{f}^{T}$ | 温度贡献项 | MeV⁴ | 每味夸克 |
| $q_{f}$ | 夸克电荷 | e | u:+2/3, d:-1/3, s:-1/3 |
| $eB$ | 磁场强度 | MeV² |  |
| $B$ | 磁感应强度 | G | 1 G ≈ 1.95×10⁻²⁰ GeV² |
| $n$ | 朗道能级指标 | 无量纲 | n=0,1,2,... |
| $\alpha_{n}$ | 简并度因子 | 无量纲 | n=0时为1，n≠0时为2 |
| $E_{f,n}$ | 准粒子能量 | MeV | 依赖于朗道能级 |
| $M_{f}$ | 组分夸克质量 | MeV | 由于征对称性破缺产生 |
| $\Phi,\bar{\Phi}$ | Polyakov圈及其共轭 | 无量纲 | 禁闭序参量 |
| $U(\Phi,\bar{\Phi},T)$ | Polyakov圈有效势 | MeV⁴ | 对数形式 |
| $G$ | 四费米子耦合常数 | MeV⁻² | 磁场依赖 |
| $K$ | 六费米子耦合常数 | MeV⁻⁵ | t'Hooft相互作用 |
| $\phi_{i}$ | 夸克凝聚 | MeV | i=u,d,s |
| $f_{\Lambda}(p)$ | 截断函数 | 无量纲 | 平滑截断用 |
| $\Lambda$ | 动量截断 | MeV | 602.3 MeV |
| $\Lambda_{QCD}$ | QCD能标 | MeV | 300 MeV |
| $Z_{f}^{\pm}$ | 统计分布的对数形式 | 无量纲 | 包含Polyakov圈耦合 |
| $f_{f,n}^{\pm}$ | 夸克/反夸克分布函数 | 无量纲 | 含Polyakov圈修正 |

---

## 📊 参数表格

**表5-1：磁场依赖耦合常数G的相关参数**
| 参数 | 数值 | 单位 | 说明 |
|------|------|------|------|
| a | 0.108805 | 无量纲 | 式(5-9)参数 |
| b | -1.0133 × 10⁻⁴ | 无量纲 | 式(5-9)参数 |
| c | 0.02228 | 无量纲 | 式(5-9)参数 |
| d | 1.84558 × 10⁻⁴ | 无量纲 | 式(5-9)参数 |
| $G_{0}$ | - | MeV⁻² | 零磁场时的耦合常数 |
| $\Lambda_{QCD}$ | 300 | MeV | QCD能标 |

**PNJL模型固定参数（第二章）**
| 参数 | 数值 | 单位 | 说明 |
|------|------|------|------|
| $\Lambda$ | 602.3 | MeV | 三动量截断 |
| $G\Lambda^{2}$ | 1.835 | 无量纲 | 四费米子耦合 |
| $K\Lambda^{5}$ | 12.56 | 无量纲 | 六费米子耦合 |
| $m_{u}$ | 5.5 | MeV | u夸克流质量 |
| $m_{d}$ | 5.5 | MeV | d夸克流质量 |
| $m_{s}$ | 140.7 | MeV | s夸克流质量 |

---

## 💡 数值实现建议

1. **朗道能级求和**：
   - 对于给定的化学势 $\mu_{f}$，填充的最高朗道能级为 $n_{\mathrm{max}} = \mathrm{Floor}\left(\frac{\mu_{f}^{2} - M_{f}^{2}}{2|q_{f}|B}\right)$
   - 求和需包含足够多的能级以确保收敛

2. **真空项处理**：
   - 使用平滑截断函数 $f_{\Lambda}(p)$ 避免截断不连续性
   - 建议采用N=10以获得良好的平滑效果

3. **低温数值处理**：
   - 低温时（$T<20$ MeV），$\Phi \approx \bar{\Phi}\approx 0$，可简化分布函数
   - 需注意数值积分在费米面附近的奇异性

4. **多重一阶相变**：
   - 低温高磁场下可能出现多个一阶相变结构
   - 求解能隙方程时需小心处理多解情况

5. **IMC效应实现**：
   - 需采用磁场依赖的耦合常数 $G(B)$ 式(5-9)
   - 参数a,b,c,d已调至重现格点QCD的IMC效应

---

## 🔗 相关公式链

1. **能隙方程**：式(2-55) - 夸克动力学质量 $M_i = m_i - 4G\phi_i + 2K\phi_j\phi_k$
2. **稳定平衡条件**：式(2-63) $\frac{\partial\Omega_{PNJL}}{\partial\phi_u} = \frac{\partial\Omega_{PNJL}}{\partial\phi_d} = \frac{\partial\Omega_{PNJL}}{\partial\phi_s} = \frac{\partial\Omega_{PNJL}}{\partial\Phi} = \frac{\partial\Omega_{PNJL}}{\partial\bar{\Phi}} = 0$
3. **Polyakov圈有效势**：式(2-50)、(2-51) - 对数形式有效势
4. **重子数密度**：$\rho_B = (\rho_u + \rho_d + \rho_s)/3$

## 📚 原文核验与符号边界

当前实现沿用外部固定磁场的物质巨势口径：式(5-1)中的外部 Maxwell 项用于定义背景场，但没有被加到 `omega/pressure` 的返回值中。因而当前压力不是包含 (B^2/2) 电磁场自能的总压力；如需总 EOS，必须单独定义 Maxwell 项及其压力各向异性约定。

此外，当前 `calculate_magnetic_pressure=-Omega` 只返回一个固定外部磁场背景下的标量物质压力；实现没有计算磁化强度，也没有给出纵向/横向压力或完整压力张量。后续若需要磁化介质 EOS，必须新增明确的方向性压力合同并单独验证。

本页所引原文为高雪艳博士论文《强相互作用物质相变与重子数涨落的研究》：第 2 章第 2.2 节（印刷页 21--24，式(2-50)--(2-70)）和第 5 章第 5.1 节（印刷页 65--68，式(5-1)--(5-13)）。第 5.1 节明确给出磁场巨势、Landau 能谱、平滑截断、Polyakov 热项、$G(eB)$ 参数化和 $\mu_u=\mu_d=\mu_s$ 的计算约定；第 2.2 节给出对数 Polyakov 势、参数表、动力学质量和五变量驻点条件。

论文版面本身存在需要保留的审计注记：式(2-63) 的等式链把 $\partial\Omega/\partial\phi_s$ 排成了重复的 $\partial\Omega/\partial\phi_u$；式(5-10) 的密度导数符号与第 2 章式(2-65)的热力学约定不一致；式(5-11) 使用了带符号的 $q_f B$，而真空/热项使用 $|q_f|eB$。因此这些字符不能在没有约定确认时直接作为实现合同；本页对反夸克分布函数按论文原页版式更正了明显的抽取错位，但没有替作者决定密度符号或电荷绝对值的最终约定。

数值接口的 `p_num`、`pz_max`、Landau 求和起点、`delta_n`/`rtol` 收敛判据以及 quadrature 实现均不是论文式(5-1)--(5-13)给出的物理公式，必须以代码/API 合同单独记录。论文给出的低温 `n_max` 是占据能级估计；它不自动等同于平滑截断真空项所需的数值求和上限。
