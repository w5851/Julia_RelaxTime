# PNJL/Mott 相变：因变量随 xi 的变化分析

## 1. 分析范围与数据来源

- 主数据（数值）：`data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv`
- 运行配置与扫描参数：
  - `data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/effective_config.json`
  - `config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml`
- 已有图像：
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a/mott_mode_a__xi-0p3.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a/mott_mode_a__xi-0p15.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a/mott_mode_a__xi0.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a/mott_mode_a__xi0p15.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a/mott_mode_a__xi0p3.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_b/mott_mode_b__M_pi.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_b/mott_mode_b__M_K.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_b/mott_mode_b__Gamma_pi.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_b/mott_mode_b__Gamma_K.png`
- 文献对比目标：`tests/validation/data/targets/relaxtime/literature/mott/relaxtime_mott_temperature_literature_targets_v1.csv`

说明：本批数据为 `muB=0`、`xi in {-0.3,-0.15,0,0.15,0.3}` 的温度扫描，实际输出覆盖 `T=100~300 MeV`（步长 1 MeV）。

## 2. 因变量定义

本报告将 **Mott 临界温度** 作为主因变量：

- `T_mott^pi`: 满足 `M_pi(T,xi)=M_u(T,xi)+M_d(T,xi)` 的温度
- `T_mott^K`: 满足 `M_K(T,xi)=M_u(T,xi)+M_s(T,xi)` 的温度

并以介子宽度与有效质量作为辅助观测量：

- `Gamma_pi(T,xi)`, `Gamma_K(T,xi)`（束缚态进入连续谱后出现非零宽度）
- `m_u(T,xi)`, `m_s(T,xi)`（反映手征恢复进程和阈值位置变化）

### 2.1 xi 的定义（来自公式文档）与实现口径

在公式文档 `docs/reference/formula/models/pnjl/distribution/PNJL_夸克有效分布函数_动量各向异性.md` 中，采用 Romatschke-Strickland 形式：

```math
f^{\text{aniso}}(\mathbf{p}) = f^{\text{iso}}\left(\sqrt{\mathbf{p}^2 + \xi(\mathbf{p}\cdot\mathbf{n})^2}\right)
```

对应到代码实现（`src/QuarkDistribution_Aniso.jl`）是：

```math
E_{\text{aniso}}=\sqrt{p^2+m^2+\xi(p\cos\theta)^2}
```

即 xi 通过修正能量壳进入夸克/反夸克分布函数。弱各向异性一阶项在代码中也显式实现为
`delta f ~ (xi/2)*(p cos(theta))^2/E * (df/dE)`。

注：仓库不同文档对“xi 正负对应的几何描述”存在表述差异（例如 `docs/api/data_contracts.md` 与 formula 文档）。

从常用定义
`xi = <p_perp^2> / (2<p_parallel^2>) - 1`
与实现中的
`E_aniso = sqrt(p^2 + m^2 + xi*p_parallel^2)`
联合看：

- `xi > 0`：对平行方向动量的“代价”增加，分布在 `n` 方向相对收缩（oblate，横向相对拉伸）。
- `-1 < xi < 0`：平行方向相对拉伸（prolate）。

因此若问“`xi>0` 在 `n` 方向是拉伸还是收缩”，在本实现口径下应理解为**沿 `n` 方向收缩、横向拉伸**。

## 3. 关键数值结果

### 3.1 xi 与 Mott 温度（主结果）

用线性插值提取阈值（在相邻温度点上对 `M_mes - M_threshold` 过零点插值）后：

| xi | T_mott^pi (MeV) | T_mott^K (MeV) |
|---:|---:|---:|
| -0.30 | 196.474 | 195.073 |
| -0.15 | 202.241 | 201.127 |
| 0.00 | 207.288 | 206.422 |
| 0.15 | 211.960 | 211.131 |
| 0.30 | 216.008 | 215.367 |

趋势：

- 两个通道都随 xi 增大而 **单调上升**。
- 线性拟合斜率：
  - `dT_mott^pi/dxi ≈ +32.52 MeV`
  - `dT_mott^K/dxi ≈ +33.73 MeV`
- 存在轻微“饱和”倾向（非线性）：
  - 二次项系数均为负（`a_pi≈-12.10`, `a_K≈-13.40`），
  - 相邻 xi 台阶增量逐步减小（例如 `T_mott^pi` 的增量从 `+5.77` 降至 `+4.05 MeV`）。

### 3.2 宽度与有效质量的辅助证据

在固定温度 `T=220 MeV`：

| xi | Gamma_pi (fm^-1) | Gamma_K (fm^-1) | m_u(T=200) (fm^-1) | m_s(T=200) (fm^-1) |
|---:|---:|---:|---:|---:|
| -0.30 | 0.7393 | 0.8821 | 0.2796 | 2.0301 |
| -0.15 | 0.5841 | 0.7611 | 0.4830 | 2.0977 |
| 0.00 | 0.4570 | 0.6513 | 0.7272 | 2.1775 |
| 0.15 | 0.3440 | 0.5357 | 0.9089 | 2.2479 |
| 0.30 | 0.2288 | 0.3899 | 1.0342 | 2.3028 |

趋势：

- `Gamma_pi`, `Gamma_K` 随 xi 增大明显下降（在同一 T 上离连续谱“更远”）。
- `m_u`, `m_s` 随 xi 增大而上升，抬升了 `q\bar q` 连续阈值（`M_u+M_d`, `M_u+M_s`）。

### 3.3 介子质量与衰变宽度随 xi 的关系（分温区）

为避免“只看单一温度点”的偏差，这里给出 `T=180/200/220 MeV` 的横向对比。

#### (a) 低于阈值区（T=180 MeV）

- `M_pi` 随 xi 升高：`0.588 -> 0.716 fm^-1`；`M_K` 也小幅升高：`2.481 -> 2.512 fm^-1`。
- 但阈值升得更快：`M_u+M_d: 2.494 -> 3.117 fm^-1`，`M_u+M_s: 3.672 -> 4.152 fm^-1`。
- 两个通道都远低于阈值，因此 `Gamma_pi = Gamma_K = 0`。

#### (b) 过渡区（T=200 MeV）

- 负 xi 一侧（`xi=-0.3`）已解束缚：
  - `M_pi-(M_u+M_d)=+0.310`，`Gamma_pi=0.195`
  - `M_K-(M_u+M_s)=+0.300`，`Gamma_K=0.408`
- 非负 xi 一侧（`xi>=-0.15`）尚未过阈：`Gamma_pi=Gamma_K=0`。

这说明 xi 改变的不只是“宽度大小”，还改变“是否进入连续谱”的相态判据。

#### (c) 高于阈值区（T=220 MeV）

- 所有 xi 都已进入连续谱，但随 xi 增大：
  - 阈值距离缩小：
    - `Delta_pi = M_pi-(M_u+M_d): 1.302 -> 0.303`
    - `Delta_K  = M_K-(M_u+M_s): 0.856 -> 0.238`
  - 宽度同步下降：
    - `Gamma_pi: 0.739 -> 0.229`
    - `Gamma_K : 0.882 -> 0.390`

结论：在同一高温点，xi 增大使系统“离阈值更近但仍在阈值上方”，相空间压缩导致衰变宽度减小。

### 3.4 与文献目标点对比（muB=0, xi=0）

- 文献目标：`T_mott^pi = 205.6 MeV`（`rtol=2%`）
- 本数据提取：`T_mott^pi ≈ 207.3 MeV`（插值）
- 相对差异：约 `+1.17%`，落在目标容差内。

### 3.5 序参量（手征 + Polyakov）数值链条

为补足“xi -> 序参量 -> Mott”的中间环节，使用
`data/outputs/results/relaxtime/scan/gap_transport_scan_step5_muB0_xi-0p3to0p3.csv`
中的平衡解（`Phi, m_u, m_s`）并在 `xi=-0.15,0.15` 处做 xi 线性插值。

#### (a) 固定温度切片（直接看 xi 影响）

`T=200 MeV` 时：

- `Phi: 0.632 (xi=-0.3) -> 0.583 (xi=+0.3)`（下降）
- `m_u: 0.2796 -> 1.0342 fm^-1`（显著上升）
- `m_s: 2.0301 -> 2.3028 fm^-1`（上升）

`T=220 MeV` 时：

- `Phi: 0.703 -> 0.675`（继续下降）
- `m_u: 0.1127 -> 0.3461 fm^-1`

这说明在相同 T 下，xi 增大对应更小的 Polyakov 环与更大的夸克有效质量，数值上与“手征软化被推迟”一致。

#### (b) 沿各自 Mott 点读取序参量（看“触发条件”而非固定 T）

在 `T_mott^pi(xi)` 上（由 `mott_phase_derived.csv` 提取）：

| xi | T_mott^pi (MeV) | Phi(T_mott^pi) | m_u(T_mott^pi) | m_s(T_mott^pi) |
|---:|---:|---:|---:|---:|
| -0.30 | 196.47 | 0.6153 | 0.3830 | 2.0746 |
| -0.15 | 202.24 | 0.6290 | 0.4139 | 2.0689 |
| 0.00 | 207.29 | 0.6408 | 0.4293 | 2.0621 |
| 0.15 | 211.96 | 0.6518 | 0.4396 | 2.0546 |
| 0.30 | 216.01 | 0.6607 | 0.4502 | 2.0495 |

观察：尽管在固定 T 下 `Phi` 随 xi 下降，但由于 `T_mott` 本身上移，沿“各自阈值轨迹”读取时 `Phi(T_mott)` 反而上升。这反映了两个效应叠加：

1. xi 的直接效应（同温抑制 `Phi`）；
2. `T_mott` 上移带来的间接效应（高温使 `Phi` 变大）。

### 3.6 介子质量 `M_pi` / `M_K` 的数值链条

#### (a) 低温束缚区（T=180 MeV）

- `M_pi` 随 xi 增加（斜率 `dM_pi/dxi ≈ +0.207 fm^-1`）
- `M_K` 也小幅上升（`dM_K/dxi ≈ +0.051 fm^-1`）
- 但阈值上升更快，因此仍处稳定束缚态，`Gamma≈0`。

#### (b) 过渡温区（T=190~210 MeV）

- `M_pi(xi)` 由正斜率转为负斜率：
  - `T=190`：`dM_pi/dxi ≈ +0.225`
  - `T=200`：`dM_pi/dxi ≈ -0.116`
  - `T=210`：`dM_pi/dxi ≈ -0.605`
- `M_K(xi)`在 `T=190` 近零斜率，`T>=200` 后转负。

这对应“阈值重排与极点重定位”竞争：在接近/跨越 Mott 区间时，介子极点（实部）对 xi 的响应方向发生切换。

#### (c) 连续谱区（T>=220 MeV）

- `M_pi`, `M_K` 均随 xi 增大而下降（例如 `T=220`：
  `dM_pi/dxi≈-0.885`, `dM_K/dxi≈-0.305`）。
- 同时 `Gamma_pi`, `Gamma_K` 也下降（`T=220`：约 `-0.84`, `-0.81`）。

结论：高温连续谱区中，xi 增大使极点实部与虚部均向“更窄、相对更低质量”方向移动，但阈值上移仍主导 Mott 点整体右移。

## 4. 物理解释（基于 PNJL 图像）

### 4.1 PNJL 中 Mott 条件的核心

Mott 转变发生在介子束缚态质量与二夸克连续阈值相交处：

- 低温：手征破缺强，夸克有效质量较大，介子通常为窄束缚态（宽度近零）。
- 升温：手征凝聚减弱、Polyakov 环上升（退禁闭增强），介质激发增多；当 `M_mes >= M_q + M_\bar q` 时，束缚态溶解进入连续谱，宽度开启。

### 4.2 数值影响链：从“改变 xi”到“改变 Mott 温度”

下面给出一条可追踪、每一步有数字支撑的链条（muB=0）。

1. **xi 进入分布函数色散项**  
   `E_aniso = sqrt(p^2 + m^2 + xi*(p cos(theta))^2)`（`src/QuarkDistribution_Aniso.jl`）。

2. **热占据与自洽 gap 解被重排**（在本扫描窗口内体现为同温质量更大）  
   例如 `T=200 MeV`：
   - `m_u: 0.2796 (xi=-0.3) -> 1.0342 (xi=+0.3)`
   - `m_s: 2.0301 (xi=-0.3) -> 2.3028 (xi=+0.3)`

3. **连续阈值抬升**（速度快于介子质量）  
   例如 `T=220 MeV`：
   - `M_u+M_d: 0.2253 -> 0.6922`
   - `M_u+M_s: 1.9137 -> 2.3486`

4. **介子-阈值间隙 Delta 变小，宽度被压缩**  
   在 `T=220 MeV`：
   - `Delta_pi: 1.302 -> 0.303`，`Gamma_pi: 0.739 -> 0.229`
   - `Delta_K : 0.856 -> 0.238`，`Gamma_K : 0.882 -> 0.390`

5. **跨阈值时刻后移，Mott 温度升高**  
   - `T_mott^pi: 196.47 -> 216.01 MeV`（`xi: -0.3 -> +0.3`）
   - `T_mott^K : 195.07 -> 215.37 MeV`（同上）
   - 斜率约 `+33 MeV / xi`。

因此，本数据支持的主结论是：**xi 增大 -> 夸克有效质量与阈值上移 -> 宽度减小/开启点后移 -> Mott 温度上升**。

### 4.3 与 Polyakov 环/退禁闭的关系

PNJL 中 Polyakov 环控制色背景对夸克热激发的抑制与释放；xi 不直接等同于 Polyakov 环，但通过改变分布函数与自洽解，间接改变了“退禁闭推进速度”与手征软化进程的相对时序。当前结果对应的是：

- 在给定 T 下，体系保持更“重”的准粒子质量谱，
- 因此需要更高温度才达到同等解束缚条件。

### 4.4 深入到极化函数层：xi 如何直接进入介子极点方程

从实现链路看：

1. `MesonMass.solve_meson_mass` 解的是复极点方程  
   `1 - 4K * Π(k0+i*Gamma/2) = 0`（`src/relaxtime/MesonMass.jl`）。

2. `Π` 由 `polarization_with_width` 给出（`src/relaxtime/PolarizationAniso.jl`）：
   - 结构为 `A1 + A2 + prefactor * B0`，其中
   - `prefactor = k^2 - lambda^2 + (m1 +/- m2)^2 + gamma^2/4`。

3. xi 进入 `B0_correction(...)`（同文件第 88-96 行），即
   `B0 -> B0 + delta_B0(xi)`；这与公式文档
   `docs/reference/formula/relaxtime/polarization/Polarization_极化函数byB0.md`
   一致。

4. 因此 xi 对 `M` 与 `Gamma` 的影响有两条主通道：
   - **间接通道**：xi 改变 `m_u,m_s,Phi`，进而改变 `A/B0/prefactor`。
   - **直接通道**：xi 直接进入 `delta_B0`，改变 `Re Π` 与 `Im Π`。

在本批数据中我们已清楚看到间接通道足够强（阈值显著上升）；而“直接通道 vs 间接通道”的定量拆分还缺少逐点 `Π` 输出。

### 4.5 基于公式的极化函数层理论分析（结合当前数据）

在不新增 `Pi/B0/A` 逐点输出的前提下，可基于现有公式与观测数据给出“方向性+约束性”分析：

1. **极点条件与结构**  
   `MesonMass.jl` 中赝标量通道满足  
   `1 - 4K * Π = 0`，且（`k_norm=0, mu1=mu2`）
   `Π ~ -(Nc/8pi^2) * [ A1+A2 + ((m1-m2)^2 - M^2 + Gamma^2/4)*B0 + width-mix ]`。

2. **A 项（单线积分）对 xi 的效应符号**  
   按 `OneLoopIntegral_A_Aniso_Collinear.md` 与 `OneLoopIntegralsAniso.jl`，
   `A = A_iso + A_corr(xi)`，其中一阶修正由 `df/dE` 驱动。结合当前数据中
   `xi` 增大导致同温 `m_u,m_s` 上升、`Phi` 下降，可推断 A 项在过渡区被重排，
   等效体现为“需要更高 T 才满足同一极点条件”。这与 `T_mott` 上移一致。

3. **B0 项（双线积分）与 xi 直接通道**  
   `B0 -> B0 + delta_B0(xi)`（`B0_correction`），公式见
   `OneLoopIntegral_B0_Aniso_Collinear.md` / `OneLoopIntegral_B0_Aniso_General.md`：
   修正项与 `xi * p^2/E * (df/dE)` 成正比，且含角积分核（对数项与主值项）。
   因此 xi 不仅通过 `m,Phi` 间接影响极点，还直接改变 `ReB0/ImB0`。

4. **为何会出现“质量与宽度同时下降、但 Mott 温度上升”**  
   数据在 `T>=220 MeV` 显示：`M_pi,M_K` 与 `Gamma_pi,Gamma_K` 随 xi 降低，
   但阈值 `M_u+M_d`, `M_u+M_s` 上升更快，导致 `Delta=M_mes-M_thr` 变小却仍>0，
   并把 `Delta=0` 的交点整体推向更高温。  
   这说明在当前参数区，**阈值通道（由手征质量控制）主导 Mott 位移**，
   而极点本身的位移（由 A/B0 共同决定）是次级修正。

5. **当前可验证与不可验证边界**  
   - 可验证：`xi -> (m_u,m_s,Phi) -> threshold -> Delta/Gamma/T_mott` 的数值链条已闭合。  
   - 尚不可完全定量：`delta_B0(xi)` 与 `A_corr(xi)` 在 `ReΠ/ImΠ` 中各占多少比例。
   这需要扫描脚本额外输出极化函数诊断列。

### 4.6 `xi` 对 `Π` 与 `K` 的分工影响（介子层）

介子极点由 `1 - 4K*Π = 0` 决定，因此可分为两条影响路径：

1. **对 `Π`（极化函数）的影响**
   - **直接项**：`xi` 显式进入 `B0_correction(…, xi)`，直接改 `ReB0/ImB0`，再进入 `ReΠ/ImΠ`。
   - **间接项**：`xi` 先改 `m_u,m_s,Phi,PhiBar`，再通过 `A` 与 `B0` 改 `Π`。
   - 当前数据支持“间接项在阈值层面更强”：`m_u,m_s` 随 xi 上升显著，阈值抬升主导 Mott 右移。

2. **对 `K`（有效耦合）的影响**
   - `K123_plus/K4567_plus/...` 并非常数，而由 `G_u,G_s` 决定；
   - `G_f = -(Nc/(4pi^2)) * m_f * A_f`（`EffectiveCouplings.jl`），故 xi 通过 `m_f` 与 `A_f` 重排影响 `G_f`，进而影响各通道有效耦合。
   - 对 π/K 通道，进入极点方程的是 `K123_plus`/`K4567_plus`，即 xi 会通过“耦合重整化”改变极点位置；但在本批数据中，
     从结果看阈值抬升效应仍是主导项，`K` 的变化更像二级修正。

简言之：`xi` 同时改 `Π` 和 `K`，但当前 `muB=0` 扫描下，决定 Mott 位移方向与量级的主导机制是“手征质量-阈值链条”，而不是单独的 `K` 漂移。

## 5. 计算流程（声明式流水线）

本批 Mott 数据可整理为如下声明式 pipeline：

1. **Profile 装载**  
   读取 `config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml`，声明扫描域：
   `xi_list, T_min/max/step, muB, p_num, t_num, max_iter`。

2. **运行时配置物化**  
   `scripts/relaxtime/run_mott_phase_scan.jl` 将最终配置写入
   `effective_config.json` 与 `run_manifest.json`（含 `config_hash/run_id`）。

3. **主循环（xi × T 网格）**  
   对每个 `(xi,T)` 调用 `solve_gap_and_meson_point`：
   - 子阶段 A：`EquilibriumFacade.solve_equilibrium_backend` 得到平衡解（`Phi, masses`）
   - 子阶段 B：`MesonMass.solve_meson_mass` 解复极点方程，返回 `M, Gamma, residual`
   - 子阶段 C：`MottTransition` 计算阈值与 gap（workflow 内部用于诊断）

4. **原始结果落盘**  
   输出 `mott_phase_scan.csv`（本批为 `reference_100_300_fine/mott_phase_scan.csv`）。

5. **派生量构造**  
   `scripts/relaxtime/run_mott_phase_derived_csv.jl` 生成
   `M_u_plus_M_d`, `M_u_plus_M_s`，得到 `mott_phase_derived.csv`。

6. **可视化阶段**  
   `run_mott_phase_plot_modes.jl` 生成 `mode_a/mode_b` 图，用于查看交点与宽度开启。

7. **分析阶段（本报告）**  
   以 `derived.csv + figures + literature target` 组合分析 xi 依赖与物理链条。

## 5. 图像读取要点

- `mode_a`（按 xi 分图）可直接观察同一 xi 下 `M_pi/M_K` 与阈值线的交点温度。
- `mode_b`（按量分图）可横向比较不同 xi：
  - `mott_mode_b__M_pi.png`, `mott_mode_b__M_K.png`：交点随 xi 右移（更高 T）。
  - `mott_mode_b__Gamma_pi.png`, `mott_mode_b__Gamma_K.png`：宽度开启位置随 xi 右移，且固定 T 宽度随 xi 降低。

## 6. 未解决问题与后续建议

1. **直接输出极化函数诊断量**：建议在 `scripts/relaxtime/run_mott_phase_scan.jl` 中追加列：
   `RePi_pi, ImPi_pi, RePi_K, ImPi_K, ReB0, ImB0, dB0_xi`，用于分离“xi 通过质量/Phi 的间接效应”与“xi 进入 B0_correction 的直接效应”。
2. **序参量-阈值-极点同点记录**：对每个 `(xi, T_mott)` 保存 `Phi, m_u, m_s, M_mes, Gamma_mes, Delta_mes`，减少跨文件插值误差。
3. **muB 依赖联动**：本报告固定 `muB=0`；建议扩展到多 `muB` 平面，检验 `dT_mott/dxi` 是否在临界区域附近增强或反号。
4. **文献覆盖度**：目前 literature target 仅含 `xi=0` 的对照点；若要做更强结论，建议补充含各向异性参数扫描的外部结果并统一单位与误差模型。

## 7. paper 侧进展回写（2026-04-11）

### 7.1 文献库可用性结论

已核查 `D:\Desktop\paper` 下文献资产，当前主题“`xi` 对介子质量与 Mott 相变的影响”可由现有库直接支撑：

- Bib 数据：`D:\Desktop\paper\bib\PNJL.bib`、`D:\Desktop\paper\bib\动量各向异性.bib`、`D:\Desktop\paper\bib\输运系数.bib`
- PDF 库：`D:\Desktop\paper\docs\papers\PNJL\`、`D:\Desktop\paper\docs\papers\动量各向异性\`、`D:\Desktop\paper\docs\papers\输运系数\`
- Mott 直连条目（已确认可用）：
  - `blaschkeMottDissociationPions2017`
  - `maslovEffectMesonicOffshell2023`

### 7.2 新论文模板目录（已创建）

已在 `D:\Desktop\paper\My_Paper\` 下创建新稿目录：

- `D:\Desktop\paper\My_Paper\pnjl_meson_mott_xi\main.tex`
- `D:\Desktop\paper\My_Paper\pnjl_meson_mott_xi\refs.bib`
- `D:\Desktop\paper\My_Paper\pnjl_meson_mott_xi\BUILD.md`
- `D:\Desktop\paper\My_Paper\pnjl_meson_mott_xi\mainNotes.bib`

并完成以下初始化：

- `main.tex` 标题/摘要/章节占位已切换至“介子质量 + Mott + xi”主题；
- `refs.bib` 已补入 Mott 核心文献键（见 7.1）；
- `BUILD.md` 编译路径已从 `pnjl_aniso_transport` 改为 `pnjl_meson_mott_xi`；
- 已通过 `pdflatex + bibtex + pdflatex + pdflatex` 生成 `main.pdf`。

### 7.3 在 Julia 项目中接入文献库的路径建议

对于“在 `Julia_RelaxTime` 根目录下以链接方式访问 `D:\Desktop\paper`”这一需求，`junction` 在 Windows 上是可行方案。

建议原则：

1. 优先使用 **目录 junction**（`mklink /J`）而不是复制文献库，减少重复与同步成本；
2. 链接名建议中性且稳定，如 `paper_lib`，避免与仓库既有目录冲突；
3. 默认不纳入 git 版本控制（通常可通过 `.gitignore` 忽略链接名）；
4. 文档中统一使用“外部文献库映射路径”表述，避免误判为仓库内源文件。

示例（在 `D:\Desktop\Julia_RelaxTime` 的父级命令行执行）：

```bat
mklink /J D:\Desktop\Julia_RelaxTime\paper_lib D:\Desktop\paper
```

若后续不再需要：

```bat
rmdir D:\Desktop\Julia_RelaxTime\paper_lib
```

注：`rmdir` 仅移除链接本身，不会删除 `D:\Desktop\paper` 实体内容。
