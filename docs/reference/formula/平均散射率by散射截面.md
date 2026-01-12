在球坐标系下，考虑轴对称各向异性分布函数 $f_i(p_i, \cos\theta_i)$ 和 $f_j(p_j, \cos\theta_j)$，其中 $\theta_i$ 和 $\theta_j$ 为动量与某一固定轴（如 $z$ 轴）的夹角。通过对两个方位角 $\phi_i$ 和 $\phi_j$ 的积分进行简化，利用被积函数仅依赖于其差值 $\phi = \phi_i - \phi_j$ 的性质，可得平均散射率的简化表达式为：

$$
\omega_{ij} = \frac{d_q^2}{\rho_i \rho_j (2\pi)^5} \int_0^\Lambda p_i^2 dp_i \int_0^\Lambda p_j^2 dp_j \int_{-1}^1 d\cos\theta_i \int_{-1}^1 d\cos\theta_j \int_0^{2\pi} d\phi \ f_i(p_i,\cos\theta_i) f_j(p_j,\cos\theta_j) \ v_{\text{rel}} \ \sigma_{ij\to cd}(s,T,\mu_q)
$$

其中：
- $d_q = 6$ 为简并度（PNJL 模型）。
- $\Lambda$ 为动量模的截断。
- $s$ 为 Mandelstam 变量，$s = (p_i + p_j)^2 = m_i^2 + m_j^2 + 2E_iE_j - 2p_i p_j \cos\Theta$，其中 $\cos\Theta = \cos\theta_i \cos\theta_j + \sin\theta_i \sin\theta_j \cos\phi$。
- **末态阈值提醒**：对给定过程 $ij\to cd$，即使 $s$ 由初态四动量计算得到，也仍需满足末态产生阈值
  $$
  s \ge (m_c+m_d)^2
  $$
  因此实践中可用一个统一的阈值条件
  $$
  s > \max\big((m_i+m_j)^2,\,(m_c+m_d)^2\big)
  $$
  来裁剪积分点（对如 $s\bar s\to u\bar u$ 这类“初态更重”的过程尤为关键）。

  - **有限截断导致的上界提醒**：当且仅当平均散射率分子动量积分采用硬截断 $p_i,p_j\in[0,\Lambda]$ 时，可达到的 Mandelstam 变量 $s$ 也被限制在一个有限区间内。记
    $$
    E_k^{\max}=\sqrt{m_k^2+\Lambda^2}\quad (k=i,j,c,d)
    $$
    则一个常用的可行上界写法为
    $$
    s \le s_{up}=\min\Big((E_i^{\max}+E_j^{\max})^2,\,(E_c^{\max}+E_d^{\max})^2\Big)
    $$
    物理含义是：在 $p\le\Lambda$ 的约束下，无论初态还是末态都不可能提供超过其各自最大总能量对应的质心系不变量。
    
    （若数值实现中对 $\sigma(s)$ 做了插值/外推，显式施加 $s\in[s_{bo},s_{up}]$ 的裁剪也能避免在“不可达的 $s$ 区域”取值带来的伪贡献；如果分子动量积分改为半无穷区间，则不需要引入此 $s_{up}$。）
  
  （数值实现里有时会在右侧再乘一个很小的安全系数 $1+\epsilon$ 来避免浮点数恰落在边界导致的开方/插值问题；如果积分节点不会精确打到边界，通常可取 $\epsilon=0$。）
- **注意**：在计算相对速度 $v_{\text{rel}}$ 时，应使用质心系下的能量 $E_i^*$ 和 $E_j^*$，而非实验室系下的能量。具体计算公式为（参考 PDF 第 4 页公式 5.11）：
  $$
  E_i^* = \frac{s + m_i^2 - m_j^2}{2\sqrt{s}}, \quad
  E_j^* = \frac{s - m_i^2 + m_j^2}{2\sqrt{s}}
  $$
- 相对速度 $v_{\text{rel}}$ 的计算公式为（参考 PDF 第 4 页公式 5.11）：
  $$
  v_{\text{rel}} = \frac{\sqrt{(E_i^* E_j^* - |\mathbf{p}_i^*||\mathbf{p}_j^*|)^2 - (m_i m_j)^2}}{E_i^* E_j^*}
  $$
  其中 $|\mathbf{p}_i^*|$ 和 $|\mathbf{p}_j^*|$ 为质心系下的动量大小，具体表达式见 PDF。
- 数密度 $\rho_i$ 和 $\rho_j$ 由下式给出：
  $$
  \rho_i = d_q \frac{1}{(2\pi)^2} \int_0^\infty p_i^2 dp_i \int_{-1}^1 d\cos\theta_i \ f_i(p_i,\cos\theta_i)
  $$
  类似地计算 $\rho_j$。

**重要说明**：数密度积分应使用半无穷积分范围 $[0, \infty)$，而非有限截断 $[0, \Lambda]$。这是因为数密度是物理可观测量，不应依赖于模型的动量截断参数。在数值实现中，可通过变量替换 $p = \text{scale} \cdot t / (1-t)$ 将半无穷积分映射到有限区间 $[0, 1)$。

该表达式已将两个方位角积分简化为一个方位角 $\phi$ 的积分，保留了必要的相对方向依赖性。

---

## $\sigma(s)$ 的快变区与插值点选取原则

当我们用“预计算 $\sigma(s)$ + 插值”来加速平均散射率积分时，插值误差的大小确实强烈依赖于 $\sigma(s)$ 随 $s$ 的变化结构。对 PNJL/介子传播子主导的散射，$\sigma(s)$ 的**快变通常来自以下几类机制**：

1. **动力学阈值附近（最重要且最普遍）**
   - 起始阈值：$s\to s_{bo}^+=\max\big((m_i+m_j)^2,(m_c+m_d)^2\big)$。
   - 在阈值附近，相空间/质心动量 $p^*(s)$ 往往呈现 $p^*(s)\propto\sqrt{s-s_{bo}}$ 的非解析行为，导致 $\sigma(s)$ 在一个很窄的 $s$ 区间内快速“爬升/转折”（即便矩阵元本身较平滑）。

2. **传播子共振/近极点结构（可能非常尖锐）**
   - 若散射振幅里出现形如
     $$
     \frac{1}{\big(s-m_M^2\big)^2+\big(m_M\Gamma_M\big)^2}
     $$
     的 Breit–Wigner 型结构（或由介子传播子/极化函数导致的等效结构），则在 $s\approx m_M^2$ 附近会出现尖峰或强烈色散。
   - 快变尺度通常由宽度控制：$|s-m_M^2|\lesssim m_M\Gamma_M$。在高温高密下 $\Gamma_M$ 可能变大，峰会变“宽”；在某些参数区间也可能出现相对尖锐的结构。

3. **环积分/极化函数的分支点与“拐点/尖点”（cusp）**
   - 一环积分与极化函数常含对数与平方根结构，可能在某些特征点附近出现导数不连续或快速变化（例如与 $\sqrt{s}$ 相关的阈值/切割）。
   - 这类结构不一定给出非常尖的峰，但会让高阶插值更容易产生振荡（Runge 现象）。

基于上述结构，一个“物理引导 + 数值自适应”的 $s$ 取点策略通常更稳健：

- **把必需的“断点”显式放进 $s$ 网格**
  - 至少包含：$s_{bo}$。
  - 若你能从模型里识别出可能的共振位置（相关介子质量 $m_M$），建议把 $s=m_M^2$（以及其左右一定范围）作为加密中心。

- **对阈值附近做非均匀加密（变量变换）**
  - 常用做法是在阈值上方使用类似
    $$
    s = s_{bo} + \alpha\,x^2\quad (x\in[0,1])
    $$
    或 $s=s_{bo}+\alpha\,x^{\gamma}$（$\gamma>1$）来把点“挤”到 $s_{bo}$ 附近。
  - 经验上，若你的误差主要来自阈值，先加密 $s_{bo}$ 附近往往比在高 $s$ 区域均匀加点更有效。

- **对共振附近做局部加密**
  - 若存在峰，建议在 $s\in[m_M^2-k\,m_M\Gamma_M,\,m_M^2+k\,m_M\Gamma_M]$（例如 $k=3\sim 5$）内用更密的点。
  - 当 $\Gamma_M$ 难以事先精确估计时，可以先用粗网格扫描 $\sigma(s)$，用“峰值位置/曲率”来自动识别需要加密的区间。

- **用自适应细分来“按需加点”（推荐）**
  - 先用一个较粗的 $s$ 网格计算 $\sigma(s)$，然后对每个区间 $[s_k,s_{k+1}]$ 用中点 $s_{mid}$（或多点）重新计算 $\sigma(s_{mid})$，并与插值得到的 $\tilde\sigma(s_{mid})$ 比较。
  - 若满足
    $$
    \frac{|\sigma(s_{mid})-\tilde\sigma(s_{mid})|}{\max(\sigma(s_{mid}),\,\sigma_0)} > \varepsilon_{\text{interp}}
    $$
    就把该区间二分（或按比例细分），迭代直到所有区间都满足容差。
  - 这里 $\sigma_0$ 是防止 $\sigma\approx 0$ 时相对误差失真的小正数。

- **插值形式上的小建议（提升稳健性）**
  - 若 $\sigma(s)$ 跨越数量级很大，可以考虑对 $\log(\sigma+\sigma_0)$ 插值，再指数还原，以更好控制相对误差。
  - 选用“形状保持”的插值（如 PCHIP/分段单调三次）可减少过冲与负值；若用样条插值，建议额外做非负截断（$\sigma\leftarrow\max(\sigma,0)$）。

以上原则的直观目标是：**让插值误差显著小于积分（GL）误差**，这样“GL5D + $\sigma(s)$ cache 插值”在同等运行时间下才会稳定地胜出。


---

## 当前数值实现方案（已敲定）：w0cdf 取点 + PCHIP 插值

在本项目的“预计算 $\sigma(s)$ + 插值”加速路线中，**不再继续对比不同取点策略**，统一采用：

- **取点（点位分布）**：`w0cdf`
- **插值（cache 查询）**：`pchip`（PCHIP，形状保持的分段三次）
- **默认点数**：`N=60`
- **默认 5D 求积节点**（用于 $\omega$ 分子积分）：`p=20, angle=4, phi=8`
- **w0cdf 取点用的低成本网格**（用于构造 $w_0(\sqrt{s})$ 的 CDF）：`p=14, angle=4, phi=8`

这里的 `w0cdf` 指：先用低成本的 GL5D 扫描得到按 $\sqrt{s}$ 分箱的权重 $w_0(\sqrt{s})$（不含截面），再按其累计分布函数（CDF）做**等质量、确定性**的采样，从而生成每个过程的 $\sqrt{s}$ 取点表。

### 为什么选 w0cdf（而不是基于 $w_\sigma$ 的 wsgrid 等）

核心考虑是“鲁棒性 + 可迁移性 + 确定性”：

1. **$w_0$ 不含 $\sigma(s)$，对参数变化更稳定**
  - $w_\sigma$（含截面）会把“截面结构”直接反馈到取点，导致不同过程/不同参数点下核心区间可能发生剧烈移动；当 $\sigma(s)$ 很小或数值噪声占主导时，基于 $w_\sigma$ 的取点更容易不稳定。
  - 采用 $w_0$ 的优点是把“相空间与分布函数带来的主要贡献区间”作为主导，避免被某些尖峰/近零区间的 $\sigma$ 牵着走。

2. **w0cdf 是确定性的，且实现简单**
  - 给定同一份 $w_0(\sqrt{s})$ 分布与点数 N，生成的网格完全确定；方便复现实验与生产部署。

3. **配合 PCHIP，更不易产生非物理振荡**
  - `pchip`（形状保持的分段三次）在处理阈值附近的快变、以及跨数量级变化的区域时，通常比普通三次样条更稳健；同时相对 `linear` 往往能以相同 N 显著降低误差。

### 评估口径（用于回归与随机点验证）

为了衡量“$\sigma(s)$ 取点/插值”对平均散射率的影响，我们使用脚本 `scan_omega_s_importance.jl` 输出的按 $\sqrt{s}$ 分箱权重，并比较每个过程的**总重要性**（默认用 $\sum w_\sigma$，若 $\sum w_\sigma = 0$ 则回退到 $\sum w_0$）。

- 由于数值噪声可能导致单箱权重出现负值，评估与网格设计统一采用
  $$
  w\leftarrow\max(w,0)
  $$
  后再求和，避免正负抵消造成相对误差失真。

（随机参数点上的对比结果与统计摘要见后续更新条目。）

### 随机参数点验证（2026-01-10）

为验证 `w0cdf + (pchip/linear)` 在参数空间内的鲁棒性，我们在
$T\in(10,350)$ MeV、$\mu_B\in(0,1500)$ MeV、$\xi\in(0,1)$ 上随机抽取 **8 个点**（seed=20260110），并对比：

- 点数：`N=120` 与 `N=240`
- 插值：`linear` 与 `pchip`
- 参考口径：用相同参数点下的 dense σ-cache（仅重算 3 个过程）生成参考 importance，再用“synthetic sparse cache”复现生产模式下的插值误差。

结果汇总 CSV：
- [data/processed/results/relaxtime/w0cdf_random_params_N120_240_summary.csv](data/processed/results/relaxtime/w0cdf_random_params_N120_240_summary.csv)

按 8 个点统计（指标为 `max_abs_rel`；越小越好）：

| N | interp | median(max_abs_rel) | max(max_abs_rel) |
|---:|---|---:|---:|
| 120 | linear | 4.71e-3 | 0.875 |
| 120 | pchip  | 4.11e-3 | 0.872 |
| 240 | linear | 1.32e-3 | 0.887 |
| 240 | pchip  | 9.73e-4 | 0.886 |

备注：绝大多数随机点上误差处于 $10^{-3}\sim 10^{-2}$ 量级；但存在少数点出现显著 outlier（例如某些过程在该参数点下 $\sigma(s)$ 接近全零导致 $w_\sigma$ 缺失/回退到 $w_0$，从而 max 指标被单点主导）。这种 outlier 更像是“该过程在该参数点下几乎不贡献/数值噪声主导”的信号，建议单独排查并在生产中用额外的稳定化策略处理（例如提高该点的 σ-cache 分辨率或为极小 σ 区域设置更明确的截断/保护）。


---

## 基于 $w_{\sigma}$ 的“重要 $\sqrt{s}$ 区域”与 17 个散射过程差异（经验观察）

上面的小节强调了“$\sigma(s)$ 本身在哪里快变”，但在平均散射率积分里，更关键的是：**哪些 $\sqrt{s}$ 区域真正贡献了主要权重**。

为此我们用一个低成本的 GL5D 网格，把积分写成对离散积分点的求和，并按 $\sqrt{s}$ 分箱累计两类权重：

- $w_0(\sqrt{s})$：不含截面（只含 $p_i^2p_j^2 f_if_j v_{rel}$ 与测度等因子）的权重；
- $w_{\sigma}(\sqrt{s})$：包含 $\sigma(s)$（来自 cache/插值）的权重。

在“按 $w_{\sigma}$ 质量集中区加密”时，我们用 $w_{\sigma}$ 的分位数（例如 $q_{0.8}$、$q_{0.95}$）来定义核心贡献区间。下面记录在 $T=150\,\mathrm{MeV},\ \mu_B=800\,\mathrm{MeV},\ \xi=0$ 时，对**全部 17 个过程**得到的显著差别（来自脚本 `scan_omega_s_importance.jl` 的 all17 输出；网格参数为 $p=8,\cos=3,\phi=6$，分箱 `bins=60`，$\sigma$ 来自 `xi0_ns1000_nt6_all` cache）。

### 1) 最显著的差别：不同过程的 $w_{\sigma}$ 核心区间位置不同

用 $[q_{0.8},\ q_{0.95}]$ 表示“80%–95% 贡献核心带”（单位 MeV）：

| 过程 | 80%–95% 核心带 $[q_{0.8},q_{0.95}]$ | 备注 |
|---|---:|---|
| `uu_to_uu` | $[779.6,\ 919.3]$ | light-light 中“更硬”（更高 $\sqrt{s}$） |
| `ud_to_ud` | $[773.4,\ 913.8]$ | 与 `uu_to_uu` 非常接近 |
| `ubarubar_to_ubarubar` | $[618.5,\ 808.8]$ | light-light 但整体更低 |
| `ubardbar_to_ubardbar` | $[601.9,\ 801.1]$ | 与 `ubarubar_to_ubarubar` 接近 |
| `uubar_to_uubar` | $[533.6,\ 741.4]$ | 分布偏低且长尾更明显（q50 更低） |
| `uubar_to_ddbar` | $[435.7,\ 633.8]$ | light-light 中偏低 |
| `udbar_to_udbar` | $[280.2,\ 579.0]$ | 更靠近阈值侧；与 `dubar_to_dubar` 完全一致 |
| `dubar_to_dubar` | $[280.2,\ 579.0]$ | 电荷共轭一致性验证 |
| `us_to_us` | $[945.0,\ 1082.2]$ | light-strange 中最“硬” |
| `usbar_to_usbar` | $[832.3,\ 986.6]$ | light-strange 中等 |
| `ubarsbar_to_ubarsbar` | $[889.9,\ 1043.4]$ | light-strange 偏高 |
| `subar_to_subar` | $[752.2,\ 890.9]$ | light-strange 偏低 |
| `uubar_to_ssbar` | $[1023.9,\ 1066.3]$ | strange-strange 中**非常窄**，集中度最高 |
| `ssbar_to_uubar` | $[1089.8,\ 1170.6]$ | strange-strange 偏低 |
| `ssbar_to_ssbar` | $[1114.2,\ 1238.0]$ | strange-strange 中等 |
| `ss_to_ss` | $[1173.0,\ 1278.3]$ | strange-strange 中最“硬” |
| `sbarsbar_to_sbarsbar` | $[1171.5,\ 1277.1]$ | 与 `ss_to_ss` 非常接近 |

这些差异说明：**“同一个 $\sqrt{s}$ 取点策略”并不适合所有过程**。即便在同一热力学点、同一网格与同一 $\sigma$-cache 下，不同过程的主要贡献带也会显著错开。

### 2) 阈值附近快变并不等于“对 $\omega$ 重要”

虽然 $\sigma(s)$ 往往在阈值附近变化最快，但从 $w_{\sigma}$ 的分布来看，许多过程的主要贡献并不在阈值极近处，而是在更高的 $\sqrt{s}$ 带宽里（尤其是 light-light 的 `uu_to_uu`/`ud_to_ud` 与 strange-strange 的 `ss_to_ss`）。因此，**加密点位应优先跟随 $w_{\sigma}$ 的质量集中区**，而不是仅仅对阈值做超密加点。

### 3) 用于 $\sigma(s)$ cache 取点的直接建议

- 每个过程至少覆盖可达区间 $[\sqrt{s}_{bo},\ \sqrt{s}_{up}]$（由有限动量截断给出）。
- 把大部分点数分配到该过程的核心带（例如 70% 点落在 $[q_{0.8},q_{0.95}]$），两端各留少量点用于覆盖尾部。
- 仍需保留少量阈值上方点用于数值稳定与插值边界控制，但不必把绝大多数点都堆在 $\sqrt{s}_{bo}$ 附近。

（后续若要进一步精炼，可在核心带内再结合 $\sigma(s)$ 的局部曲率/误差估计做自适应细分，实现“物理权重 + 插值误差”的双重驱动加密。）


---

**注意**：文档中先前给出的将角度积分简化到 $\cos\theta \in [0,1]$、$\phi \in [0,\pi]$ 的做法不正确。原因在于

$$
\cos\Theta = \cos\theta_i\cos\theta_j + \sin\theta_i\sin\theta_j\cos\phi
$$

并不在联合变换 $(\cos\theta_i,\cos\theta_j,\phi) \to (-\cos\theta_i,-\cos\theta_j,-\phi)$ 下保持不变，因此无法通过对称性将积分区域简单折半。数值实现必须保留原始积分范围：$\cos\theta_i,\cos\theta_j \in [-1,1]$ 和 $\phi \in [0,2\pi]$。

数密度的常数因子可等价写为：

$$
\rho_i = d_q \frac{1}{(2\pi)^2} \int_0^{\infty} p_i^2 dp_i \int_{-1}^1 d\cos\theta_i\ f_i(p_i,\cos\theta_i)
= \frac{d_q}{2\pi^2} \int_0^{\infty} p_i^2 dp_i \int_{0}^1 d\cos\theta_i\ f_i(p_i,\cos\theta_i)
$$

（已移除）此前文档中提到的“3D 各向同性特化/Fortran 风格参考实现”相关内容容易造成混淆，当前实现统一走通用公式路径。

