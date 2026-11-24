## 📘 计算弛豫时间时对平均散射率的求和说明

**散射过程来源说明**  
本文档中涉及的散射过程列表见 [ScatteringProcesses_AllPossible.md](../scattering/ScatteringProcesses_AllPossible.md) 中的表1（夸克-夸克散射）、表2（夸克-反夸克散射）与表3（反夸克-反夸克散射，表1 的电荷共轭）。平均散射率的积分定义见 [AverageScatteringRate_FromCrossSection.md](../scattering/AverageScatteringRate_FromCrossSection.md)。以下弛豫时间表达式利用同位旋对称性（\(u\) 和 \(d\) 夸克质量相同）进行简化，并与代码实现保持一致。

在驰豫时间近似（RTA）下，第 \( i \) 类夸克的弛豫时间 \( \tau_i \) 由下式给出：

\[
\tau_i^{-1} = \sum_j \rho_j \, \bar{w}_{ij}
\]

其中：
- \( \rho_j \) 是第 \( j \) 类粒子的数密度；
- \( \bar{w}_{ij} \) 是 \( i \) 和 \( j \) 粒子间的平均散射率，由式 (5.10) 计算。

### 📊 各粒子的弛豫时间表达式

根据论文中给出的表达式（考虑同位旋对称 \( u \leftrightarrow d \)），整理如下：

| 粒子类型 | 弛豫时间 \( \tau_i^{-1} \) 表达式（求和项） |
|----------|----------------------------------------|
| **u夸克** | \( n_u (\bar{w}_{uu \to uu} + \bar{w}_{ud \to ud}) \) <br> + \( n_{\bar{u}} (\bar{w}_{u\bar{u} \to u\bar{u}} + \bar{w}_{u\bar{u} \to d\bar{d}} + \bar{w}_{u\bar{u} \to s\bar{s}} + \bar{w}_{u\bar{d} \to u\bar{d}}) \) <br> + \( n_s \bar{w}_{us \to us} \) <br> + \( n_{\bar{s}} \bar{w}_{u\bar{s} \to u\bar{s}} \) |
| **d夸克** | 同位旋对称下 \(\tau_d^{-1}=\tau_u^{-1}\)；实现中直接返回与 **u夸克** 相同的结果（\(d=u\)）。 |
| **s夸克** | \( 2n_u \bar{w}_{us \to us} \) <br> + \( 2n_{\bar{u}} \bar{w}_{s\bar{u} \to s\bar{u}} \) <br> + \( n_s \bar{w}_{ss \to ss} \) <br> + \( n_{\bar{s}} (\bar{w}_{s\bar{s} \to s\bar{s}} + 2\bar{w}_{s\bar{s} \to u\bar{u}}) \) |
| **反u夸克** | \( n_u (\bar{w}_{u\bar{u} \to u\bar{u}} + \bar{w}_{u\bar{u} \to d\bar{d}} + \bar{w}_{u\bar{u} \to s\bar{s}} + \bar{w}_{d\bar{u} \to d\bar{u}}) \) <br> + \( n_{\bar{u}} (\bar{w}_{\bar{u}\bar{d} \to \bar{u}\bar{d}} + \bar{w}_{\bar{u}\bar{u} \to \bar{u}\bar{u}}) \) <br> + \( n_s \bar{w}_{s\bar{u} \to s\bar{u}} \) <br> + \( n_{\bar{s}} \bar{w}_{\bar{u}\bar{s} \to \bar{u}\bar{s}} \) |
| **反d夸克** | 同位旋对称下 \(\tau_{\bar{d}}^{-1}=\tau_{\bar{u}}^{-1}\)；实现中直接返回与 **反u夸克** 相同的结果（\(\bar{d}=\bar{u}\)）。 |
| **反s夸克** | \( 2n_u \bar{w}_{u\bar{s} \to u\bar{s}} \) <br> + \( 2n_{\bar{u}} \bar{w}_{\bar{u}\bar{s} \to \bar{u}\bar{s}} \) <br> + \( n_{\bar{s}} \bar{w}_{\bar{s}\bar{s} \to \bar{s}\bar{s}} \) <br> + \( n_s (\bar{w}_{s\bar{s} \to s\bar{s}} + 2\bar{w}_{s\bar{s} \to u\bar{u}}) \) |

---

## 🔎 代码实现对照（omega 与 process key）

代码里实际求和发生在 `RelaxationTime.relaxation_rates`，并使用平均散射率字典/NamedTuple 的 key（`Symbol`）取值。为了避免文档与实现漂移，这里给出逐项对应关系：

| 目标 | 实现中的密度因子 | 实现中的 rate key（process） | 对应符号写法 |
|------|------------------|----------------------------|--------------|
| \(\omega_u\) | \(n_u\) | `:uu_to_uu`, `:ud_to_ud` | \(\bar{w}_{uu\to uu}\), \(\bar{w}_{ud\to ud}\) |
| \(\omega_u\) | \(n_{\bar{u}}\) | `:uubar_to_uubar`, `:uubar_to_ddbar`, `:uubar_to_ssbar`, `:udbar_to_udbar` | \(\bar{w}_{u\bar{u}\to u\bar{u}}\), \(\bar{w}_{u\bar{u}\to d\bar{d}}\), \(\bar{w}_{u\bar{u}\to s\bar{s}}\), \(\bar{w}_{u\bar{d}\to u\bar{d}}\) |
| \(\omega_u\) | \(n_s\) | `:us_to_us` | \(\bar{w}_{us\to us}\) |
| \(\omega_u\) | \(n_{\bar{s}}\) | `:usbar_to_usbar` | \(\bar{w}_{u\bar{s}\to u\bar{s}}\) |
| \(\omega_s\) | \(2n_u\) | `:us_to_us` | \(\bar{w}_{us\to us}\)（同位旋下合并 \(ds\to ds\)） |
| \(\omega_s\) | \(2n_{\bar{u}}\) | `:subar_to_subar` | \(\bar{w}_{s\bar{u}\to s\bar{u}}\)（同位旋下合并 \(s\bar{d}\to s\bar{d}\)） |
| \(\omega_s\) | \(n_s\) | `:ss_to_ss` | \(\bar{w}_{ss\to ss}\) |
| \(\omega_s\) | \(n_{\bar{s}}\) | `:ssbar_to_ssbar`, `:ssbar_to_uubar` | \(\bar{w}_{s\bar{s}\to s\bar{s}}\), \(\bar{w}_{s\bar{s}\to u\bar{u}}\)（并代表 \(d\bar{d}\)，所以有系数 2） |
| \(\omega_{\bar{u}}\) | \(n_u\) | `:uubar_to_uubar`, `:uubar_to_ddbar`, `:uubar_to_ssbar`, `:dubar_to_dubar` | \(\bar{w}_{u\bar{u}\to u\bar{u}}\), \(\bar{w}_{u\bar{u}\to d\bar{d}}\), \(\bar{w}_{u\bar{u}\to s\bar{s}}\), \(\bar{w}_{d\bar{u}\to d\bar{u}}\) |
| \(\omega_{\bar{u}}\) | \(n_{\bar{u}}\) | `:ubardbar_to_ubardbar`, `:ubarubar_to_ubarubar` | \(\bar{w}_{\bar{u}\bar{d}\to \bar{u}\bar{d}}\), \(\bar{w}_{\bar{u}\bar{u}\to \bar{u}\bar{u}}\) |
| \(\omega_{\bar{u}}\) | \(n_s\) | `:subar_to_subar` | \(\bar{w}_{s\bar{u}\to s\bar{u}}\) |
| \(\omega_{\bar{u}}\) | \(n_{\bar{s}}\) | `:ubarsbar_to_ubarsbar` | \(\bar{w}_{\bar{u}\bar{s}\to \bar{u}\bar{s}}\) |
| \(\omega_{\bar{s}}\) | \(2n_u\) | `:usbar_to_usbar` | \(\bar{w}_{u\bar{s}\to u\bar{s}}\)（同位旋下合并 \(d\bar{s}\to d\bar{s}\)） |
| \(\omega_{\bar{s}}\) | \(2n_{\bar{u}}\) | `:ubarsbar_to_ubarsbar` | \(\bar{w}_{\bar{u}\bar{s}\to \bar{u}\bar{s}}\)（同位旋下合并 \(\bar{d}\bar{s}\to \bar{d}\bar{s}\)） |
| \(\omega_{\bar{s}}\) | \(n_{\bar{s}}\) | `:sbarsbar_to_sbarsbar` | \(\bar{w}_{\bar{s}\bar{s}\to \bar{s}\bar{s}}\) |
| \(\omega_{\bar{s}}\) | \(n_s\) | `:ssbar_to_ssbar`, `:ssbar_to_uubar` | \(\bar{w}_{s\bar{s}\to s\bar{s}}\), \(\bar{w}_{s\bar{s}\to u\bar{u}}\)（并代表 \(d\bar{d}\)，所以有系数 2） |

## 📌 说明

1. **符号约定**：
   - \( n_i \)：第 \( i \) 类粒子的数密度；
   - \( \bar{w}_{ij \to cd} \)：表示 \( i + j \to c + d \) 散射过程的平均散射率；
   - 上标“¯”表示反粒子。

2. **散射过程类型**：
   - **夸克-夸克散射**（如 \( uu \to uu \)）；
   - **夸克-反夸克散射**（如 \( u\bar{u} \to u\bar{u} \)）；
   - **夸克-不同味夸克散射**（如 \( us \to us \)）；
   - **夸克-反不同味夸克散射**（如 \( u\bar{s} \to u\bar{s} \)）。

3. **同位旋对称性**：
   - 本项目默认采用同位旋对称近似：\(n_u=n_d\)、\(n_{\bar{u}}=n_{\bar{d}}\)，且相应散射率满足同位旋等价关系。
   - 因此实现中直接返回 \(\tau_d=\tau_u\)、\(\tau_{\bar{d}}=\tau_{\bar{u}}\)（不再重复计算）。

4. **为什么会出现系数 \(2\)**：
   - **\(2n_u\,\bar{w}_{us\to us}\)**：表示 \(s\) 与 \(u\) 和 \(d\) 两种轻夸克的散射贡献相同（\(\bar{w}_{us\to us}=\bar{w}_{ds\to ds}\) 且 \(n_d=n_u\)），因此 \(n_u\bar{w}_{us\to us} + n_d\bar{w}_{ds\to ds} = 2n_u\bar{w}_{us\to us}\)。
   - **\(2n_{\bar{u}}\,\bar{w}_{s\bar{u}\to s\bar{u}}\)** 与 **\(2n_{\bar{u}}\,\bar{w}_{\bar{u}\bar{s}\to \bar{u}\bar{s}}\)**：同理，\(\bar{u}\) 与 \(\bar{d}\) 对 \(s\) 或 \(\bar{s}\) 的散射在同位旋下等价并且密度相同。
   - **\(2\bar{w}_{s\bar{s}\to u\bar{u}}\)**：\(s\bar{s}\) 湮灭/产生轻味对时，\(u\bar{u}\) 与 \(d\bar{d}\) 两个末态在同位旋下等价且贡献相同，因此用一个过程的平均速率代表两者并乘以 2。
   - 反过来，像 **\(\bar{u}\)** 的 \(n_u(\cdots + \bar{w}_{d\bar{u}\to d\bar{u}})\) 这种写法，是因为 \(\bar{u}\) 与轻夸克散射时允许的独立过程并不完全相同（例如 \(u\bar{u}\) 与 \(d\bar{u}\)），因此不能简单合并成一个 \(2\times\) 系数，需要显式保留电荷共轭/不同味过程。

4. **涉及的介子传播子**：
   - 各散射过程的散射截面 \( \sigma_{ij \to cd} \) 需通过介子传播子计算，具体形式见论文第5.1.4节。

5. **与独立散射过程表格的对应关系**：
   - \(\bar{u}\) 与 \(\bar{s}\) 的弛豫时间中会出现反夸克-反夸克过程，因此需同时参考表3。
   - 表2中额外出现的 \(d\bar{u}\to d\bar{u}\) 与 \(s\bar{u}\to s\bar{u}\) 属于电荷共轭/换味补充过程，用来与 Fortran 的过程集合保持一致。

6. **数值计算注意事项**：
   - 在实现中通常直接返回 \(u=d\)、\(\bar{u}=\bar{d}\) 的结果（同位旋对称）。
---

