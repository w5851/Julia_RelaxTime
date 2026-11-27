# 桌面旧程序备份分诊记录（2026-11-27）

范围：`D:\Desktop\程序备份` 原 7 个压缩包（5.8 MB，2025-09-09 至 2025-11-29）。
本文件记录删除决定、保留项的内容事实，以及是否吸纳进本项目的判断依据。

## 删除（4 个，3.2 MB）

| 文件 | 内容 | 依据 |
| --- | --- | --- |
| `Julia.zip` | `PNJL_Physics_Simulation` + `Rotation_PNJL`，含 `.git` | 本仓库的前身（同样的 `agent/`、`docs/api_*`、`pnjl_aniso_refactor` 报告）。rotation 现由 `src/models/variants/rotation` 承担 |
| `Console2.zip` | Visual Fortran 工程骨架 + Debug 产物 | 见下：另两个 Fortran 包的共同母本，本身被 `src/models/pnjl_physics` 取代 |
| `PNJL磁场tmu.zip` | Fortran 磁场变体 | magnetic 路线由 `src/models/pnjl_magnetic` 承担 |
| `magnetic_pnjl.zip` | Julia 磁场版（8 个 `.jl` + 相图 CSV） | 文件集与 `PNJL_tasalis` 同骨架；磁场部分已取代 |

`Console2`、`PNJL磁场tmu`、`tsallis磁场PNJL t_mu` 三个包的 `.f90` 前五个文件完全同名
（`AA Main project`、`BB constants`、`BB equation`、`BB quantity`、`BB setting constants`），
是同一个 Visual Fortran 工程的三个变体，不是三份独立程序。

## 保留（3 个，2.6 MB）

### `tsallis磁场PNJL t_mu.zip`（Fortran）

真实源码 14 个 `.f90`（另有 Debug 下的 `__genmod` 产物）。是**唯一完整的 Tsallis + 磁场实现**：
`q = 1.05d0`，`exp_q`/`log_q` 与 q 变形的 Polyakov 因子 `AA`/`AA_bar` 出现在
`BB equation.f90` 的三种味道分布中，含 Landau 能级求和。

### `PNJL_tasalis.zip`（Julia）

上面那份的 Julia 移植（8 个 `.jl` + 2 个 `.dat`）。**读码后发现实质缺陷，不能直接吸纳**：

- `exp_q(x, q)` 用 `if p<0` 判断一个**数组**（Julia 中直接抛错），且两个分支代数上完全等价，
  即 Tsallis 的截断语义从未实现；`log_q` 同样两支等价；
- `main_rho.jl` 里 `-dOmgea_dmu_B / rho0` 调用了函数但没传参；
- 同时存在 `Omega` 与 `Omega2` 两套不一致的自由能构造（`Omega2` 用 `fphi`/`fphibar`，
  `Omega` 用 `log_q(AA)`），未注释取舍；
- `get_nodes(128)` 把动量截断硬编码为 `30.0 fm^-1`，与现行热平均采样 profile 体系不兼容。

### `PNJL_AXION.zip`（Python/JAX）

**包名有误导性：其中没有任何 axion 物理**。全包 `axion` 字面命中 0；`theta` 是**各向异性参数**，
喂给 `Energy(sigma, eta, theta, …)`。实际内容是一套 JAX 加速的各向异性 PNJL，
带 `find_cep`、spinodal 与 `Trho`/`Tmu` 扫描，与本项目现行 PNJL 主线功能重叠。

## 吸纳判断

**三个包都不作为代码导入。** 理由分三类：

1. **`PNJL_AXION`**：功能重叠且没有 axion 物理。若将来关心的是 JAX 梯度/并行带来的速度，
   那是工程借鉴问题，与是否吸纳这套代码无关。
2. **Tsallis**：这是现行仓库的**真实能力空白**（`src`/`config`/`tests`/`docs` 中 `tsallis`
   命中 0，git 历史从未出现；论文仓仅 `bib/*.bib` 有条目，正文未写）。但吸纳不是复制文件，
   而是一项研究工作，至少要先定：
   - q 的约定（一致/不一致 Tsallis，以及 `AA^-q` 这个 Polyakov 因子的来处是否与所选约定自洽）；
   - 参数面（备份里两份实现的 `hc` 都是 `197.33`，与现行 `config/physics/default.toml` 的
     `hbarc = 197.3269804` 不一致，直接搬会引入第三套常数）；
   - 载体（须按 `src/models/` 的 entrypoint + config + provenance + 测试层重建，而不是沿用
     脚本式入口）。

   若决定做，建议以 Fortran 版为物理参照（它是原版，Julia 版带上述缺陷），重新实现而非移植。

## 常数对照（移植风险）

| 量 | 备份中的取值 | 现行生产取值 |
| --- | --- | --- |
| `hc` / `hbarc` | `197.33` | `197.3269804` |
| `T0`（Tsallis 包） | `210/hc` | 见 `config/physics/default.toml` |
| `T0`（JAX 包） | `187/hc` | 同上 |

两份备份自己就不一致（210 对 187），说明它们是不同时期、不同参数面的探索，不能当作
一套既定约定。
