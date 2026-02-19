# PNJL 扫描输出契约（T-μ / T-ρ）

本文档定义 `PNJL.run_tmu_scan` 与 `PNJL.run_trho_scan` 的 CSV 输出契约，作为下游分析、作图与回归脚本的统一依据。

## 版本与范围

- 契约名：`pnjl_scan_csv_v1`
- 适用模块：`src/pnjl/scans/TmuScan.jl`、`src/pnjl/scans/TrhoScan.jl`
- 默认输出：
  - `data/outputs/results/pnjl/scan/tmu/tmu_scan.csv`
  - `data/outputs/results/pnjl/scan/trho/trho_scan.csv`

## 文件级约定

- 编码：UTF-8
- 分隔符：英文逗号 `,`
- 首行必须为表头（header）
- 布尔字段 `converged` 使用 `true/false`
- 数值字段使用十进制浮点文本（可含科学计数法）
- 失败点也写入行记录（`converged=false`），并通过 `message` 标注原因

## 公共字段（两类扫描共享）

以下字段在 T-μ 与 T-ρ 两类扫描中均存在：

- `T_MeV`
- `xi`
- `pressure_fm4`
- `entropy_fm3`
- `energy_fm4`
- `phi_u`
- `phi_d`
- `phi_s`
- `Phi1`
- `Phi2`
- `M_u_MeV`
- `M_d_MeV`
- `M_s_MeV`
- `iterations`
- `residual_norm`
- `converged`
- `message`

## T-μ 扫描字段（`run_tmu_scan`）

`run_tmu_scan` 的完整 header：

```text
T_MeV,mu_MeV,xi,pressure_fm4,rho,entropy_fm3,energy_fm4,phi_u,phi_d,phi_s,Phi1,Phi2,M_u_MeV,M_d_MeV,M_s_MeV,iterations,residual_norm,converged,message
```

其中 `mu_MeV` 与 `rho` 为 T-μ 模式特有字段。

## T-ρ 扫描字段（`run_trho_scan`）

`run_trho_scan` 的完整 header：

```text
T_MeV,rho,xi,mu_u_MeV,mu_d_MeV,mu_s_MeV,mu_avg_MeV,mu_B_MeV,mu_Q_MeV,mu_S_MeV,pressure_fm4,entropy_fm3,energy_fm4,rho_u_fm3,rho_d_fm3,rho_s_fm3,phi_u,phi_d,phi_s,Phi1,Phi2,M_u_MeV,M_d_MeV,M_s_MeV,iterations,residual_norm,converged,message
```

其中 `rho`、`mu_*_MeV`、`mu_B/Q/S_MeV`、`rho_*_fm3` 为 T-ρ 模式特有字段。

## 兼容性规则

- 下游读取方应按“列名”而非“列位置”解析。
- 若将来新增列，默认只允许追加在末尾，保持既有列名不变。
- 若必须变更既有列名或语义，需升级契约版本（例如 `pnjl_scan_csv_v2`）并在文档中并存说明。
