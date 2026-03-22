# Legacy Fortran 介子量导出规范（Mott 验证对接）

## 1. 目的

为 Julia_RelaxTime 的 `Mott` 相变温度与介子质量验证建立跨语言可追溯输入。

本规范定义：

- 老程序（Fortran）最小导出字段与单位；
- 文件命名、分隔符、缺失值约定；
- 最小样例点与验收规则；
- 与 Julia 侧映射/测试的对接方式。

## 2. 适用范围与非目标

### 2.1 适用范围

- 语言：Fortran legacy 程序。
- 物理范围：本期先支持 `mu_B=0`，并覆盖 `xi=0`（各向同性）与 `xi` 扫描（各向异性）。
- 通道范围（按理论可计算能力）：
  - 赝标量：`pi`, `K`, `eta`, `eta_prime`
  - 标量：`sigma_pi`, `sigma_K`, `sigma`, `sigma_prime`

### 2.2 非目标

- 本规范不要求老程序重构求解器。
- 本规范不定义图表样式。
- 本规范不覆盖 `mu_B != 0` 的全参数扫描。

## 3. 参考路径（当前快照）

### 3.1 Fortran

- 项目根：`D:/Desktop/Fortran/NJL-ReTime/PNJL-shear-mu-T-method2/Meson-mass-0703`
- 相关源码：
  - `D:/Desktop/Fortran/NJL-ReTime/PNJL-shear-mu-T-method2/Meson-mass-0703/propagator.f90`
  - `D:/Desktop/Fortran/NJL-ReTime/PNJL-shear-mu-T-method2/Meson-mass-0703/amplitude2.f90`
  - `D:/Desktop/Fortran/NJL-ReTime/PNJL-shear-mu-T-method2/Meson-mass-0703/polarization.f90`
  - `D:/Desktop/Fortran/NJL-ReTime/PNJL-shear-mu-T-method2/Meson-mass-0703/main.f90`

## 4. 导出文件规范

### 4.1 文件格式

- 编码：UTF-8
- 分隔符：逗号 `,`
- 小数：科学计数法或常规十进制均可
- 首行：必须是 header
- 缺失值：使用 `NaN`（不要留空字符串）

### 4.2 文件命名

建议命名：

- `legacy_meson_scan_fortran_muB0_v1.csv`

若包含各向异性扫描可加后缀：

- `legacy_meson_scan_fortran_muB0_xi_v1.csv`

## 5. 字段定义（强制）

每行表示一个 `(source_impl, T_MeV, muB_MeV, xi, meson)` 的计算结果。

### 5.1 必选字段

- `record_id`: 字符串唯一 id
- `source_impl`: `fortran`
- `T_MeV`: 温度（MeV）
- `muB_MeV`: 重子化学势（MeV）
- `xi`: 各向异性参数
- `meson`: 介子名（`pi`, `K`, `eta`, `eta_prime`, `sigma_pi`, `sigma_K`, `sigma`, `sigma_prime`）
- `mass_MeV`: 介子质量（MeV）
- `threshold_MeV`: 阈值质量（MeV）
- `gap_MeV`: `mass_MeV - threshold_MeV`
- `mott_flag`: `0`/`1`（是否满足 Mott 条件）
- `solver_status`: `success`/`fallback`/`fail`

### 5.2 建议字段

- `gamma_MeV`: 介子宽度（MeV）
- `residual_norm`: 求解残差
- `notes`: 备注（如分支切换、初值回退）

## 6. Mott 判定口径（导出侧）

- 主判据：`abs(gap_MeV) <= mott_gap_tol_MeV`
- 默认容差建议：`mott_gap_tol_MeV = 1.0`
- 对应写法：
  - `mott_flag = 1` 当 `abs(gap_MeV) <= 1.0`
  - 否则 `mott_flag = 0`

说明：Julia 侧将使用统一后处理提取 `T_Mott`（线性插值主路径 + 二分精化可选），导出侧无需直接输出 `T_Mott`。

## 7. 最小样例点（首轮）

### 7.1 各向同性样例（`xi=0`）

- `muB_MeV = 0`
- `T_MeV ∈ {120, 140, 160, 180, 200, 220, 240, 260}`

### 7.2 各向异性样例（`muB=0`）

- `T_MeV ∈ {160, 180, 200, 220}`
- `xi ∈ [-0.4, 0.4], step=0.05`

> 兼容性说明：对于不支持各向异性内核的 legacy 程序，`xi` 可作为“导出标签维度”用于对齐 Julia 扫描网格；此时同一 `(T, meson)` 在不同 `xi` 下的 `mass/threshold/gap` 允许相同，需在 `solver_status` 或文档备注中明确该口径。

## 8. 验收规则

- 文件可被 Julia CSV 解析器完整读取。
- 必选字段无缺失。
- `gap_MeV` 与 `mass_MeV - threshold_MeV` 一致（误差不超过 `1e-6`）。
- `mott_flag` 与容差判据一致。
- 每个目标 `meson` 在样例点均至少有一行有效记录。

## 9. 与 Julia 侧对接

本规范对应 Julia_RelaxTime 下游任务：

- 映射层：`scripts/analysis/mott_reference_mapping.jl`（待实现）
- 验证测试：`tests/validation/relaxtime/test_mott_meson_validation_isotropic.jl`（待实现）
- 扫描脚本：`scripts/analysis/scan_mott_meson_vs_xi_mu0.jl`（待实现）

## 10. 风险与回退

- 若 Fortran 暂无法同时导出全部通道：
  - 先导出已可稳定求解的通道；
  - 在 `solver_status=fail` 的行保留上下文，而不是删点。
- 若两端单位实现不一致：
  - 先冻结单位换算表，再进入误差比较。

## 11. 执行人检查清单

- [ ] Fortran 导出文件符合字段与单位规范。
- [ ] 目标通道集合完整。
- [ ] `muB=0`、`xi` 样例网格覆盖完成。
- [ ] Julia 映射脚本可直接读取并产出标准结构。

## 12. 仓库内模板文件（已生成）

为减少老程序侧手工拼字段出错，仓库内已生成一份 schema+网格模板：

- `tests/validation/data/targets/relaxtime/legacy/meson/legacy_meson_scan_fortran_muB0_v1.csv`

模板特征：

- 已包含本规范要求的全部必选字段；
- 已覆盖本规范定义的首轮 `T` 与 `xi` 网格；
- 当前数值列为占位 `NaN`、`solver_status=fail`，需要被老程序真实结果替换。
